// SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

//! @file This file defines functions for one-hot atom and bond features,
//!       declared in one_hot.h and called from features.cpp

#include "one_hot.h"

#include <GraphMol/RingInfo.h>
#include <GraphMol/ROMol.h>
#include <RDGeneral/types.h>
#include <stdint.h>
#include <string.h>

#include <string>
#include <type_traits>

#include "features.h"
#include "float_features.h"

// Helper class to automatically generates a reverse lookup table at compile time,
// with `MAX_OUT` used as a sentinel to indicate that a value wasn't present
// in the original list.
template <size_t NUM_IN, size_t MAX_OUT> class OneHotLookup {
  size_t indices[NUM_IN];

 public:
  constexpr OneHotLookup(const size_t list[MAX_OUT]) : indices() {
    std::fill(indices, indices + NUM_IN, MAX_OUT);
    for (size_t i = 0; i < MAX_OUT; ++i) {
      indices[list[i]] = i;
    }
  }
  constexpr size_t operator[](size_t i) const { return (i < NUM_IN) ? indices[i] : MAX_OUT; }
};

// +2 to avoid negative numbers
constexpr size_t formalChargeList[] = {
  size_t(-1 + 2),
  size_t(-2 + 2),
  size_t(+1 + 2),
  size_t(+2 + 2),
  size_t(0 + 2),
};
constexpr size_t                             formalChargeCount = std::extent<decltype(formalChargeList)>::value;
constexpr OneHotLookup<5, formalChargeCount> formalChargeLookup(formalChargeList);

constexpr size_t atomicNumCount = 100;

constexpr size_t degreeCount       = 5;
constexpr size_t totalDegreeCount  = 6;
constexpr size_t valenceCount      = 7;
constexpr size_t chiralityCount    = 4;
constexpr size_t numHydrogensCount = 5;
constexpr size_t ringSizeCount     = 6;

constexpr size_t hybridizationList[] = {
  RDKit::Atom::HybridizationType::SP,
  RDKit::Atom::HybridizationType::SP2,
  RDKit::Atom::HybridizationType::SP3,
  RDKit::Atom::HybridizationType::SP3D,
  RDKit::Atom::HybridizationType::SP3D2,
};
constexpr size_t                              hybridizationCount = std::extent<decltype(hybridizationList)>::value;
constexpr OneHotLookup<8, hybridizationCount> hybridizationLookup(hybridizationList);

constexpr size_t hybridizationExpandedList[] = {
  RDKit::Atom::HybridizationType::S,
  RDKit::Atom::HybridizationType::SP,
  RDKit::Atom::HybridizationType::SP2,
  RDKit::Atom::HybridizationType::SP2D,
  RDKit::Atom::HybridizationType::SP3,
  RDKit::Atom::HybridizationType::SP3D,
  RDKit::Atom::HybridizationType::SP3D2,
};
constexpr size_t hybridizationExpandedCount = std::extent<decltype(hybridizationExpandedList)>::value;
constexpr OneHotLookup<8, hybridizationExpandedCount> hybridizationExpandedLookup(hybridizationExpandedList);

constexpr size_t hybridizationOrganicList[] = {
  RDKit::Atom::HybridizationType::S,
  RDKit::Atom::HybridizationType::SP,
  RDKit::Atom::HybridizationType::SP2,
  RDKit::Atom::HybridizationType::SP3,
};
constexpr size_t hybridizationOrganicCount = std::extent<decltype(hybridizationOrganicList)>::value;
constexpr OneHotLookup<8, hybridizationOrganicCount> hybridizationOrganicLookup(hybridizationOrganicList);

// First 4 rows of periodic table and Iodine
constexpr size_t atomicNumCommonList[] = {
  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18,
  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
  53,  // Iodine
};
constexpr size_t atomicNumCommonCount = std::extent<decltype(atomicNumCommonList)>::value;
constexpr OneHotLookup<118, atomicNumCommonCount> atomicNumCommonLookup(atomicNumCommonList);

// Organic chemistry elements only
constexpr size_t atomicNumOrganicList[] = {
  1,
  5,
  6,
  7,
  8,
  9,
  14,
  15,
  16,
  17,
  35,
  53,
};
constexpr size_t atomicNumOrganicCount = std::extent<decltype(atomicNumOrganicList)>::value;
constexpr OneHotLookup<118, atomicNumOrganicCount> atomicNumOrganicLookup(atomicNumOrganicList);

// RDKit::Bond::BondType are not in bondTypeList and RDKit::Bond::BondType::AROMATIC
// will be indicated at index 3.
constexpr size_t bondTypeList[] = {
  RDKit::Bond::BondType::SINGLE,
  RDKit::Bond::BondType::DOUBLE,
  RDKit::Bond::BondType::TRIPLE,
};
constexpr size_t                          bondTypeCount = std::extent<decltype(bondTypeList)>::value;
constexpr OneHotLookup<22, bondTypeCount> bondTypeLookup(bondTypeList);

constexpr size_t bondStereoList[] = {
  RDKit::Bond::BondStereo::STEREONONE,
  RDKit::Bond::BondStereo::STEREOANY,
  RDKit::Bond::BondStereo::STEREOZ,
  RDKit::Bond::BondStereo::STEREOE,
  RDKit::Bond::BondStereo::STEREOCIS,
  RDKit::Bond::BondStereo::STEREOTRANS,
};
constexpr size_t                           bondStereoCount = std::extent<decltype(bondStereoList)>::value;
constexpr OneHotLookup<6, bondStereoCount> bondStereoLookup(bondStereoList);

// Returns the number of values per atom, required by `feature` in `get_one_hot_atom_feature`'s
// `data` argument.
size_t get_one_hot_atom_feature_size(AtomOneHotFeature feature) {
  switch (feature) {
    case AtomOneHotFeature::ATOMIC_NUM:
      return atomicNumCount + 1;
    case AtomOneHotFeature::ATOMIC_NUM_COMMON:
      return atomicNumCommonCount + 1;
    case AtomOneHotFeature::ATOMIC_NUM_ORGANIC:
      return atomicNumOrganicCount + 1;
    case AtomOneHotFeature::DEGREE:
      return degreeCount + 1;
    case AtomOneHotFeature::TOTAL_DEGREE:
      return totalDegreeCount + 1;
    case AtomOneHotFeature::VALENCE:
      return valenceCount + 1;
    case AtomOneHotFeature::IMPLICIT_VALENCE:
      return valenceCount + 1;
    case AtomOneHotFeature::HYBRIDIZATION:
      return hybridizationCount + 1;
    case AtomOneHotFeature::HYBRIDIZATION_EXPANDED:
      return hybridizationExpandedCount + 1;
    case AtomOneHotFeature::HYBRIDIZATION_ORGANIC:
      return hybridizationOrganicCount + 1;
    case AtomOneHotFeature::CHIRALITY:
      return chiralityCount + 1;
    case AtomOneHotFeature::GROUP:
      return groupCount + 1;
    case AtomOneHotFeature::PERIOD:
      return periodCount + 1;
    case AtomOneHotFeature::FORMAL_CHARGE:
      return formalChargeCount + 1;
    case AtomOneHotFeature::NUM_HYDROGENS:
      return numHydrogensCount + 1;
    case AtomOneHotFeature::RING_SIZE:
      return ringSizeCount;

    default:
      // Missing implementation
      assert(0);
      return 0;
  }
}

// Returns the 0-based index within the atomic-number one-hot block for a given atomicNum.
// Mirrors the logic in each ATOMIC_NUM* case of get_one_hot_atom_feature.
// For unrecognized atomic numbers, returns the "other" slot (last valid index).
size_t get_atomic_num_onehot_index(uint8_t atomicNum, AtomOneHotFeature feature) {
  switch (feature) {
    case AtomOneHotFeature::ATOMIC_NUM: {
      // 1-indexed, atomicNum in [1,100] → index [0,99]; unknown → 100
      size_t idx = size_t(atomicNum);
      --idx;
      return (idx >= atomicNumCount) ? atomicNumCount : idx;
    }
    case AtomOneHotFeature::ATOMIC_NUM_COMMON:
      return atomicNumCommonLookup[size_t(atomicNum)];  // returns atomicNumCommonCount for unknowns
    case AtomOneHotFeature::ATOMIC_NUM_ORGANIC:
      return atomicNumOrganicLookup[size_t(atomicNum)];  // returns atomicNumOrganicCount for unknowns
    default:
      assert(0 && "get_atomic_num_onehot_index called with non-atomic-num feature");
      return 0;
  }
}

// Fills in a particular atom `feature`'s one-hot encoding into `data`, for all atoms.
// See the declaration in one_hot.h for more details.
template <typename T>
size_t get_one_hot_atom_feature(const GraphData& graph, T* data, AtomOneHotFeature feature, size_t stride) {
  const size_t        num_atoms          = graph.num_atoms;
  const RDKit::ROMol& mol                = *graph.mol.get();
  const size_t        feature_size       = get_one_hot_atom_feature_size(feature);
  const size_t        total_feature_size = feature_size * num_atoms;
  if (total_feature_size == 0) {
    return feature_size;
  }
  {
    T* current_data = data;
    for (size_t i = 0; i < num_atoms; ++i) {
      memset(current_data, 0, sizeof(data[0]) * feature_size);
      current_data += stride;
    }
  }
  switch (feature) {
    case AtomOneHotFeature::ATOMIC_NUM:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        size_t atomicNum = graph.atoms[atomIndex].atomicNum;
        --atomicNum;
        data[(atomicNum >= atomicNumCount) ? atomicNumCount : atomicNum] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::ATOMIC_NUM_COMMON:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        size_t atomicNum                       = graph.atoms[atomIndex].atomicNum;
        data[atomicNumCommonLookup[atomicNum]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::ATOMIC_NUM_ORGANIC:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        size_t atomicNum                        = graph.atoms[atomIndex].atomicNum;
        data[atomicNumOrganicLookup[atomicNum]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::DEGREE:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto   degree    = mol.getAtomWithIdx(atomIndex)->getDegree();
        size_t dataIndex = (degree < degreeCount) ? degree : degreeCount;
        data[dataIndex]  = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::TOTAL_DEGREE:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto   total_degree = mol.getAtomWithIdx(atomIndex)->getTotalDegree();
        size_t dataIndex    = (total_degree < totalDegreeCount) ? total_degree : totalDegreeCount;
        data[dataIndex]     = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::VALENCE:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto   valence   = mol.getAtomWithIdx(atomIndex)->getTotalValence();
        size_t dataIndex = (size_t(valence) < valenceCount) ? size_t(valence) : valenceCount;
        data[dataIndex]  = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::IMPLICIT_VALENCE:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto   valence   = mol.getAtomWithIdx(atomIndex)->getImplicitValence();
        size_t dataIndex = (size_t(valence) < valenceCount) ? size_t(valence) : valenceCount;
        data[dataIndex]  = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::HYBRIDIZATION:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto hybridization                       = mol.getAtomWithIdx(atomIndex)->getHybridization();
        data[hybridizationLookup[hybridization]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::HYBRIDIZATION_EXPANDED:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto hybridization                               = mol.getAtomWithIdx(atomIndex)->getHybridization();
        data[hybridizationExpandedLookup[hybridization]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::HYBRIDIZATION_ORGANIC:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        auto hybridization                              = mol.getAtomWithIdx(atomIndex)->getHybridization();
        data[hybridizationOrganicLookup[hybridization]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::CHIRALITY:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        std::string        chirality;
        const RDKit::Atom* atom           = mol.getAtomWithIdx(atomIndex);
        int                chiralTagValue = static_cast<int>(atom->getChiralTag());
        data[(chiralTagValue >= chiralityCount) ? chiralityCount : chiralTagValue] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::GROUP:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        size_t atomicNum = graph.atoms[atomIndex].atomicNum;
        size_t dataIndex = groupCount;
        if (atomicNum - 1 < std::extent<decltype(atomicNumToGroupTable)>::value) {
          uint8_t group = atomicNumToGroupTable[atomicNum - 1];
          // Group numbers are 1-based, but the array indices aren't.
          dataIndex     = group - 1;
        }
        data[dataIndex] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::PERIOD:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        size_t atomicNum = graph.atoms[atomIndex].atomicNum;
        size_t dataIndex = periodCount;
        if (atomicNum - 1 < std::extent<decltype(atomicNumToPeriodTable)>::value) {
          uint8_t period = atomicNumToPeriodTable[atomicNum - 1];
          // Period numbers are 1-based, but the array indices aren't.
          dataIndex      = period - 1;
        }
        data[dataIndex] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::FORMAL_CHARGE:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        int formalCharge                                   = graph.atoms[atomIndex].formalCharge;
        data[formalChargeLookup[size_t(formalCharge + 2)]] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::NUM_HYDROGENS:
      for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
        int numHydrogens = graph.atoms[atomIndex].totalNumHs;
        data[(numHydrogens >= numHydrogensCount) ? numHydrogensCount : numHydrogens] = FeatureValues<T>::one;
      }
      return feature_size;
    case AtomOneHotFeature::RING_SIZE: {
      {
        // RING_SIZE feature is not strictly one-hot; If an atom is present in multiple rings of different sizes, it
        // will be encoded as multiple ones.
        RDKit::RingInfo& ringInfo = *mol.getRingInfo();
        if (!ringInfo.isInitialized()) {
          throw std::runtime_error("RingInfo is not initialized");
        }
        for (size_t atomIndex = 0; atomIndex < num_atoms; ++atomIndex, data += stride) {
          std::vector<int> ringSizes = ringInfo.atomRingSizes(atomIndex);
          for (size_t i = 0; i < ringSizes.size(); ++i) {
            // Ring sizes from 3 - 8 marked
            if (ringSizes[i] <= 8) {
              data[ringSizes[i] - 3] = FeatureValues<T>::one;
            }
          }
        }
      }
      return feature_size;

      default:
        // Missing implementation
        assert(0);
        return feature_size;
    }
  }
}

// Explicit instantiations, so that the function can be templated
// but still be used from other cpp files.
template size_t get_one_hot_atom_feature<int16_t>(const GraphData&  graph,
                                                  int16_t*          data,
                                                  AtomOneHotFeature feature,
                                                  size_t            stride);
template size_t get_one_hot_atom_feature<float>(const GraphData&  graph,
                                                float*            data,
                                                AtomOneHotFeature feature,
                                                size_t            stride);
template size_t get_one_hot_atom_feature<double>(const GraphData&  graph,
                                                 double*           data,
                                                 AtomOneHotFeature feature,
                                                 size_t            stride);

// Returns the number of values per bond, required by `feature` in `get_one_hot_bond_feature`'s
// `data` argument.
size_t get_one_hot_bond_feature_size(BondFeature feature) {
  switch (feature) {
    case BondFeature::TYPE_ONE_HOT:
      return bondTypeCount + 1;
    case BondFeature::STEREO_ONE_HOT:
      return bondStereoCount + 1;
    default:
      break;
  }
  // Missing implementation
  assert(0);
  return 0;
}

// Fills in a particular bond `feature`'s one-hot encoding into `data`, for all bonds.
// See the declaration in one_hot.h for more details.
template <typename T>
size_t get_one_hot_bond_feature(const GraphData& graph, T* data, BondFeature feature, size_t stride) {
  const size_t num_bonds          = graph.num_bonds;
  const size_t feature_size       = get_one_hot_bond_feature_size(feature);
  const size_t total_feature_size = feature_size * num_bonds;
  if (total_feature_size == 0) {
    return 0;
  }
  {
    T* current_data = data;
    for (size_t i = 0; i < num_bonds; ++i) {
      memset(current_data, 0, sizeof(data[0]) * feature_size);
      current_data += stride;
    }
  }
  switch (feature) {
    case BondFeature::TYPE_ONE_HOT:
      for (size_t i = 0; i < num_bonds; ++i, data += stride) {
        auto type                  = graph.bonds[i].bondType;
        data[bondTypeLookup[type]] = FeatureValues<T>::one;
      }
      return feature_size;
    case BondFeature::STEREO_ONE_HOT:
      for (size_t i = 0; i < num_bonds; ++i, data += stride) {
        auto stereo                    = graph.bonds[i].stereo;
        data[bondStereoLookup[stereo]] = FeatureValues<T>::one;
      }
      return feature_size;
    default:
      // Missing implementation
      assert(0);
      return feature_size;
  }
}

// Explicit instantiations, so that the function can be templated
// but still be used from other cpp files.
template size_t get_one_hot_bond_feature<int16_t>(const GraphData& graph,
                                                  int16_t*         data,
                                                  BondFeature      feature,
                                                  size_t           stride);
template size_t get_one_hot_bond_feature<float>(const GraphData& graph,
                                                float*           data,
                                                BondFeature      feature,
                                                size_t           stride);
template size_t get_one_hot_bond_feature<double>(const GraphData& graph,
                                                 double*          data,
                                                 BondFeature      feature,
                                                 size_t           stride);
