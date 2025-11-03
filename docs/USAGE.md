# How to use cuik-molmaker
Ensure that you have installed `cuik-molmaker` from [NVIDIA PyPI](https://pypi.nvidia.com) or built it from source. See [README.md](../README.md) for more details.

## Generate atom and bond features
### For a single molecule
#### Form arrays of required features
```python
import cuik_molmaker

# List all available atom onehot features
print(cuik_molmaker.list_all_atom_onehot_features())

# List all available atom float features
print(cuik_molmaker.list_all_atom_float_features())

# List all available bond features
print(cuik_molmaker.list_all_bond_features())

atom_onehot_feature_array = cuik_molmaker.atom_onehot_feature_names_to_array(['atomic-number', 'total-degree', 'formal-charge'])
atom_float_feature_array = cuik_molmaker.atom_float_feature_names_to_array(['mass', 'aromatic'])
bond_feature_array = cuik_molmaker.bond_feature_names_to_array(['bond-type-onehot', 'conjugated'])
```
If any of the features are not needed, pass an empty array by setting the array to `np.array([])`.

#### Set parameters for generation
```python
smiles = "CC(=O)O"
# Include explicit hydrogens in molecular graph
explicit_h = False

# For some float features, this substracts the corresponding feature value for carbon atom
offset_carbon = False
# If true, bond features will be duplicated. This is useful for GNNs that use directed edges.
# In small molecule cases, (most) bonds are undirected and the forward/backward edge features are the same.
duplicate_edges = True

# Adds an edge connecting an atom to itself. This is useful for GNNs that use self-loops.
add_self_loop = False
```

#### Generate atom and bond features
```python
all_features =cuik_molmaker.mol_featurizer(smiles, atom_onehot_feature_array, atom_float_feature_array, bond_feature_array, explicit_h, offset_carbon, duplicate_edges, add_self_loop)

# This returns a list of NumPy arrays.
# First index contains atom features as a NumPy array
# Atom features are concatencated from all one-hot features followed by all float features
print(all_features[0].shape) # (num_atoms, atom_feature_dim)

# Second index contains bond features as a NumPy array
print(all_features[1].shape) # (2*num_bonds, bond_feature_dim)

# Third index contains edge indices in COO format as a NumPy array
print(all_features[2].shape) # (2, 2*num_bonds)
```

### For a batch of molecules
```python
smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", # aspirin
               "CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1", # diphenhydramine
]
batch_features = cuik_molmaker.batch_mol_featurizer(smiles_list, atom_onehot_feature_array, atom_float_feature_array, bond_feature_array, explicit_h, offset_carbon, duplicate_edges, add_self_loop)

# Atom features from all molecules are concatenated along dimension 0
print(batch_features[0].shape) # (total_num_atoms, atom_feature_dim)

# Bond features from all molecules are concatenated along dimension 0
print(batch_features[1].shape) # (2*total_num_bonds, bond_feature_dim)

# Edge indices of different molecules are concatenated along dimension 1
print(batch_features[2].shape) # (2, 2*total_num_bonds)

# Reverse edge index: Reverse of the edge index
print(batch_features[3].shape) # (2*total_num_bonds,)

# Associate node index: Indicates the molecule idx each node belongs to
print(batch_features[4].shape) # (total_num_atoms,)
```


## Generate molecule features
#### Generate RDKit 2D descriptors for a list of molecules 
```python
from cuik_molmaker.mol_features import MoleculeFeaturizer

featurizer = MoleculeFeaturizer(molecular_descriptor_type="rdkit2D", rdkit2D_normalization_type="fast")

smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", # aspirin
               "CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1", # diphenhydramine
]
rdkit2D_descriptors = featurizer.featurize(smiles_list)

# Print the shape of the descriptors
# num_rdkit2D_descriptors depends on the version of RDKit used. It is 217 for RDKit 2025.03.2
print(rdkit2D_descriptors.shape) # (num_molecules, num_rdkit2D_descriptors)

```

#### Generate RDKit 2D descriptors and normalize them
Normalization is required for use in GNNs. Three types of normalization are supported:
- `descriptastorus`: Normalization parameters are borrowed from [Descriptastorus](https://github.com/bp-kelley/descriptastorus) package
- `best`: Best fitting normalization functions for a sample of molecules from [ChEMBL](https://www.ebi.ac.uk/chembl/)
- `fast`: Fast normalization functions for a sample of molecules from [ChEMBL](https://www.ebi.ac.uk/chembl/). These normalization functions deviate from the `best` ones by a small present tolerance value.
```python

featurizer = MoleculeFeaturizer(molecular_descriptor_type="rdkit2D", rdkit2D_normalization_type="fast")

smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", # aspirin
               "CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1", # diphenhydramine
]
rdkit2D_descriptors = featurizer.featurize(smiles_list)

# Verify normalization
print(rdkit2D_descriptors.min(), rdkit2D_descriptors.max()) # (0.0, 1.0)
```


