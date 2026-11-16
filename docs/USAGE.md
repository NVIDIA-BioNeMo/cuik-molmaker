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

#### Controlling atom order (optional)
Both `mol_featurizer` and `batch_mol_featurizer` accept an optional `ordered` argument, which
defaults to `True`. When it is true and the SMILES string carries atom map numbers that form a
complete ordering of the atoms — numbered either `0..n-1` or `1..n` — the atoms are reordered
accordingly, so that features are emitted in the order you asked for:

```python
# Atoms come out in the order O, C rather than the order written in the SMILES.
cuik_molmaker.mol_featurizer("[CH3:2][OH:1]", ..., ordered=True)
```

Atom maps that do not form such an ordering (partial mapping, duplicates, gaps, or a range not
starting at 0 or 1) are ignored, leaving the atoms in the order the SMILES parser produced them.
Pass `ordered=False` to ignore atom map numbers entirely. Either way the map numbers are stripped
before featurization and never appear in the features.

#### Keeping explicit hydrogens and ignoring stereochemistry (optional)
`mol_featurizer` and `batch_mol_featurizer` also accept `keep_h` and `ignore_stereo`,
both defaulting to `False`:

```python
# keep_h retains explicit hydrogens already written in the SMILES; it does not add
# any that aren't there (use explicit_H for that).
cuik_molmaker.mol_featurizer("[H]C([H])([H])O", ..., keep_h=True)  # 5 atoms, not 2

# ignore_stereo clears R/S and cis/trans stereochemistry from the molecule
# (RDKit::MolOps::removeStereochemistry) before featurizing.
cuik_molmaker.mol_featurizer("N[C@@H](C)C(=O)O", ..., ignore_stereo=True)
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

### For a batch of reactions (Condensed Graph of Reaction)
```python
# Reactant and product SMILES must be atom-mapped: corresponding atoms carry the
# same `:N` map number on both sides (see the Daylight SMILES/SMIRKS spec).
# Providing a correct, unique mapping is the caller's responsibility; duplicate
# map numbers on a side are not validated and give an arbitrary correspondence.
reac_smiles_list = ["[CH3:1][Br:2].[OH-:3]", "[CH3:1][CH2:2][Cl:3].[F-:4]"]
prod_smiles_list = ["[CH3:1][OH:3].[Br-:2]", "[CH3:1][CH2:2][F:4].[Cl-:3]"]

# Reaction mode controls how reactant/product features are combined into the CGR.
# Available modes: REAC_DIFF (default), REAC_PROD, PROD_DIFF and their _BALANCE variants.
mode = cuik_molmaker.reaction_mode_to_int("REAC_DIFF")

# keep_h keeps hydrogens that are explicitly written in the input SMILES; it does
# not add any, it only keeps the ones already specified (e.g. a mapped `[H:3]`).
# add_h adds new hydrogens via RDKit AddHs; these are unmapped, so for a CGR they
# become reactant- or product-only phantom atoms and are usually not what you want.
keep_h = True
add_h = False

# ignore_stereo clears R/S and cis/trans stereochemistry from both the reactant and
# product before featurizing (RDKit::MolOps::removeStereochemistry). Defaults to False.
ignore_stereo = False

rxn_features = cuik_molmaker.batch_reaction_featurizer(
    reac_smiles_list, prod_smiles_list,
    atom_onehot_feature_array, atom_float_feature_array, bond_feature_array,
    keep_h, add_h, offset_carbon, mode, ignore_stereo=ignore_stereo)

# CGR atom features from all reactions, concatenated along dimension 0.
# The feature dimension is doubled relative to a single molecule (reactant + product).
print(rxn_features[0].shape) # (total_num_atoms, cgr_atom_feature_dim)

# CGR bond features from all reactions, concatenated along dimension 0
print(rxn_features[1].shape) # (total_num_directed_edges, cgr_bond_feature_dim)

# Edge indices in COO format (row 0 = source atom, row 1 = destination atom)
print(rxn_features[2].shape) # (2, total_num_directed_edges)

# Reverse edge index: for each directed edge, the index of its opposite-direction edge
print(rxn_features[3].shape) # (total_num_directed_edges,)

# Associate node index: indicates the reaction idx each node belongs to
print(rxn_features[4].shape) # (total_num_atoms,)
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


