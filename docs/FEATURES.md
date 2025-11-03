# Atom One-Hot Features

The following features are available for atom featurization in one-hot encoding format and are available in the `atom_onehot_feature_names_to_array` function. They are also listed in the `list_all_atom_onehot_features` function.

| Feature | Description |
|---------|------------|
| `atomic-number` | All atomic numbers from 1 - 100 |
| `atomic-number-common` | First 4 rows of periodic table and Iodine |
| `atomic-number-organic` | Organic chemistry elements only |
| `degree` | Number of explicit neighboring atoms |
| `total-degree` | Number of neighboring atoms including hydrogens |
| `valence` | Total valence of the atom |
| `implicit-valence` | Implicit valence of the atom |
| `hybridization` | Hybridization type of the atom (SP, SP2, SP3, SP3D, SP3D2) |
| `hybridization-organic` | Hybridization type of the atom for organic elements (S, SP, SP2, SP3) |
| `chirality` | Chirality of the atom  |
| `group` | Periodic table group of the atom |
| `period` | Periodic table period of the atom |
| `formal-charge` | Formal charge on atom |
| `num-hydrogens` | Total number of hydrogens (explicit and implicit) on an atom |
| `ring-size` | Size of the ring the atom is present in (ring size 3 - 8 are multi-hot encoded) |


# Atom Float Features

The following features are available for atom featurization in floating-point format and are available in the `atom_float_feature_names_to_array` function. They are also listed in the `list_all_atom_float_features` function.

| Feature | Description |
|---------|------------|
| `atomic-number` | Atomic number of the atom |
| `mass` | Atomic mass of the atom |
| `valence` | Total valence of the atom |
| `implicit-valence` | Implicit valence of the atom |
| `hybridization` | Numerical representation of hybridization state |
| `chirality` | Numerical encoding of atom chirality |
| `aromatic` | Whether the atom is aromatic (1.0) or not (0.0) |
| `in-ring` | Whether the atom is in a ring (1.0) or not (0.0) |
| `min-ring` | Size of the smallest ring containing this atom |
| `max-ring` | Size of the largest ring containing this atom |
| `num-ring` | Number of rings the atom is part of |
| `degree` | Number of bonded neighbors |
| `radical-electron` | Number of radical electrons |
| `formal-charge` | Formal charge on the atom |
| `group` | Periodic table group number |
| `period` | Periodic table period number |
| `single-bond` | Number of single bonds |
| `aromatic-bond` | Number of aromatic bonds |
| `double-bond` | Number of double bonds |
| `triple-bond` | Number of triple bonds |
| `is-carbon` | Whether the atom is carbon (1.0) or not (0.0) |
| `hydrogen-bond-donor` | Whether the atom can donate hydrogen bonds (1.0) or not (0.0) determined by SMARTS match |
| `hydrogen-bond-acceptor` | Whether the atom can accept hydrogen bonds (1.0) or not (0.0) determined by SMARTS match |
| `acidic` | Whether the atom is acidic (1.0) or not (0.0) determined by SMARTS match |
| `basic` | Whether the atom is basic (1.0) or not (0.0) determined by SMARTS match |

# Bond Features

The following features are available for bond featurization and are available in the `bond_feature_names_to_array` function. They are also listed in the `list_all_bond_features` function.

| Feature | Description |
|---------|------------|
| `is-null` | Binary indicator: 1 if a bond is present, 0 otherwise |
| `bond-type-float` | Numerical representation of bond type (e.g., 2.0 for double bond, 1.5 for aromatic) |
| `bond-type-onehot` | One-hot encoding of selected bond types  |
| `in-ring` | Binary indicator: 1.0 if the bond is in at least one ring, 0.0 otherwise |
| `conjugated` | Binary indicator: 1.0 if the bond is conjugated, 0.0 otherwise |
| `stereo` | One-hot encoding of bond stereo configurations |
| `conformer-bond-length` | Length of the bond from a conformer (either first or computed) |
| `estimated-bond-length` | Length of the bond estimated using a fast heuristic |
