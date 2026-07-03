# Reaction featurizer test data

Golden reference files for `tests/python/test_reaction_features.py`. The test is
a plain positional regression check: it runs `cuik_molmaker.batch_reaction_featurizer`
on the sample reactions and compares against the committed `.xz` references.

The goldens are generated **from Chemprop, not from cuik-molmaker**, so the test
checks cuik-molmaker against an independent reference rather than against itself.
The recipe below is fully self-contained and reproducible.

## Files

| File | Description |
|------|-------------|
| `sample_rxns_100.csv` | 100 atom-mapped reactions sampled from RDB7 (column `rxn_smiles`, format `reactant>>product`). This is the frozen input to the recipe below. |
| `sample_rxns_100_<VERSION>_<MODE>_ref.xz` | 24 goldens: 4 atom featurizer versions (`V1`, `V2`, `ORGANIC`, `RIGR`) x 6 reaction modes (`REAC_DIFF`, `REAC_PROD`, `PROD_DIFF`, and their `_BALANCE` variants). |

Each `.xz` is a pickled dict with keys `reac_smiles`, `prod_smiles`, `V` (atom
features), `E` (bond features), `edge_index`, `rev_edge_index`, `batch`.

## Provenance and license

The reactions are sampled from **RDB7** (CC BY 4.0). Only the atom-mapped reaction
SMILES are used (no activation energies or other columns). Full attribution and
license text are in [`../../LICENSE/third-party.txt`](../../LICENSE/third-party.txt):

> Spiekermann, K.; Pattanaik, L.; Green, W. H. "High Accuracy Barrier Heights,
> Enthalpies, and Rate Coefficients for Chemical Reactions." *Scientific Data*
> 2022, 9, 417. Zenodo: https://doi.org/10.5281/zenodo.6618262

RDB7 is **not** redistributed in this repository. `sample_rxns_100.csv` was built
by downloading the source `wb97xd3.csv` from Zenodo, forming the forward
(`rsmi>>psmi`) and reverse (`psmi>>rsmi`) reaction for every row in file order,
shuffling with `seed=7`, and keeping the first 100 reactions that RDKit parses
with all atoms mapped. That file is committed and frozen; regenerating the
goldens below does not require re-sampling.

## Regenerating the goldens

The goldens are produced by Chemprop's `CondensedGraphOfReactionFeaturizer` and
`BatchMolGraph` collation. cuik-molmaker is never imported, so the references are
independent of the code under test. Pin Chemprop to a released version for
reproducibility:

```bash
conda create -n chemprop_2_2_4 python=3.11
conda activate chemprop_2_2_4
pip install chemprop==2.2.4
python gen_goldens_chemprop.py   # the script below, run from this directory
```

The committed goldens were generated with:
`python 3.11.15`, `chemprop 2.2.4`, `rdkit 2026.03.3`, `numpy 2.4.6`, `torch 2.12.1`.

`gen_goldens_chemprop.py`:

```python
#!/usr/bin/env python
"""Regenerate the reaction golden references from Chemprop.

Independent of cuik-molmaker: golden values are produced by Chemprop's own
CondensedGraphOfReactionFeaturizer and BatchMolGraph collation. The committed
test compares cuik-molmaker's batch_reaction_featurizer output against them.
"""
import csv
import lzma
import os
import pickle

from chemprop.data.collate import BatchMolGraph
from chemprop.featurizers.atom import MultiHotAtomFeaturizer, RIGRAtomFeaturizer
from chemprop.featurizers.bond import MultiHotBondFeaturizer, RIGRBondFeaturizer
from chemprop.featurizers.molgraph import CondensedGraphOfReactionFeaturizer
from chemprop.utils import make_mol

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_CSV = os.path.join(DATA_DIR, "sample_rxns_100.csv")

MODES = [
    "REAC_DIFF",
    "REAC_PROD",
    "PROD_DIFF",
    "REAC_DIFF_BALANCE",
    "REAC_PROD_BALANCE",
    "PROD_DIFF_BALANCE",
]


def featurizers(version):
    """Map an atom-featurizer version to Chemprop's (atom, bond) featurizers."""
    if version == "V1":
        return MultiHotAtomFeaturizer.v1(), MultiHotBondFeaturizer()
    if version == "V2":
        return MultiHotAtomFeaturizer.v2(), MultiHotBondFeaturizer()
    if version == "ORGANIC":
        return MultiHotAtomFeaturizer.organic(), MultiHotBondFeaturizer()
    if version == "RIGR":
        return RIGRAtomFeaturizer(), RIGRBondFeaturizer()
    raise ValueError(f"unknown version: {version}")


def main():
    with open(SAMPLE_CSV) as f:
        rows = [r["rxn_smiles"] for r in csv.DictReader(f)]
    reac_smiles = [s.split(">>")[0] for s in rows]
    prod_smiles = [s.split(">>")[1] for s in rows]

    # keep_h=True: RDB7 reactions carry mapped explicit hydrogens.
    reac_mols = [make_mol(s, keep_h=True, add_h=False) for s in reac_smiles]
    prod_mols = [make_mol(s, keep_h=True, add_h=False) for s in prod_smiles]

    for version in ["V1", "V2", "ORGANIC", "RIGR"]:
        atom_f, bond_f = featurizers(version)
        for mode in MODES:
            cgr = CondensedGraphOfReactionFeaturizer(
                atom_featurizer=atom_f, bond_featurizer=bond_f, mode_=mode
            )
            bmg = BatchMolGraph([cgr((r, p)) for r, p in zip(reac_mols, prod_mols)])
            ref = {
                "reac_smiles": reac_smiles,
                "prod_smiles": prod_smiles,
                "V": bmg.V.numpy(),
                "E": bmg.E.numpy(),
                "edge_index": bmg.edge_index.numpy(),
                "rev_edge_index": bmg.rev_edge_index.numpy(),
                "batch": bmg.batch.numpy(),
            }
            out = os.path.join(DATA_DIR, f"sample_rxns_100_{version}_{mode}_ref.xz")
            with lzma.open(out, "wb") as f:
                pickle.dump(ref, f)
            print(
                f"wrote {os.path.basename(out)}  V={ref['V'].shape}  E={ref['E'].shape}"
            )


if __name__ == "__main__":
    main()
```

On the pinned versions above, `cuik_molmaker.batch_reaction_featurizer` reproduces
these Chemprop references exactly (zero difference across all 24 files); the
regression test therefore uses an exact positional comparison and does not import
Chemprop.
