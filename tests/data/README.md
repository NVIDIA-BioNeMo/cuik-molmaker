# Reaction featurizer test data

Golden reference files for `tests/python/test_reaction_features.py`. The test is
a plain positional regression check: it runs `cuik_molmaker.batch_reaction_featurizer`
on the sample reactions and compares against the committed `.xz` references.

## Files

| File | Description |
|------|-------------|
| `sample_rxns_100.csv` | 100 atom-mapped reactions sampled from RDB7 (column `rxn_smiles`, format `reactant>>product`) |
| `sample_rxns_100_<VERSION>_<MODE>_ref.xz` | 24 goldens: 4 atom featurizer versions (`V1`, `V2`, `ORGANIC`, `RIGR`) x 6 reaction modes (`REAC_DIFF`, `REAC_PROD`, `PROD_DIFF`, and their `_BALANCE` variants) |
| `gen_reaction_refs.py` | End-to-end generator: downloads RDB7, samples `sample_rxns_100.csv`, and writes the 24 `.xz` goldens |

Each `.xz` is a pickled dict with keys `reac_smiles`, `prod_smiles`, `V` (atom
features), `E` (bond features), `edge_index`, `rev_edge_index`, `batch`.

## Provenance and license

The reactions are sampled from **RDB7** (CC BY 4.0). Only the atom-mapped reaction
SMILES are used (no activation energies or other columns). Full attribution and
license text are in [`../../LICENSE/third-party.txt`](../../LICENSE/third-party.txt):

> Spiekermann, K.; Pattanaik, L.; Green, W. H. "High Accuracy Barrier Heights,
> Enthalpies, and Rate Coefficients for Chemical Reactions." *Scientific Data*
> 2022, 9, 417. Zenodo: https://doi.org/10.5281/zenodo.6618262

RDB7 is **not** redistributed in this repository; `gen_reaction_refs.py` downloads
the source `wb97xd3.csv` directly from Zenodo into a local cache (not committed).

## Regenerating

Run from a conda env in which `cuik_molmaker` is importable. A single command
reproduces everything -- it downloads `wb97xd3.csv` from Zenodo into a local
cache, samples the 100 reactions deterministically (seed=7), and writes the 24
golden `.xz` files:

```bash
python tests/data/gen_reaction_refs.py
```

The goldens are produced by the C++ `batch_reaction_featurizer`. Its output was
verified against Chemprop's `CondensedGraphOfReactionFeaturizer` on these 100
reactions across all 4 featurizer versions x 6 modes before the references were
frozen: bond features and edge indices match exactly; atom features match to
within float32 round-off (max abs diff ~6e-9, from the atomic-mass feature --
cuik computes in float32, chemprop in float64). The committed test therefore
does not import Chemprop.
