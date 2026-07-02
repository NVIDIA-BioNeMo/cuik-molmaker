#!/usr/bin/env python
"""Regenerate the reaction test fixtures for test_reaction_features.py.

Fully self-contained and reproducible end to end: downloads the public RDB7
data file from Zenodo, samples 100 reactions deterministically, then generates
the golden reference files -- all in one command.

Run from a conda env that has ``cuik_molmaker`` importable::

    python tests/data/gen_reaction_refs.py

Pipeline (all outputs written next to this script):
  1. download ``wb97xd3.csv`` from RDB7 (Zenodo 6618262) into a local cache;
  2. build the reaction pool: forward (``rsmi>>psmi``) then reverse
     (``psmi>>rsmi``) for every row, preserving file order;
  3. shuffle with ``seed=7`` and keep the first 100 reactions that RDKit parses
     and ``cuik_molmaker`` featurizes successfully -> ``sample_rxns_100.csv``
     (column ``rxn_smiles`` = "reactant>>product");
  4. write ``sample_rxns_100_<VERSION>_<MODE>_ref.xz`` for all 4 atom featurizer
     versions x 6 reaction modes = 24 files.

The goldens are produced by the C++ ``batch_reaction_featurizer`` itself. Its
output was verified against Chemprop's ``CondensedGraphOfReactionFeaturizer``
on these 100 reactions across all 4 featurizer versions x 6 modes before the
references were frozen: bond features and edge indices match exactly; atom
features match to within float32 round-off (max abs diff ~6e-9, from the
atomic-mass feature -- cuik computes in float32, chemprop in float64). The
committed test therefore does not import Chemprop.

Data source and license (CC BY 4.0; see ../../LICENSE/third-party.txt):
  Spiekermann, K.; Pattanaik, L.; Green, W. H. "High Accuracy Barrier Heights,
  Enthalpies, and Rate Coefficients for Chemical Reactions." Scientific Data
  2022, 9, 417. Zenodo: https://doi.org/10.5281/zenodo.6618262
"""

import csv
import lzma
import os
import pickle
import random
import urllib.request

from rdkit import Chem, RDLogger

import cuik_molmaker

RDLogger.DisableLog("rdApp.*")

# Direct download link for the RDB7 wb97xd3 data file (Zenodo record 6618262).
ZENODO_URL = "https://zenodo.org/api/records/6618262/files/wb97xd3.csv/content"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
# Local cache for the downloaded source (not committed).
SRC_CSV = os.environ.get("RDB7_WB97XD3_CSV", os.path.join(OUT_DIR, "wb97xd3.csv"))
SAMPLE_CSV = os.path.join(OUT_DIR, "sample_rxns_100.csv")
N_SAMPLE = 100
SEED = 7

REACTION_MODES = [
    "REAC_DIFF",
    "REAC_PROD",
    "PROD_DIFF",
    "REAC_DIFF_BALANCE",
    "REAC_PROD_BALANCE",
    "PROD_DIFF_BALANCE",
]

# V1/V2/ORGANIC share the same bond features; RIGR uses a reduced set.
FEATURIZER_CONFIGS = {
    "V1": {
        "atom_onehot": [
            "atomic-number",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization",
        ],
        "atom_float": ["aromatic", "mass"],
        "bond": ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"],
    },
    "V2": {
        "atom_onehot": [
            "atomic-number-common",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization-expanded",
        ],
        "atom_float": ["aromatic", "mass"],
        "bond": ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"],
    },
    "ORGANIC": {
        "atom_onehot": [
            "atomic-number-organic",
            "total-degree",
            "formal-charge",
            "chirality",
            "num-hydrogens",
            "hybridization-organic",
        ],
        "atom_float": ["aromatic", "mass"],
        "bond": ["is-null", "bond-type-onehot", "conjugated", "in-ring", "stereo"],
    },
    "RIGR": {
        "atom_onehot": ["atomic-number-common", "total-degree", "num-hydrogens"],
        "atom_float": ["mass"],
        # RIGR uses reduced bond features -- only 2 vs 14 for V1/V2/ORGANIC.
        "bond": ["is-null", "in-ring"],
    },
}


def download_source():
    """Download the RDB7 wb97xd3.csv from Zenodo if it is not already cached."""
    if not os.path.exists(SRC_CSV):
        print(f"downloading RDB7 source from {ZENODO_URL}")
        urllib.request.urlretrieve(ZENODO_URL, SRC_CSV)
    return SRC_CSV


def build_reaction_pool(path):
    """Forward + reverse reactions in deterministic (file) order."""
    fwd, rev = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            rs, ps = row["rsmi"].strip(), row["psmi"].strip()
            fwd.append(f"{rs}>>{ps}")
            rev.append(f"{ps}>>{rs}")
    return fwd + rev


def valid(rxn):
    """parse+sanitize both sides; no atom at map 0; cuik featurizes OK."""
    try:
        r, p = rxn.split(">>")
        pr = Chem.SmilesParserParams()
        pr.removeHs = False
        mr = Chem.MolFromSmiles(r, pr)
        mp = Chem.MolFromSmiles(p, pr)
        if mr is None or mp is None:
            return False
        if any(a.GetAtomMapNum() == 0 for a in mr.GetAtoms()):
            return False
        if any(a.GetAtomMapNum() == 0 for a in mp.GetAtoms()):
            return False
        cfg = FEATURIZER_CONFIGS["V2"]
        ao = cuik_molmaker.atom_onehot_feature_names_to_array(cfg["atom_onehot"])
        af = cuik_molmaker.atom_float_feature_names_to_array(cfg["atom_float"])
        bf = cuik_molmaker.bond_feature_names_to_array(cfg["bond"])
        mi = cuik_molmaker.reaction_mode_to_int("REAC_DIFF")
        cuik_molmaker.batch_reaction_featurizer(
            [r], [p], ao, af, bf, True, False, False, mi
        )
        return True
    except Exception:
        return False


def sample_reactions():
    """Download RDB7, sample N_SAMPLE valid reactions (seed=SEED), write the CSV."""
    src = download_source()
    rxns = build_reaction_pool(src)
    print(f"RDB7 reaction pool (forward + reverse): {len(rxns)}")
    random.seed(SEED)
    random.shuffle(rxns)

    picked = []
    for s in rxns:
        if len(picked) >= N_SAMPLE:
            break
        if valid(s):
            picked.append(s)
    print(f"picked {len(picked)} valid RDB7 reactions")

    with open(SAMPLE_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rxn_smiles"])
        for s in picked:
            w.writerow([s])
    print(f"wrote {SAMPLE_CSV}: {len(picked)} RDB7 reactions")


def generate_goldens():
    """Featurize sample_rxns_100.csv and write the 24 golden .xz files."""
    with open(SAMPLE_CSV) as f:
        rows = [row["rxn_smiles"] for row in csv.DictReader(f)]
    # rxn_smiles column format: "reactant>>product"
    reac_smiles = [s.split(">>")[0] for s in rows]
    prod_smiles = [s.split(">>")[1] for s in rows]

    for version, cfg in FEATURIZER_CONFIGS.items():
        atom_onehot = cuik_molmaker.atom_onehot_feature_names_to_array(
            cfg["atom_onehot"]
        )
        atom_float = cuik_molmaker.atom_float_feature_names_to_array(cfg["atom_float"])
        bond_feats = cuik_molmaker.bond_feature_names_to_array(cfg["bond"])

        for mode in REACTION_MODES:
            mode_int = cuik_molmaker.reaction_mode_to_int(mode)
            V, E, edge_index, rev_edge_index, batch = (
                cuik_molmaker.batch_reaction_featurizer(
                    reac_smiles,
                    prod_smiles,
                    atom_onehot,
                    atom_float,
                    bond_feats,
                    True,
                    False,
                    False,
                    mode_int,
                )
            )
            ref = {
                "reac_smiles": reac_smiles,
                "prod_smiles": prod_smiles,
                "V": V,
                "E": E,
                "edge_index": edge_index,
                "rev_edge_index": rev_edge_index,
                "batch": batch,
            }
            out_path = os.path.join(OUT_DIR, f"sample_rxns_100_{version}_{mode}_ref.xz")
            with lzma.open(out_path, "wb") as f:
                pickle.dump(ref, f)
            print(f"Wrote {out_path}  V={V.shape}  E={E.shape}")


def main():
    sample_reactions()
    generate_goldens()


if __name__ == "__main__":
    main()
