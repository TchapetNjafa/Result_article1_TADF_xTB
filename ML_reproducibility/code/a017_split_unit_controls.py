#!/usr/bin/env python3
"""A-017: does the NTO-over-Morgan rank advantage depend on the holdout unit?

Pre-registered (L3 beat 1, 2026-10-06) before any result was seen:
  Q1 fold-assignment variance: for each grouping (DOI, scaffold, molecule),
     50 random group->fold assignments (5 folds), OOF RF predictions for Morgan
     and NTO; report median / 2.5-97.5 percentile of delta_rho and the
     fraction of seeds with delta_rho > 0.
  Q2 matched control for "19/20": the A-016 repeated holdout (GroupShuffleSplit,
     20 splits, test 0.2, seed 17) rerun with scaffold and molecule grouping.
  Q3 scaffold-split delta_rho CI: scaffold-cluster bootstrap (5000) of the fixed
     scaffold GroupKFold OOF predictions (A-006 T1 folds).
Decision rule fixed in advance: if the DOI distribution in Q1 sits clearly
above scaffold/molecule (medians differ by > 0.05 and DOI fraction > 0 exceeds
the others), the holdout unit drives the contrast; otherwise the A-016 point
estimate is a single-split artefact. Retrospective, same corpus, no lab transfer.
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from _dataset import NTO, rf  # noqa: E402

N_FOLDS = 5
N_SEEDS = 50
N_BOOT = 5000
BOOT_SEED = 54321
REPEAT_SEED = 17


def oof_from_folds(features, target, fold_of):
    predictions = np.zeros(len(target))
    for fold in range(N_FOLDS):
        test = fold_of == fold
        predictions[test] = rf().fit(features[~test], target[~test]).predict(features[test])
    return predictions


def random_group_folds(groups, seed):
    # ponytail: round-robin over shuffled groups; fold sizes balanced only to max group size (4)
    unique = np.random.default_rng(seed).permutation(sorted(set(groups)))
    fold = {group: index % N_FOLDS for index, group in enumerate(unique)}
    return np.array([fold[group] for group in groups])


def rho(y, p):
    return float(spearmanr(y, p).statistic)


def summarize(values):
    values = np.asarray(values)
    return {"median": float(np.median(values)),
            "p2_5_p97_5": [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))],
            "min_max": [float(values.min()), float(values.max())],
            "fraction_positive": float(np.mean(values > 0))}


def main():
    data = pd.read_csv(ROOT / "data/a016_jmgm_doi_group_input.csv", float_precision="round_trip")
    y = data.exp_est.values
    mols = [Chem.MolFromSmiles(s) for s in data.smiles]
    assert all(m is not None for m in mols)
    morgan = np.vstack([np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048), dtype=float)
                        for m in mols])
    nto = data[NTO].values
    groupings = {
        "doi": data.primary_doi.values,
        "scaffold": np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=m) or "none" for m in mols]),
        "molecule": np.arange(len(y)).astype(str),
    }
    assert len(y) == 231 and len(set(groupings["doi"])) == 173

    # Self-check: fixed GroupKFold folds must reproduce A-006 (T1 scaffold, T2 DOI).
    prior = json.loads((ROOT / "data/a006_equivalence_ceiling.json").read_text())
    fixed = {}
    for name, key in (("scaffold", "T1_morgan_vs_nto_equivalence"), ("doi", "T2_doi_split_rho")):
        fold_of = np.zeros(len(y), dtype=int)
        for fold, (_, test) in enumerate(GroupKFold(N_FOLDS).split(morgan, y, groupings[name])):
            fold_of[test] = fold
        pm, pn = oof_from_folds(morgan, y, fold_of), oof_from_folds(nto, y, fold_of)
        assert abs(rho(y, pm) - prior[key]["morgan"]["rho"]) < 1e-10, name
        assert abs(rho(y, pn) - prior[key]["nto"]["rho"]) < 1e-10, name
        fixed[name] = (pm, pn)
    n_scaffolds = len(set(groupings["scaffold"]))

    # Q1
    q1 = {}
    for name, groups in groupings.items():
        rows = []
        for seed in range(N_SEEDS):
            fold_of = random_group_folds(groups, seed)
            pm, pn = oof_from_folds(morgan, y, fold_of), oof_from_folds(nto, y, fold_of)
            rows.append({"seed": seed, "rho_morgan": rho(y, pm), "rho_nto": rho(y, pn),
                         "mae_morgan": float(np.mean(np.abs(pm - y))),
                         "mae_nto": float(np.mean(np.abs(pn - y)))})
        frame = pd.DataFrame(rows)
        q1[name] = {"delta_rho": summarize(frame.rho_nto - frame.rho_morgan),
                    "rho_morgan_median": float(frame.rho_morgan.median()),
                    "rho_nto_median": float(frame.rho_nto.median()),
                    "delta_mae_morgan_minus_nto": summarize(frame.mae_morgan - frame.mae_nto),
                    "n_groups": len(set(groups))}
        print(name, json.dumps(q1[name]["delta_rho"]), flush=True)

    # Q2
    q2 = {}
    for name, groups in groupings.items():
        deltas = []
        for train, test in GroupShuffleSplit(n_splits=20, test_size=0.2,
                                             random_state=REPEAT_SEED).split(morgan, y, groups):
            deltas.append(rho(y[test], rf().fit(nto[train], y[train]).predict(nto[test])) -
                          rho(y[test], rf().fit(morgan[train], y[train]).predict(morgan[test])))
        q2[name] = {"n_nto_higher_of_20": int(np.sum(np.array(deltas) > 0)),
                    "median_delta_rho": float(np.median(deltas))}
    assert q2["doi"]["n_nto_higher_of_20"] == 19, "must reproduce A-016"

    # Q3
    pm, pn = fixed["scaffold"]
    groups = groupings["scaffold"]
    index = [np.flatnonzero(groups == g) for g in sorted(set(groups))]
    rng = np.random.default_rng(BOOT_SEED)
    boot = []
    for _ in range(N_BOOT):
        sample = np.concatenate([index[i] for i in rng.integers(0, len(index), len(index))])
        boot.append(rho(y[sample], pn[sample]) - rho(y[sample], pm[sample]))
    q3 = {"delta_rho": rho(y, pn) - rho(y, pm),
          "scaffold_cluster_ci_95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
          "fraction_le_zero": float(np.mean(np.array(boot) <= 0))}

    result = {"method": "A-017 holdout-unit controls; see module docstring for pre-registration",
              "scope": "retrospective, same 231-molecule corpus; no new-lab transfer",
              "n_molecules": len(y), "n_doi_groups": 173, "n_scaffolds": n_scaffolds,
              "n_seeds": N_SEEDS, "self_check_reproduces_a006": True,
              "Q1_random_fold_assignment": q1, "Q2_repeated_holdout_20": q2,
              "Q3_fixed_scaffold_split": q3}
    (ROOT / "data/a017_split_unit_controls.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"Q2": q2, "Q3": q3}, indent=2))


if __name__ == "__main__":
    main()
