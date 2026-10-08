#!/usr/bin/env python3
"""A-019: do Morgan fingerprints exploit a publication/series shortcut?

Pre-registered (L1, 2026-10-06, after A-017, before running this script):
Hypothesis H: under random (molecule) holdout, Morgan predictions borrow the
label level of same-paper siblings present in training; under paper holdout
that shortcut is removed, which explains Morgan's rho drop (A-017).
Same 50 random fold assignments per unit as A-017 (molecule, doi).
  T1 subset drop: rho(random) - rho(doi) on molecules from multi-molecule papers
     (103) vs single-molecule papers (128), per representation.
  T2 sibling pull: under random holdout, Spearman between the OOF prediction and
     the leave-self-out mean label of same-paper siblings (multi-molecule papers
     only), minus the same quantity under paper holdout; per representation.
Decision rule (fixed): H SUPPORTED if (a) Morgan's T1 drop in multi-paper
molecules is >= 2x its drop in single-paper molecules AND exceeds NTO's
multi-paper drop, and (b) Morgan's T2 sibling pull (random minus doi) exceeds
NTO's. Both medians over seeds. Otherwise NOT SUPPORTED (report which part failed).
Retrospective, same corpus; a supported H is an association, not proof of cause.
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from _dataset import NTO  # noqa: E402
from a017_split_unit_controls import N_SEEDS, oof_from_folds, random_group_folds, rho  # noqa: E402


def main():
    data = pd.read_csv(ROOT / "data/a016_jmgm_doi_group_input.csv", float_precision="round_trip")
    y = data.exp_est.values
    doi = data.primary_doi.values
    mols = [Chem.MolFromSmiles(s) for s in data.smiles]
    reps = {"morgan": np.vstack([np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048),
                                          dtype=float) for m in mols]),
            "nto": data[NTO].values}
    size = pd.Series(doi).map(pd.Series(doi).value_counts()).values
    multi = size > 1
    assert multi.sum() == 103 and (~multi).sum() == 128
    sib_mean = np.array([y[(doi == doi[i]) & (np.arange(len(y)) != i)].mean() if multi[i] else np.nan
                         for i in range(len(y))])
    units = {"molecule": np.arange(len(y)).astype(str), "doi": doi}

    rows = []
    for seed in range(N_SEEDS):
        for unit, groups in units.items():
            fold_of = random_group_folds(groups, seed)
            for rep, X in reps.items():
                p = oof_from_folds(X, y, fold_of)
                rows.append({"seed": seed, "unit": unit, "rep": rep,
                             "rho_multi": rho(y[multi], p[multi]),
                             "rho_single": rho(y[~multi], p[~multi]),
                             "sibling_pull": rho(sib_mean[multi], p[multi]),
                             "mae_multi": float(np.mean(np.abs(p[multi] - y[multi]))),
                             "mae_single": float(np.mean(np.abs(p[~multi] - y[~multi])))})
    frame = pd.DataFrame(rows)
    wide = frame.pivot_table(index=["seed", "rep"], columns="unit").reset_index()

    out = {}
    for rep in reps:
        w = wide[wide.rep == rep]
        diff = {k: (w[(k, "molecule")] - w[(k, "doi")]).values
                for k in ("rho_multi", "rho_single", "sibling_pull", "mae_multi", "mae_single")}
        out[rep] = {k + "_random_minus_doi": {"median": float(np.median(v)),
                                              "p2_5_p97_5": [float(np.percentile(v, 2.5)),
                                                             float(np.percentile(v, 97.5))]}
                    for k, v in diff.items()}
        out[rep]["medians_by_unit"] = {f"{k}_{u}": float(w[(k, u)].median())
                                       for k in ("rho_multi", "rho_single", "sibling_pull")
                                       for u in units}
    med = lambda rep, k: out[rep][k + "_random_minus_doi"]["median"]  # noqa: E731
    a = (med("morgan", "rho_multi") >= 2 * med("morgan", "rho_single")
         and med("morgan", "rho_multi") > med("nto", "rho_multi"))
    b = med("morgan", "sibling_pull") > med("nto", "sibling_pull")
    result = {"method": "A-019 series-shortcut test; pre-registration in docstring",
              "n_multi_paper_molecules": 103, "n_single_paper_molecules": 128, "n_seeds": N_SEEDS,
              "results": out, "criterion_a_subset_drop": bool(a), "criterion_b_sibling_pull": bool(b),
              "verdict": "SUPPORTED" if a and b else "NOT SUPPORTED"}
    (ROOT / "data/a019_series_shortcut_test.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
