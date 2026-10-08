#!/usr/bin/env python3
"""A-020: three controls requested by L3 beat 2 (pre-registered before running).

Q1 (pass-2 H-B) RDKit control: the 217 RDKit 2D descriptors through the A-017/A-019
   protocol (molecule and doi holdout, 50 seeds). Rule: FINGERPRINT-LIKE if median
   multi-paper rho drop (random - doi) >= 0.15 AND single-paper rho (doi) < 0.17;
   NTO-LIKE if drop <= 0.08 AND single-paper rho >= 0.17; else INTERMEDIATE.
Q2 (pass-1 H-new-1) within-paper ranking under doi holdout: Spearman between
   paper-centred label and paper-centred prediction on the 103 multi-paper
   molecules, per representation. Rule: "ranks within papers" only if median
   > 0.2; otherwise the across-unseen-paper wording is required.
Q3 (pass-1 H-new-2) single-paper contrast interval: OOF predictions averaged
   over the 50 doi-holdout seeds; paired molecule bootstrap (5000, seed 54321)
   of rho_NTO - rho_Morgan on the 128 single-paper molecules. Rule: keep the
   contrast in the abstract only if the 95% interval excludes 0.
Retrospective, same corpus.
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
from _dataset import NTO, load_dataset, make_features  # noqa: E402
from a017_split_unit_controls import N_SEEDS, oof_from_folds, random_group_folds, rho  # noqa: E402

N_BOOT, BOOT_SEED = 5000, 54321


def centred(values, groups):
    return values - pd.Series(values).groupby(groups).transform("mean").values


def main():
    data = pd.read_csv(ROOT / "data/a016_jmgm_doi_group_input.csv", float_precision="round_trip")
    ds = load_dataset()
    assert (ds.df.smiles.values == data.smiles.values).all() and np.array_equal(ds.y, data.exp_est.values)
    y, doi = data.exp_est.values, data.primary_doi.values
    mols = [Chem.MolFromSmiles(s) for s in data.smiles]
    reps = {"morgan": np.vstack([np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048),
                                          dtype=float) for m in mols]),
            "nto": data[NTO].values,
            "rdkit": make_features(ds, "RDKit")[0]}
    size = pd.Series(doi).map(pd.Series(doi).value_counts()).values
    multi = size > 1
    assert multi.sum() == 103
    sib = np.array([y[(doi == doi[i]) & (np.arange(len(y)) != i)].mean() if multi[i] else np.nan
                    for i in range(len(y))])
    yc = centred(y[multi], doi[multi])
    units = {"molecule": np.arange(len(y)).astype(str), "doi": doi}

    rows, doi_preds = [], {r: [] for r in reps}
    for seed in range(N_SEEDS):
        for unit, groups in units.items():
            folds = random_group_folds(groups, seed)
            for rep, X in reps.items():
                p = oof_from_folds(X, y, folds)
                if unit == "doi":
                    doi_preds[rep].append(p)
                rows.append({"seed": seed, "unit": unit, "rep": rep, "rho_all": rho(y, p),
                             "rho_multi": rho(y[multi], p[multi]), "rho_single": rho(y[~multi], p[~multi]),
                             "sibling_pull": rho(sib[multi], p[multi]),
                             "within_paper_rho": rho(yc, centred(p[multi], doi[multi]))})
    frame = pd.DataFrame(rows)
    med = frame.groupby(["rep", "unit"]).median(numeric_only=True).drop(columns="seed")

    # Q1
    w = frame[frame.rep == "rdkit"].pivot(index="seed", columns="unit")
    drop = (w["rho_multi"]["molecule"] - w["rho_multi"]["doi"]).values
    d_multi, single = float(np.median(drop)), float(med.loc[("rdkit", "doi"), "rho_single"])
    q1 = ("FINGERPRINT-LIKE" if d_multi >= 0.15 and single < 0.17 else
          "NTO-LIKE" if d_multi <= 0.08 and single >= 0.17 else "INTERMEDIATE")

    # Q2
    within = {r: float(med.loc[(r, "doi"), "within_paper_rho"]) for r in reps}

    # Q3
    avg = {r: np.mean(doi_preds[r], axis=0)[~multi] for r in reps}
    ys = y[~multi]
    rng = np.random.default_rng(BOOT_SEED)
    boot = []
    for _ in range(N_BOOT):
        b = rng.integers(0, len(ys), len(ys))
        boot.append(rho(ys[b], avg["nto"][b]) - rho(ys[b], avg["morgan"][b]))
    ci = np.percentile(boot, [2.5, 97.5])
    q3 = {"rho_morgan_avgpred": rho(ys, avg["morgan"]), "rho_nto_avgpred": rho(ys, avg["nto"]),
          "rho_rdkit_avgpred": rho(ys, avg["rdkit"]),
          "delta_nto_minus_morgan": rho(ys, avg["nto"]) - rho(ys, avg["morgan"]),
          "ci95": [float(ci[0]), float(ci[1])], "fraction_le_zero": float(np.mean(np.array(boot) <= 0)),
          "keep_in_abstract": bool(ci[0] > 0)}

    result = {"method": "A-020 controls; pre-registration in docstring", "n_seeds": N_SEEDS,
              "medians": {f"{r}|{u}": v for (r, u), v in med.round(6).to_dict(orient="index").items()},
              "Q1_rdkit": {"multi_drop_median": d_multi,
                           "multi_drop_p2_5_p97_5": [float(np.percentile(drop, 2.5)),
                                                     float(np.percentile(drop, 97.5))],
                           "single_paper_rho_doi": single, "verdict": q1},
              "Q2_within_paper_rho_doi": within,
              "Q2_ranks_within_papers": {r: v > 0.2 for r, v in within.items()},
              "Q3_single_paper_contrast": q3}
    (ROOT / "data/a020_controls.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
