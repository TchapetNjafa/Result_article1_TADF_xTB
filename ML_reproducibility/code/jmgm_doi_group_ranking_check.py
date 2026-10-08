#!/usr/bin/env python3
"""Recheck A-006's fixed DOI-OOF ranking contrast using paper-level resampling.

This is conditional on the existing DOI GroupKFold predictions. It does not
retrain for each bootstrap sample or test transfer to an independent laboratory.
"""
import ast
import argparse
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, GroupShuffleSplit

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from _dataset import NTO, load_dataset, make_features, rf  # noqa: E402

SEED = 54321
N_BOOT = 5000
N_FOLDS = 5
N_REPEATS = 20
REPEAT_SEED = 17


def names(value):
    return ([str(name).strip() for name in ast.literal_eval(value)]
            if str(value).startswith("[") else [str(value)])


def value_number(value):
    try:
        parsed = ast.literal_eval(value)
        return float(parsed[0]) if isinstance(parsed, list) else float(parsed)
    except Exception:
        try:
            return float(re.findall(r"[-\d.]+", str(value))[0])
        except Exception:
            return np.nan


def oof(features, target, groups):
    predictions = np.zeros(len(target))
    cv = GroupKFold(n_splits=N_FOLDS)
    for train, test in cv.split(features, target, groups):
        predictions[test] = rf().fit(features[train], target[train]).predict(features[test])
    return predictions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=pathlib.Path,
                        help="Use the prepared, machine-readable supplementary dataset")
    args = parser.parse_args()
    if args.input_csv:
        prepared = pd.read_csv(args.input_csv, float_precision="round_trip")
        assert prepared[NTO + ["exp_est", "primary_doi", "smiles"]].notna().all().all()
        y = prepared.exp_est.values
        groups = prepared.primary_doi.values
        molecules = [Chem.MolFromSmiles(smiles) for smiles in prepared.smiles]
        assert all(molecule is not None for molecule in molecules)
        morgan = np.vstack([
            np.array(AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=2048),
                     dtype=float) for molecule in molecules])
        nto = prepared[NTO].values
    else:
        ds = load_dataset()
        y = ds.y
        morgan, _ = make_features(ds, "Morgan")
        nto, _ = make_features(ds, "NTO")
        exp = pd.read_csv(ROOT / "MCA_submission/tables/experimental_delta_est_extracted.csv")
        exp = exp[exp["specifier"].astype(str).str.contains("EST", case=False, na=False)].copy()
        exp["nl"] = exp["compound.names"].map(names)
        exp["v"] = exp["standard_value"].map(value_number)
        exp = exp[np.isfinite(exp.v) & exp.v.between(-0.3, 1.5)]
        molecules = set(ds.df.molecule.astype(str))
        rows = [(name, row.v, row.doi) for _, row in exp.iterrows()
                for name in row.nl if name in molecules]
        long = pd.DataFrame(rows, columns=["molecule", "v", "doi"])
        primary = long.groupby("molecule").doi.agg(lambda doi: doi.value_counts().index[0])
        mapped = ds.df.molecule.astype(str).map(primary)
        mask = mapped.notna().values
        groups = mapped[mask].values
        y = y[mask]
        morgan = morgan[mask]
        nto = nto[mask]
        prepared = ds.df.loc[mask, ["molecule", "smiles", "exp_est"] + NTO].copy()
        prepared.insert(3, "primary_doi", groups)
        prepared.to_csv(ROOT / "data/a016_jmgm_doi_group_input.csv", index=False,
                        float_format="%.17g")
    pred_m = oof(morgan, y, groups)
    pred_n = oof(nto, y, groups)
    rho_m = float(spearmanr(y, pred_m).statistic)
    rho_n = float(spearmanr(y, pred_n).statistic)
    observed = rho_n - rho_m

    prior = json.loads((ROOT / "data/a006_equivalence_ceiling.json").read_text())["T2_doi_split_rho"]
    assert len(y) == prior["n"] and len(set(groups)) == prior["n_doi_groups"]
    assert abs(rho_m - prior["morgan"]["rho"]) < 1e-10, (
        f"Morgan rho changed: {rho_m} != {prior['morgan']['rho']}")
    assert abs(rho_n - prior["nto"]["rho"]) < 1e-10, (
        f"NTO rho changed: {rho_n} != {prior['nto']['rho']}")

    indices = [np.flatnonzero(groups == group) for group in sorted(set(groups))]
    rng = np.random.default_rng(SEED)
    delta = []
    for _ in range(N_BOOT):
        chosen = rng.integers(0, len(indices), size=len(indices))
        sample = np.concatenate([indices[index] for index in chosen])
        diff = spearmanr(y[sample], pred_n[sample]).statistic - spearmanr(
            y[sample], pred_m[sample]).statistic
        if np.isfinite(diff):
            delta.append(float(diff))
    assert len(delta) == N_BOOT, "A degenerate bootstrap draw requires explicit handling"
    ci = np.percentile(delta, [2.5, 97.5])
    leave_one = []
    for excluded in indices:
        keep = np.ones(len(y), dtype=bool)
        keep[excluded] = False
        leave_one.append(float(spearmanr(y[keep], pred_n[keep]).statistic -
                               spearmanr(y[keep], pred_m[keep]).statistic))
    repeated = []
    splitter = GroupShuffleSplit(n_splits=N_REPEATS, test_size=0.2,
                                 random_state=REPEAT_SEED)
    for train, test in splitter.split(morgan, y, groups):
        m_test = rf().fit(morgan[train], y[train]).predict(morgan[test])
        n_test = rf().fit(nto[train], y[train]).predict(nto[test])
        m_rho = float(spearmanr(y[test], m_test).statistic)
        n_rho = float(spearmanr(y[test], n_test).statistic)
        assert np.isfinite(m_rho) and np.isfinite(n_rho)
        repeated.append({"n_test_molecules": len(test),
                         "n_test_papers": len(set(groups[test])),
                         "morgan_rho": m_rho, "nto_rho": n_rho,
                         "delta_rho": n_rho - m_rho})
    split_deltas = np.array([item["delta_rho"] for item in repeated])
    result = {
        "method": "DOI-group bootstrap of paired fixed GroupKFold OOF predictions",
        "scope": "conditional on this one split and fitted models; no new-lab validation",
        "seed": SEED, "n_boot": N_BOOT, "n_molecules": len(y),
        "n_primary_doi_groups": len(indices),
        "rho_morgan": rho_m, "rho_nto": rho_n,
        "delta_rho_nto_minus_morgan": observed,
        "doi_group_bootstrap_ci_95": [float(ci[0]), float(ci[1])],
        "fraction_bootstrap_delta_le_zero": float(np.mean(np.array(delta) <= 0)),
        "leave_one_doi_delta_range": [min(leave_one), max(leave_one)],
        "prior_molecule_bootstrap_ci_95": prior["delta_rho_CI"],
        "repeated_doi_holdout": {
            "n_splits": N_REPEATS, "seed": REPEAT_SEED,
            "test_size_paper_fraction": 0.2,
            "median_delta_rho": float(np.median(split_deltas)),
            "range_delta_rho": [float(np.min(split_deltas)), float(np.max(split_deltas))],
            "n_splits_nto_higher_rho": int(np.sum(split_deltas > 0)),
            "note": "overlapping paper partitions; descriptive stability, not independent samples or a confidence interval",
            "splits": repeated,
        },
        "comparison_matches_prior_points": True,
    }
    output = ROOT / "data/a016_jmgm_doi_group_ranking.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
