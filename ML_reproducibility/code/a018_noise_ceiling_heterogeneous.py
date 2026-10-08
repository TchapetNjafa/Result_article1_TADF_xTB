#!/usr/bin/env python3
"""A-018: R2_max with per-molecule label noise (replaces A-015's uniform 0.049 eV).

Pre-registered (L1, 2026-10-06), from L3 pass-3 H4: A-015 applied the RMS SEM of
the median (computed on the 102 multi-report molecules) to all 231 molecules.
Single-report molecules have no replicate, so their label noise is one report's
noise. Here:
  noise_var_i = (SD_i / sqrt(n_i))**2  for multi-report molecules (as A-015/reconcile)
  noise_var_i = sigma_report**2        for single-report molecules
  sigma_report in {0.0856 pooled within-molecule SD (primary), 0.0905, 0.072}
  R2_max = 1 - mean(noise_var) / var(y);  molecule bootstrap 5000, seed 54321.
Assumptions kept from A-015 (stated, not tested): additive noise independent of
structure; reports independent (58% zero-spread pairs suggest re-citation, which
would make the multi-report SEM too small). SD/sqrt(n) is the SEM of the mean;
the median's SEM is larger, so multi-report noise is a lower bound.
"""
import ast
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from _dataset import load_dataset  # noqa: E402

SIGMA_REPORT = {"pooled_within_sd_0.0856": 0.0856, "single_report_about_median_0.0905": 0.0905,
                "molecule_weighted_0.072": 0.072}
N_BOOT, SEED = 5000, 54321
R2_OBSERVED = 0.2567  # A-015 / A-006 scaffold-OOF, unchanged


def value(v):
    try:
        p = ast.literal_eval(v)
        return float(p[0]) if isinstance(p, list) else float(p)
    except Exception:
        try:
            return float(re.findall(r"[-\d.]+", str(v))[0])
        except Exception:
            return np.nan


def main():
    ds = load_dataset()
    names = ds.df.molecule.astype(str).values
    exp = pd.read_csv(ROOT / "MCA_submission/tables/experimental_delta_est_extracted.csv")
    exp = exp[exp["specifier"].astype(str).str.contains("EST", case=False, na=False)].copy()
    exp["nl"] = exp["compound.names"].map(
        lambda x: [str(n).strip() for n in ast.literal_eval(x)] if str(x).startswith("[") else [str(x)])
    exp["v"] = exp["standard_value"].map(value)
    exp = exp[np.isfinite(exp.v) & exp.v.between(-0.3, 1.5)]
    reports = {}
    for _, r in exp.iterrows():
        for n in r["nl"]:
            reports.setdefault(n, []).append(r.v)

    n_rep = np.array([len(reports.get(n, [])) for n in names])
    sem2 = np.array([np.var(reports[n], ddof=1) / len(reports[n]) if len(reports.get(n, [])) > 1
                     else np.nan for n in names])
    multi = n_rep > 1
    assert multi.sum() == 102, multi.sum()
    assert abs(np.sqrt(np.mean(sem2[multi])) - 0.049) < 5e-4, "must reproduce A-015 SEM"
    y = ds.y
    assert len(y) == 231
    var_y = float(np.var(y))

    rng = np.random.default_rng(SEED)
    boots = [rng.integers(0, len(y), len(y)) for _ in range(N_BOOT)]
    out = {}
    for label, sigma in SIGMA_REPORT.items():
        noise = np.where(multi, sem2, sigma ** 2)
        r2max = 1 - noise.mean() / var_y
        ci = np.percentile([1 - noise[b].mean() / np.var(y[b]) for b in boots], [2.5, 97.5])
        out[label] = {"rms_noise_eV": float(np.sqrt(noise.mean())), "r2_max": float(r2max),
                      "r2_max_ci95": [float(ci[0]), float(ci[1])],
                      "pct_ceiling_captured": float(100 * R2_OBSERVED / r2max)}
    result = {"method": "A-018 heterogeneous per-molecule noise ceiling; see docstring",
              "n_molecules": 231, "n_multi_report": int(multi.sum()),
              "n_single_report": int((~multi).sum()), "n_zero_reports_matched": int((n_rep == 0).sum()),
              "sigma_y_eV": float(np.sqrt(var_y)), "r2_observed": R2_OBSERVED,
              "a015_uniform_sem_r2_max": 0.9139, "variants": out}
    (ROOT / "data/a018_noise_ceiling_heterogeneous.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
