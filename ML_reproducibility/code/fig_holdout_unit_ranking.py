"""Figure: rank correlation by holdout unit (A-017), paper multiplicity (A-019 medians, as stored in A-020), RDKit control (A-020).
Reads data/a017_split_unit_controls.json and
data/a020_controls.json only (RDKit has no scaffold-holdout run: plotted as n.t.)."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent.parent
a17 = json.load(open(root / "data/a017_split_unit_controls.json"))
a20 = json.load(open(root / "data/a020_controls.json"))["medians"]
out = root / "next_submission_manuscript/jmgm/figures/holdout_unit_ranking.pdf"

MORGAN, NTO, RDKIT = "#D55E00", "#0072B2", "#009E73"  # Okabe-Ito
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9), constrained_layout=True)
w = 0.26

# panel a
units = [("molecule", "Molecule\n(random)"), ("scaffold", "Scaffold"), ("doi", "Source\npaper")]
q1 = a17["Q1_random_fold_assignment"]
x = np.arange(3)
rdkit_a = [a20["rdkit|molecule"]["rho_all"], np.nan, a20["rdkit|doi"]["rho_all"]]
series = [(-w, [q1[u]["rho_morgan_median"] for u, _ in units], MORGAN, "Morgan"),
          (0.0, rdkit_a, RDKIT, "RDKit 2D"),
          (w, [q1[u]["rho_nto_median"] for u, _ in units], NTO, "NTO")]
for off, v, col, lab in series:
    ax1.bar(x + off, v, w, color=col, label=lab)
    for xi, vi in zip(x + off, v):
        txt, y = ("n.t.", 0.008) if np.isnan(vi) else (f"{vi:.2f}", vi + 0.008)
        ax1.text(xi, y, txt, ha="center", va="bottom", fontsize=6.5)
ax1.set_xticks(x, [l for _, l in units])
ax1.set_ylabel(r"Median Spearman $\rho$ (50 fold assignments)")
ax1.set_xlabel("Holdout unit")
ax1.set_ylim(0, 0.45)
ax1.legend(frameon=False, loc="upper right", ncol=3, fontsize=6.5)
ax1.set_title("a", loc="left", fontweight="bold")

# panel b
groups = [("Multi-paper molecules\n(n = 103)", "rho_multi"), ("Single-paper molecules\n(n = 128)", "rho_single")]
conds = [("molecule", "random", 0.0), ("doi", "paper", 0.0)]
bw = 0.13
xs = np.arange(2)
spec = [("morgan", MORGAN, "molecule", "Morgan, random", None), ("morgan", MORGAN, "doi", "Morgan, paper holdout", "//"),
        ("rdkit", RDKIT, "molecule", "RDKit, random", None), ("rdkit", RDKIT, "doi", "RDKit, paper holdout", "//"),
        ("nto", NTO, "molecule", "NTO, random", None), ("nto", NTO, "doi", "NTO, paper holdout", "//")]
for i, (m, col, u, lab, hatch) in enumerate(spec):
    v = [a20[f"{m}|{u}"][k] for _, k in groups]
    pos = xs + (i - 2.5) * bw
    ax2.bar(pos, v, bw, color=col, hatch=hatch, edgecolor="white" if hatch else col, label=lab, linewidth=0)
    for p, vi in zip(pos, v):
        ax2.text(p, vi + 0.008, f"{vi:.2f}", ha="center", va="bottom", fontsize=5, rotation=90)
ax2.set_xticks(xs, [g for g, _ in groups])
ax2.set_ylabel(r"Median Spearman $\rho$ (50 fold assignments)")
ax2.legend(frameon=False, fontsize=6, loc="upper right", ncol=2)
ax2.set_ylim(0, 0.85)
ax2.set_title("b", loc="left", fontweight="bold")

out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out)
print(out)
