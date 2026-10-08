# Analysis Ledger — Structure-Based ML for Experimental ΔE_ST

Append-only. No number reaches the manuscript except through an entry here.

---

### A-001 .. A-003 — ARCHIVED 2026-09-03
Full entries: `analysis-ledger-archive/ledger-A001-A003-2026-09.md` (grep target, do not read to restore context).
- **A-001** sTDA gap head-to-head, 7 predictors, scaffold GroupKFold → `data/ref1_gap_baseline.json`.
  Raw sTDA gap MAE 0.251 / R² −2.22 / bias +0.221; calibrated 0.104; as a feature, no change.
- **A-002** P1 paired baseline test, P2 DOI-grouped CV (173 papers), P3 split geometry → `data/ref3_p1p2p3.json`.
  P2b: Morgan MAE 0.102 / R² 0.184 / ρ 0.155 · NTO 0.099 / 0.142 / **ρ 0.346**; Morgan-vs-mean CI crosses zero.
  P3: 92.9%% singleton scaffolds, scaffold CV ≈ random CV, DOI CV harder.
- **A-003** label-noise reconciliation → `data/label_noise_reconciled.json`.
  Inter-report RMS SD 0.072 eV; **RMS standard error of the median label 0.049 eV** — the number to
  quote for how precisely the target is known. States it does NOT support 'label noise alone explains the error'.

### A-004 — ARCHIVED 2026-09-04
Full entry: `analysis-ledger-archive/ledger-A004-2026-09.md` (grep target).
Referee-3 points 4-7 verdicts. Key numbers still cited: **P4** expansion to 728 raises MAE 0.098->0.124,
R2 0.18->-0.01, rho 0.34->0.16 (caveat i: different structure pipeline, confounded). **P5** every one of the
top-100 conformal intervals spans a negative gap; **mean ratio of interval WIDTH to predicted value = 7.6**
(half-width ratio is 3.8 — do not conflate). **P6** vertical vs adiabatic n=27/14 mol: MAE 0.195, RMSE 0.245,
bias +0.071, r2 0.704, rho 0.711, slope 1.002, intercept 0.070, max dev 0.548. **P7** 66 IST records/61 names,
n_with_smiles=0; 5 IST molecules in corpus.

### A-005 — ARCHIVED 2026-09-04
Full entry: `analysis-ledger-archive/ledger-A005-2026-09.md` (grep target).
TD-DFT reference, wB97X-D4/def2-SVP RIJCOSX TDA, 231/231 converged, 604.5 core-h. Key numbers:
raw TD-DFT gap MAE **0.654** / R² **-16.3** / signed bias **+0.654** (worst predictor in the study);
calibrated in-fold **0.099**; Morgan+gap **0.087** vs Morgan 0.091 (ΔMAE +0.0039, CI [+0.0004,+0.0071],
Wilcoxon p 0.016); NTO+gap 0.092 vs 0.096. Pearson r vs experiment 0.43 (not >0.95). Claim C9.

### A-006 — ARCHIVED 2026-09-05
Full entry: `analysis-ledger-archive/ledger-A006-2026-09.md` (grep target).
T1 Morgan-vs-NTO paired equivalence: dMAE -0.0043, CI [-0.0171,+0.0092], Wilcoxon p 0.657,
**margin 0.017 eV**. T2 DOI split: NTO rho 0.346 vs Morgan 0.155, **drho +0.191 CI [+0.037,+0.344]**
(caveat iv: resamples molecules not DOI groups). T3 enrichment -> superseded by A-006b.
T4 ceiling audit: stored rss_floor 0.1282 > rf_mae 0.0956 (incoherent); defensible label
precision **0.049 eV**. T5 clean-label -> WITHDRAWN by A-008 (constant predictor reproduces it).

### A-006b — enrichment reconciliation (archived during verification B)
Full entry: `analysis-ledger-archive/ledger-before-verification-B.md` (search target).
Threshold is <=0.1 eV, base rate 0.4935. Morgan top-10/top-50 molecules: precision
0.600/0.680, enrichment 1.22/1.38. Morgan/NTO top-20%: 1.37/1.59.
Interpretation: counts and fractions are different operating points; do not relabel them.
Claim C13; no evidence defect under the original conventions.

### A-006c — library provenance (archived during verification B)
Full entry: `analysis-ledger-archive/ledger-before-verification-B.md` (search target).
PubChem substructure filtering, not fragment enumeration; 22194/22194 and 100/100
pass the code's actual donor/acceptor SMARTS. Top-100: 97 halogen, 53 stereocentre.
Interpretation: reviewer used a narrower printed list; library definition passes but
chemical breadth and composition must remain disclosed. Claim C14.

### A-007 — leakage and nested CV (archived during verification B)
Full entry: `analysis-ledger-archive/ledger-before-verification-B.md` (search target).
Source: code/a007_leakage_and_nested_cv.py, data/a007_leakage_nested_cv.json.
Interpretation: computed-energy difference leakage depends on learner capacity to subtract;
nested model selection must be distinguished from flat comparisons. Claims C15/C16.
Numbers and methods remain in the complete archived entry; no claim superseded here.

### A-008 — Two control checks that INVALIDATE A-006 T5 and correct the 0.072 eV label figure
- **Date:** 2026-09-04. **Trigger:** L3 beat-2 review (`review-10-critical-thinking-beat2.md`, CRITICAL 1
  and HIGH 4). Both challenged results this project had already written into the manuscript. Both are
  confirmed by direct recomputation, so both of my earlier readings were wrong.
- **Source:** recomputation on the A-006 folds (seed 0, scaffold GroupKFold(5), NTO features);
  `data/label_noise_reconciled.json` read directly.

**(1) The clean-label ordering is reproduced by a feature-free predictor — A-006 T5 is not evidence.**

| subset | n | model MAE | constant-predictor MAE | SD of y in subset |
|---|---|---|---|---|
| replicated, clean (sd $\le$ 0.05) | 79 | 0.0741 | **0.0739** | 0.1008 |
| replicated, noisy (sd $>$ 0.05) | 23 | 0.1062 | **0.1070** | 0.1525 |
| singletons (no replicate) | 129 | 0.1069 | **0.1099** | 0.1952 |

  The random forest ties a constant predictor on every subset. The 0.074 / 0.106 / 0.107 ordering
  tracks the dispersion of the target within each subset (SD 0.10 / 0.15 / 0.20), which is what a
  mean-absolute-error must do, and carries no information about label quality.

**(2) The 0.072 eV label-scatter figure does not reproduce.** `data/label_noise_reconciled.json`
  (n = 102 multi-report molecules) contains: RMS per-molecule SD (ddof=1) **0.0856**; (ddof=0) **0.0701**;
  RMS single report about the median **0.0905**; RMS single report about the mean 0.0875;
  RMS standard error of the median **0.0490**; mean absolute deviation about median 0.0312.
  **0.072 is not among them.** The manuscript attaches 0.072 to the sentence "a single report sits
  X (RMS) from its molecule's median", whose computed value is **0.0905** — a 26 % understatement.

- **Interpretation.**
  (1) A-006 T5 must be withdrawn as support for the label-limited account, and every manuscript
  sentence drawing on it removed. What the control actually shows is broader and worth stating: on
  this task the forest does not beat a constant predictor *within* any label-quality stratum, so the
  0.096 eV headline rests on between-stratum variation. This strengthens the paper's own cautionary
  thesis while destroying the specific argument I had built. The honest residual claim about the
  ceiling is A-006 T4's 0.049 eV precision of the median, which reproduces exactly.
  (2) Quote **0.0905 eV** for single-report scatter about the median, or **0.049 eV** for the precision
  of the median actually regressed on, and name the estimator in the sentence. 0.072 must not appear.
- **Caveats.** (i) The constant-predictor comparison is within-subset and unpaired; it shows the model
  adds nothing *within* a stratum, not that it adds nothing overall (A-002 P1 shows it beats the global
  mean baseline on paired folds, CI excluding zero). (ii) n = 23 for the noisy stratum throughout.
- **Claims:** **withdraws C12** (A-006 T5). Supersedes the 0.072 figure wherever it appears.
  New **C17**: within every label-quality stratum the model ties a constant predictor.

### A-009 — Three further claims checked against the deposited data (beat-2 open findings)
- **Date:** 2026-09-04. **Trigger:** remaining OPEN/CRITICAL items from L3 beat 2
  (`review-9-adversarial-beat2.md`, `review-10-critical-thinking-beat2.md`, `review-8-verification-beat2.md`).
- **Sources read directly:** `data/ref3_p1p2p3.json` (P2c), `MCA_submission/tables/experimental_delta_est_extracted.csv`.

**(1) The 0.080 eV multi-report MAE has a matched baseline that was not reported.**
  `P2c_report_subsets.multi_report`: $n=102$, Morgan MAE **0.08013**, mean-baseline MAE **0.08908**,
  $\Delta$MAE **0.00896** with a CI spanning zero. The subset's low absolute error reflects its
  narrower label spread, not better prediction — the same artefact that invalidated A-006 T5 (A-008).
  Single-report subset for contrast: $n=129$, MAE 0.1006 vs baseline 0.1229, $\Delta$ 0.0223, CI excludes zero.

**(2) "The DOI split holds out a measurement protocol" has no metadata behind it.**
  In the source table, `phase` is non-null for **0/1520** rows and `atmosphere` for **0/1520**.
  `record_method` is non-null for 1520/1520 but its values are text-parser names
  (`AutoTableParser`, `QuantityModelTemplateParser`, …), not experimental protocols.
  The DOI split holds out a source paper, which is a *proxy* for shared laboratory practice.
  This confirms A-002 caveat (i) and the manuscript now says so wherever the split is described.

**(3) The 79 % vs 30 % SHAP contrast is a dimensionality artefact.**
  Ten descriptors are \qty{29}{\percent} of the 35-dimensional NTO set; ten bits are
  \qty{0.5}{\percent} of a 2048-bit fingerprint. Concentration of attribution follows from feature
  count and is not evidence of superior interpretability. What survives is that NTO variables are
  named physical quantities while Morgan bits are hashed substructures — a property of the
  representation, not of the attribution totals.

- **Interpretation.** All three were stated in the manuscript in a form stronger than the data
  supports, and all three are now corrected in place rather than deleted: (1) the matched baseline
  and the zero-spanning CI are reported alongside the 0.080 figure; (2) "measurement protocol" is
  replaced by "a laboratory's shared practice … a proxy, since the corpus records no measurement
  conditions per value"; (3) the attribution totals are given with the denominators that make them
  comparable, and the interpretability case rests on named variables instead.
- **Caveats.** (i) P2c subsets are not independent of A-002 P2a — the same 102 molecules define both.
  (ii) The parser-name finding means no protocol-controlled analysis is possible in this corpus at all;
  a solvent-controlled study would need a different data source.
- **Claims:** qualifies C5 (no-transfer) — the split tests source-paper transfer, not protocol transfer.
  Withdraws the interpretability-superiority reading of the SHAP contrast.

### A-008b — CORRECTION: A-008 claim (2) was wrong; 0.072 eV does reproduce
- **Date:** 2026-09-04. **Supersedes A-008 (2)**, which asserted that 0.072 eV "is not among"
  the computed estimators and had it replaced by 0.0905 in 14 places and in `qa_score.py`.
  That was my error, and it was already answered by A-003 before I made it.
- **Source:** direct recomputation over the 102 multi-report molecules, and the reconciliation
  table in the archived A-003 (`analysis-ledger-archive/ledger-A001-A003-2026-09.md`).
- **Result.** Both figures are correct estimators of the same quantity under different weightings:

| estimator | value |
|---|---|
| RMS of a single report about its molecule's median, **each molecule weighted equally** | **0.0724** |
| the same quantity, **report-pooled** (each report weighted equally) | 0.0905 |

  A-003 selected the molecule-weighted figure and recorded "the published 0.072 stands".
- **Interpretation.** The published value is right, and the defect was never the number: it was
  that the manuscript printed "0.072 (RMS)" without naming the weighting, so a recomputation that
  reasonably chose report-pooling appeared to contradict it. Two independent reviewers and I all
  landed on 0.0905 for exactly that reason. The fix is to state the estimator, which the manuscript
  now does at both primary sites; report-pooling over-weights molecules carrying many reports, while
  the regression target is one median per molecule, so molecule-weighting is the apt choice.
  All 14 substitutions have been reverted and the `qa_score.py` A4 assertion restored to 0.072.
- **Lesson recorded for this project:** a number that "fails to reproduce" should be checked against
  the ledger's own reconciliation entries before the manuscript is changed. A-003 existed and I did
  not consult it.
- **Claims:** restores C3 as originally stated, now with the estimator named. A-008 claim (1) — the
  constant-predictor control that withdrew the clean-label result — is unaffected and still stands.

### A-010 — CORRECTION: the TD-DFT cost figures were a projection, not a measurement
- **Date:** 2026-09-05. **Trigger:** L3 beat-4 numerical audit (`review-15-numerical-audit.md`, M1–M3).
- **Corrects archived A-005**, which recorded "604.5 core-hours / median 18.2 min / slowest 1.13 h"
  as the measured cost of the TD-DFT run. Those three values are the contents of
  `data/tddft_reference_timing.json → summary.projection_231`, a block whose own name says it is a
  power-law extrapolation (`wall_seconds = 0.2634 · atoms^1.89`). The manuscript then stated them as
  observed run statistics in **12 places across 6 files**.
- **Measured values**, from the 231 per-molecule `wall_seconds` in `data/tddft_reference_231.csv`
  and `summary.elapsed_hours`:

| quantity | published (projection) | measured |
|---|---|---|
| total cost | 604.5 core-hours | **609** core-hours (8 processes × 76.15 h elapsed); Σ per-molecule × 8 = 631.9 |
| median per molecule | 18.2 min | **16.05 min** |
| slowest molecule | 1.13 h | **2.087 h** (4F-m-ν-DABNA) |

- **Interpretation.** The projection was close on the total (0.8 % low) but wrong by 1.85× on the
  slowest molecule, and the manuscript presented all three as measurements. The corrected text now
  quotes 609 core-hours with its basis stated (elapsed wall time × process count) and says explicitly
  that these are measured times rather than a scaling extrapolation. The scientific conclusion is
  unchanged — the point was ever that a several-hundred-core-hour reference buys 0.004 eV — but the
  cost of a computation is exactly the kind of number a referee may try to reproduce.
- **Caveats.** (i) 609 (elapsed × procs) and 632 (Σ per-molecule × 8) differ because the two measure
  occupancy and summed job time; 609 is the resource footprint and is what is published.
  (ii) The power-law fit itself remains useful for extrapolating to larger molecules and is retained
  in the deposited timing file, now clearly distinguished from the measured values.
- **Claims:** amends C9's cost clause. No accuracy result depends on these numbers.

### A-011 — The nearest available condition-controlled test (answers Referee 3, point 2)
- **Date:** 2026-09-05. **Trigger:** the one Digital Discovery objection still answered only by
  explanation. Referee 3 asked for performance "on a strictly solvent- and condition-controlled
  subset rather than attributing model failure entirely to external data noise".
- **Why a proxy:** `phase` and `atmosphere` are non-null for **0/1520** source rows (A-009), so a
  strictly controlled subset cannot be constructed from this corpus at all. The closest available
  proxy is the source paper: molecules whose reported values all come from one DOI were measured
  under a single laboratory's conditions; molecules drawing on two or more DOIs certainly mix them.
- **Source:** `code/a011_condition_controlled_proxy.py` → `data/a011_condition_controlled.json`.
  Same scaffold GroupKFold(5), seed 0, 2000 bootstrap resamples.

| group | n | Morgan MAE | NTO MAE | constant predictor | target SD | advantage over constant |
|---|---|---|---|---|---|---|
| single-source (1 DOI) | 218 | 0.0933 | 0.0965 | 0.1004 | 0.1703 | **+0.0071** |
| multi-source (≥2 DOIs) | 13 | 0.0584 | 0.0802 | 0.0608 | 0.0814 | **+0.0024** |

  Difference in advantage-over-constant (single − multi) = **+0.0047 eV**,
  95% CI **[−0.0138, +0.0345]** — **does not exclude zero**.

- **Interpretation.** The test returns a null. Restricting to molecules measured under one
  laboratory's conditions does not measurably improve what the model adds over a feature-free
  predictor. Two details matter for reading it honestly. First, the multi-source group has the
  *lower* raw MAE (0.058 vs 0.093) — the opposite of what a condition-noise account predicts — but
  its target dispersion is also half as large (SD 0.081 vs 0.170) and its own constant predictor
  scores 0.061, so the raw difference is dispersion, not label quality. This is precisely the
  artefact that invalidated A-006 T5, and the comparison was therefore run against a per-group
  constant predictor from the outset (A-008). Second, only **13** molecules draw on more than one
  source paper, so the test has little power; it constrains the effect rather than excluding it.
  The defensible statement is that we ran the closest test the data permits and found no evidence
  that mixed measurement conditions drive the error — not that conditions are irrelevant.
- **Caveats.** (i) n = 13 in the multi-source group. (ii) One DOI means one laboratory, not one
  solvent: a single paper may itself report several media. (iii) Single-source includes the 129
  single-report molecules, whose condition consistency is untestable rather than verified.
- **Claims:** new **C18** (nearest-available condition-controlled test is null, CI includes zero).
  Converts Referee 3 point 2 from "cannot be tested" to "tested by proxy, null, with the reason the
  strict test is impossible stated".

### A-012 — Verification B diagnostics: A-011 control and ceiling arithmetic
- Status: diagnostic verification only; no replacement manuscript result authorised.
- Question: is A-011's constant comparator held out, and is the stored RSS arithmetic valid?
- Source/method: outputs/critical-reviews/verification-B-evidence/parent_check.py;
  code/_dataset.py; original experimental table; GroupKFold(5), same scaffold groups.
  Output: outputs/critical-reviews/verification-B-evidence/parent-check.json.
- Result: A-011 group counts reproduce (218 single-source, 13 multi-source).
  In-sample group-median MAEs reproduce as 0.1004288991 and 0.0608076923 eV.
  Fitting each group median on training-fold members only gives held-out MAEs
  0.1006444954 and 0.0692692308 eV respectively. Deterministic diagnostics; no CI estimated.
  RSS of stored 0.072, 0.195, 0.2964 eV components is 0.3620248058 eV, not stored
  0.1282 eV; stored forest MAE is 0.0956 eV.
- Interpretation: the published comparator is the empirical target-dispersion minimum,
  not a held-out predictor. It is usable as a descriptive statistic if named as such,
  but its comparison with out-of-fold forest errors is not a matched predictive benchmark.
  This matters most in the small multi-source group. These calculations do not establish
  a new model advantage or confidence interval, nor an effect of measurement conditions.
  The RSS field is arithmetically inconsistent; even the recomputed sum would not be a
  justified experimental noise floor because its components mix different quantities.
- Caveats: no model retraining, no bootstrap rerun, no group-correlation correction;
  source-paper membership does not identify a laboratory or controlled conditions.
- Claim links: flags C18/A-011 and the inherited ceiling decomposition; not a prose repair.
- Independent interpretation check: BLOCKED (Hermes credit exhaustion; Codex token failure).
  No independent report exists; see verification-B-summary.md. A-012 remains diagnostic only.

### A-013 — Pretrained molecular-transformer baseline (DA-5, L3 pass 2)
- **Date:** 2026-09-10. **Trigger:** L3 novelty death-angle DA-5
  (`divergence-death-angle-L3-novelty-2026-09-10.md`): the claim "data quality, not
  model capacity, is the limit" had been tested only with tabular learners on
  hand-built vectors; no pretrained encoder or GNN was run on this task. Author
  authorised one frozen-encoder baseline.
- **Source:** `code/a013_pretrained_encoder_baseline.py` -> `data/a013_pretrained_encoder.json`.
  Encoder `DeepChem/ChemBERTa-77M-MLM` (RoBERTa, 384-dim), frozen, mean-pooled over the
  attention mask. Downstream: the canonical `rf()` (400 trees, seed 0) and a RidgeCV head
  fitted inside each fold. Same `GroupKFold(5)` on Bemis-Murcko scaffolds, seed 0, as A-006.

| model | MAE (eV) | R2 | Spearman rho |
|---|---|---|---|
| Morgan-RF (reproduced) | 0.0913 | +0.274 | 0.306 |
| NTO-RF (reproduced) | 0.0956 | +0.253 | 0.356 |
| ChemBERTa-77M frozen + RF | **0.1128** | +0.062 | 0.206 |
| ChemBERTa-77M frozen + RidgeCV | **0.1086** | +0.104 | 0.222 |

  Paired (identical folds): ChemBERTa+RF minus Morgan-RF, delta MAE **+0.0215 eV**
  (95% CI **[+0.0077, +0.0344]**, excludes zero; Wilcoxon p reported in JSON) -- the
  pretrained representation is significantly *worse*. Its best head (0.109 eV) does not
  beat the mean-value baseline (0.107 eV).
- **Interpretation.** A pretrained SMILES-transformer embedding, the representation a
  2025-26 reader would reach for first, carries less usable signal for experimental
  dEST on this corpus than a 2048-bit Morgan fingerprint. The Morgan/NTO reference
  numbers reproduce the published 0.091/0.096 exactly, so the CV recipe matches the
  rest of the study. This removes DA-5: the "capacity is not the missing ingredient"
  statement now holds against a pretrained encoder, not only tabular learners.
- **Caveats.** (i) Frozen-embedding protocol. Fine-tuning a 77M-parameter encoder on
  231 labels is guaranteed to overfit and is not a meaningful capacity test at this
  scale; not run. (ii) One encoder checkpoint (ChemBERTa-77M-MLM); a larger model
  (MoLFormer-XL) was not tested. (iii) A message-passing GNN trained end-to-end on 231
  molecules is also expected to overfit and was not run; the claim is about
  representations usable at this data scale, not about deep learning in general.
- **Claims:** strengthens C-capacity (no higher-capacity learner or pretrained
  representation beats the forest). New wording in results subsec:ceilings and
  Limitations bullet (6). Manuscript number: 0.109-0.113 eV for the pretrained encoder.

### A-014 — Bootstrap CI on triage/enrichment precision (audit item 10.3, AUDIT_REPORT.md 2026-09-12)
Full entry: `analysis-ledger-archive/ledger-A014-2026-09.md` (grep target). Summary: triage enrichment bootstrap CI — top-10/top-20 not distinguishable from 1.0×; top-50/top-20% exclude 1.0×. Triage now in SI.

### A-015 — Max achievable R² (Crusius) — SUPERSEDED by A-018
Full entry: `analysis-ledger-archive/ledger-A015-2026-09.md` (grep target). Values 0.71–0.91 no longer used; see A-018 (0.80–0.86).
### A-016 — Paper-unit uncertainty and repeated source-paper partitions for ranking
- **Date:** 2026-10-06. **Source:** `code/jmgm_doi_group_ranking_check.py` -> `data/a016_jmgm_doi_group_ranking.json`. Reconstructs A-006's primary-DOI assignments, models, and five GroupKFold out-of-fold predictions exactly; assertions reproduce $n=231$, 173 DOI groups, and the original Morgan/NTO rank correlations to 1e-10. `data/a016_jmgm_doi_group_input.csv` holds the 231 row-level modelling inputs, 35 NTO features, SMILES and primary DOI. The standalone `jmgm_manuscript/supplementary_data/` copy reruns with `--input-csv` and produced an output JSON byte-identical to the original; its CSV reader needs `float_precision="round_trip"` because the default parser introduced 1e-16 differences in 62 labels and shifted the fitted rank correlation. The three original evidence files are staged in Git but not committed or uploaded.
- **Paper-unit interval:** 5000 paired percentile bootstrap draws of **173 primary DOI groups**, resampling each selected paper's molecules together and recomputing both Spearman correlations from the fixed out-of-fold predictions. Seed 54321. Morgan $\rho=0.154717$, NTO $\rho=0.345913$, NTO−Morgan $\Delta\rho=+0.191196$; DOI-cluster 95% percentile interval **[+0.014317,+0.366034]**. The A-006 molecule-level interval was [+0.037444,+0.344206]. The cluster interval barely excludes zero and does not propagate model refitting or alternative fold allocation.
- **Partition stability probe:** 20 GroupShuffleSplit train/test divisions (seed 17; 20% of paper groups held out each time), refitting each representation on each training division. NTO $\rho$ exceeds Morgan in **19/20** divisions; median $\Delta\rho=+0.13355$, range [−0.01418,+0.40378]. Repeated tests overlap; 19/20 is a descriptive stability count, **not** a probability, independent confirmation, calibrated interval, external source transfer, or evidence about laboratories.
- **Robustness bound:** leaving out each primary DOI group from the original fixed out-of-fold predictions one at a time gives $\Delta\rho$ range [+0.15347,+0.21651]. This is a sensitivity check, not model retraining.
- **Interpretation:** evidence for a corpus-conditional accuracy/ranking disagreement survives paper-cluster resampling and is mostly positive under repeated source-paper partitions, but it remains a retrospective analysis of the same 231 emitters. No external or protocol-controlled data and no prospective top-k ranking validation; do not make a universal or laboratory-transfer claim. The paper-unit interval does not turn the finding into established molecular-design utility.

### A-018 — Noise ceiling with per-molecule label noise (supersedes A-015's 0.91 upper value)
- **Date:** 2026-10-06. **Trigger:** L3 beat 1, pass-3 H4: A-015 applied the RMS SEM of the median (0.049 eV, computed on the 102 multi-report molecules only) to all 231 molecules; A-015 also has no source script in the repo or the 2026-09-10 backup.
- **Method.** `code/a018_noise_ceiling_heterogeneous.py`. noise_var_i = SD_i²/n_i for 102 multi-report molecules (reproduces 0.049 eV RMS, asserted); σ_report² for 129 single-report molecules. R²_max = 1 − mean(noise_var)/var(y). Molecule bootstrap 5000, seed 54321. Pre-registered in docstring before running. **Source:** `data/a018_noise_ceiling_heterogeneous.json`, log `outputs/extract/a018_run.txt`.

| σ_report for single-report molecules | RMS noise | R²_max | 95% CI | % captured (R²=0.257) |
|---|---|---|---|---|
| 0.0856 eV pooled within-molecule SD (**primary**) | 0.072 eV | **0.81** | [0.68, 0.88] | 31.5% |
| 0.0905 eV single report about median | 0.075 eV | 0.80 | [0.65, 0.87] | 32.2% |
| 0.072 eV molecule-weighted scatter | 0.063 eV | 0.86 | [0.75, 0.91] | 29.9% |

- **Interpretation.** Treating single-report labels as noisy as one report lowers the ceiling from A-015's 0.91 to 0.80–0.86. A-015's lower value (0.71) applied single-report noise to every molecule, including the 102 with replicates, so it is too pessimistic; 0.91 applied replicate precision to molecules without replicates, so it is too optimistic. The bracket the manuscript should carry is **R²_max ≈ 0.80–0.86** (primary 0.81, CI 0.68–0.88). The qualitative conclusion stands: the model captures about 30–32% of the ceiling, so label noise is not the dominant constraint.
- **Caveats.** (i) Additive, structure-independent noise assumed. (ii) Report independence assumed; 57.8% of multi-report molecules have zero spread (possible re-citation), which would make their SEM too small and R²_max too high. (iii) SD/√n is the SEM of a mean; the median's SEM is larger, so the multi-report term is a lower bound. (iv) σ_y here 0.1666 eV from labels (A-015 used 0.1670 from the predictions file); the difference does not move any reported digit. (v) CI resamples molecules, not the noise estimate itself.
- **Consequence for manuscript.** Replace "R²_max ≈ 0.71–0.91" with "≈ 0.80–0.86" in abstract, introduction, results, and conclusions; add an SI paragraph stating assumptions (i)–(v). "27–35%" becomes "about 30–32%".
- **Claims:** C-019 updated (A-015 → A-018).

### A-017 — Holdout-unit controls for the NTO−Morgan rank contrast
- **Date:** 2026-10-06. **Trigger:** L3 beat 1 pass-2 H1 (Δρ 0.191 from one fixed DOI split) and pass-1 H1/H2 (split vs criterion confounded; no matched control for 19/20).
- **Method.** `code/a017_split_unit_controls.py` (pre-registration + decision rule in docstring, written before running). Input `data/a016_jmgm_doi_group_input.csv` (231 molecules; 173 DOI, 212 scaffold groups). Self-check reproduces A-006 T1/T2 ρ to 1e-10 and A-016's 19/20. Q1: 50 random group→fold assignments (5 folds) per holdout unit {DOI, scaffold, molecule}, RF(400) OOF for Morgan and NTO. Q2: A-016 repeated holdout (GroupShuffleSplit 20, test 0.2, seed 17) per unit. Q3: scaffold-cluster bootstrap (5000) on fixed A-006 scaffold folds. **Source:** `data/a017_split_unit_controls.json`; log `outputs/extract/a017_run.txt`.

| holdout unit | Δρ median (50 seeds) | 2.5–97.5% of seeds | Δρ>0 | median ρ Morgan | median ρ NTO | ΔMAE Morgan−NTO median [2.5,97.5] | Q2 NTO higher /20 |
|---|---|---|---|---|---|---|---|
| source paper (DOI) | **0.160** | [0.073, 0.246] | 50/50 | 0.172 | 0.337 | −0.0004 [−0.0061, 0.0068] eV | 19 |
| scaffold | 0.073 | [0.007, 0.144] | 49/50 | 0.269 | 0.351 | −0.0045 [−0.0091, 0.0010] eV | 11 |
| molecule (random) | 0.055 | [−0.033, 0.126] | 45/50 | 0.305 | 0.352 | −0.0081 [−0.0157, −0.0002] eV | 13 |

  Q3 fixed scaffold split: Δρ 0.051, scaffold-cluster 95% CI [−0.099, 0.202], 25.7% of draws ≤ 0.
- **Interpretation.** Pre-registered rule met (DOI median exceeds scaffold by 0.087 and molecule by 0.105; DOI 50/50 positive vs 49/50 and 45/50): **the holdout unit, not the error-vs-ranking criterion, drives the size of the contrast.** The mechanism is asymmetric: NTO rank correlation is nearly invariant to the holdout unit (median ρ 0.337–0.352), while Morgan's falls from 0.305 (random) to 0.269 (scaffold) to 0.172 (source paper). Fingerprints carry rank information that partly does not survive holding out whole papers (consistent with series- or paper-specific structure–label association; cause not tested); the semi-empirical descriptors' rank information does. A-016's 0.191 lies in the upper part of the 50-seed DOI distribution; the fold-assignment-robust value is ≈0.16 (seed range 0.07–0.25). Under scaffold holdout the contrast is small (≈0.07) and, on the fixed split, not resolved (CI spans 0). For error, Morgan is slightly better under random holdout (−0.008 eV, upper bound −0.0002) and indistinguishable under scaffold or paper holdout.
- **Caveats.** Seed percentiles describe fold-assignment variability, not a confidence interval; all units use the same 231-molecule corpus; 128/173 DOI groups are single molecules, so paper holdout differs from random mainly through the 45 multi-molecule papers; no lab-transfer test; "series-specific" reading is a hypothesis.
- **Consequence for manuscript.** Central question and conclusions should name the holdout unit (and the Morgan-degradation / NTO-invariance asymmetry) as the finding; report Δρ ≈ 0.16 (seed range) alongside or instead of 0.191; add matched 11/20 and 13/20 next to 19/20; "Accuracy and ranking answer different questions" (conclusions:15) must be rewritten. Re-scoping decision → author + `bmad-correct-course`.
- **Derived (L2-B, 2026-10-06):** Morgan median ρ falls by 43.6% from random to source-paper holdout ((0.305 − 0.172)/0.305); used as "44%" in highlight 2.
- **Claims:** new C-020 (holdout-unit dependence; NTO invariance vs Morgan degradation).

### A-019 — Series/publication shortcut test (Morgan vs NTO)
- **Date:** 2026-10-06. **Trigger:** A-017 asymmetry (Morgan ρ falls under paper holdout, NTO does not); author asked to test the shortcut hypothesis.
- **Method.** `code/a019_series_shortcut_test.py` (hypothesis + decision rule in docstring, fixed before running). Same 50 fold assignments as A-017 for molecule (random) and DOI holdout. T1: ρ on 103 molecules from multi-molecule papers vs 128 from single-molecule papers. T2 "sibling pull": Spearman between OOF prediction and leave-self-out mean label of same-paper siblings (103 molecules). **Source:** `data/a019_series_shortcut_test.json`; log `outputs/extract/a019_run.txt`.

| median over 50 seeds | Morgan random | Morgan paper | NTO random | NTO paper |
|---|---|---|---|---|
| ρ, multi-paper molecules (103) | 0.538 | **0.264** | 0.483 | 0.446 |
| ρ, single-paper molecules (128) | **0.097** | 0.099 | 0.243 | 0.245 |
| sibling pull (103) | 0.676 | 0.255 | 0.471 | 0.418 |

  Random-minus-paper medians [2.5–97.5% of seeds]: Morgan ρ_multi 0.269 [0.174, 0.413], ρ_single −0.003 [−0.121, 0.139], sibling pull 0.415 [0.276, 0.606], MAE_multi −0.025 [−0.034, −0.014] eV; NTO ρ_multi 0.040 [−0.026, 0.111], ρ_single −0.015, sibling pull 0.048 [−0.014, 0.138], MAE_multi −0.007 [−0.013, 0.001] eV.
- **Interpretation.** Verdict **SUPPORTED** (criteria a and b both met). Morgan's rank information sits almost entirely on molecules that have same-paper siblings: ρ 0.54 when siblings are in training, 0.26 when the whole paper is held out, and ≈0.10 on single-paper molecules under either holdout. Its predictions track the siblings' mean label (pull 0.68) and lose that tracking under paper holdout (0.26). Its MAE advantage on multi-paper molecules (0.025 eV) disappears under paper holdout. NTO shows none of these shifts (ρ_multi 0.48→0.45; single-paper ρ 0.24 under both) and ranks single-paper molecules better than Morgan in every setting. The fingerprint model's apparent ranking skill under random or scaffold holdout is largely a within-publication (series/lab) association; the semi-empirical NTO descriptors carry rank information that does not depend on having seen the publication.
- **Caveats.** Association, not demonstrated cause: siblings share both close structure and measurement protocol, and this test cannot separate structural near-duplication from lab/protocol bias. Same corpus, retrospective. Subset ρ on 103/128 molecules is noisier than whole-set ρ. Seed ranges are fold-assignment variability, not CIs.
- **Consequence for manuscript.** Supports the F3-family central question (NTO retains rank information that fingerprints lose on unseen publications) with a tested mechanism. Report single-paper-molecule ρ (Morgan ≈0.10 vs NTO ≈0.24) as the cleanest contrast. Warn modellers that random/scaffold splits on literature-mined TADF data overstate fingerprint ranking skill.
- **Claims:** new C-021 (publication-shortcut association for Morgan; NTO invariance).

### A-020 — RDKit control, within-paper ranking, single-paper contrast interval
- **Date:** 2026-10-06. **Trigger:** L3 beat 2: pass-2 H-B (no non-quantum dense-descriptor control), pass-1 H-new-1 (ranking within vs across held-out papers), H-new-2 (single-paper contrast without interval).
- **Method.** `code/a020_controls.py` (three pre-registered rules in docstring, fixed before running). Same 50 fold assignments as A-017/A-019 (molecule, DOI holdout); Morgan, NTO, 217 RDKit 2D descriptors. Within-paper ρ = Spearman of paper-centred label vs paper-centred prediction (103 multi-paper molecules). Single-paper interval: OOF predictions averaged over the 50 DOI-holdout seeds, paired molecule bootstrap 5000 (seed 54321). **Source:** `data/a020_controls.json`; log `outputs/extract/a020_run.txt`.

| median (50 seeds) | ρ all, random | ρ all, paper | ρ multi, random → paper | ρ single, paper | within-paper ρ, random | within-paper ρ, paper |
|---|---|---|---|---|---|---|
| Morgan | 0.305 | 0.172 | 0.538 → 0.264 | 0.099 | −0.354 | −0.063 |
| RDKit 2D | 0.303 | 0.209 | 0.518 → 0.350 | 0.102 | −0.223 | −0.027 |
| NTO | 0.352 | 0.337 | 0.483 → 0.446 | 0.245 | −0.041 | 0.029 |

  RDKit multi-paper drop 0.176 [0.050, 0.281]. Single-paper contrast (seed-averaged predictions): ρ NTO 0.266, Morgan 0.108, RDKit 0.114; Δ NTO−Morgan 0.158, 95% CI [−0.052, 0.374], 8.0% of draws ≤ 0.
- **Interpretation.** Q1 verdict **FINGERPRINT-LIKE**: the dense non-quantum RDKit set loses rank information under paper holdout like Morgan does (drop 0.18, single-paper ρ 0.10); among the three tested sets, holdout-unit invariance is specific to the excited-state NTO descriptors, so dimensionality or smoothness alone does not explain it. Q2: no representation ranks molecules *within* a held-out paper (within-paper ρ −0.06 to 0.03); the NTO advantage is ranking *across* unseen publications/series. Under random holdout, Morgan's within-paper ρ is −0.35 — its predictions pull towards the sibling mean, which inverts the order inside a paper. Q3: the single-paper contrast (0.27 vs 0.11) is not resolved (CI spans 0); it must leave the abstract and appear in Results only with its interval.
- **Caveats.** Seed ranges = fold-assignment variability; CI in Q3 conditional on seed-averaged fixed predictions; same corpus; within-paper ρ pools 45 small papers (2–4 molecules), low power.
- **Consequence for manuscript.** Wording: "ranks molecules from held-out papers" → "orders emitters across publications absent from training; neither representation orders molecules within a held-out paper". Add RDKit control sentence + table row (strengthens NTO specificity). Remove 0.10 vs 0.24 from abstract; report Q3 with CI in Results.
- **Claims:** C-021 refined (across-publication, not within); new C-022 (NTO invariance not shared by RDKit 2D).

### A-021 — Provenance entry for pre-existing JSON results cited in the JMGM text (L4 number sweep)
- **Date:** 2026-10-06. **Trigger:** L4 audit number sweep (`outputs/l4/number_sweep.csv`, `orphan_triage.csv`): 12 manuscript numbers traced only to `data/*.json`, no ledger ID. No recomputation; values read from the files below.

| manuscript value | quantity | source (key) |
|---|---|---|
| R² 95% CI RDKit −0.11 to 0.51; Morgan −0.19 to 0.53 | scaffold-CV R² bootstrap CI | `data/final_model_metrics.json` feature_set_comparison.{RDKit,Morgan}.R2_CI [−0.107, 0.513], [−0.192, 0.526] |
| 99.5% | S1_S_he rank-1 SHAP stability, reliable subset, 200 bootstrap refits | `data/conditional_shap.json` rank1_stability_reliable_n200 0.995 (full corpus 0.795) |
| 152 | reliable-prediction subset size | `data/conditional_shap.json` n_reliable |
| 19.6% | fraction of molecules whose label noise exceeds model MAE | `data/ceiling_decomposition.json` A_frac_exceeding_model_mae 0.196 |
| 0.70 | sTDA vs ORCA adiabatic Pearson r² | `data/ceiling_decomposition.json` B_pearson_r2 0.704 |
| 0.169 eV | MLP (100,50) MAE, capacity check | `data/ceiling_tests.json` MLP (100,50).MAE |
| 70 / 38 | coarsest Butina clustering: groups / largest cluster | `data/cluster_cv.json` n_clusters 70, largest_cluster 38 |
| shuffle-null mean 0.121 eV | 1000 target shuffles, scaffold CV | `data/scaffold_split_analysis.json` null_mae_mean 0.121 (manuscript results:39 corrected from 0.120 at L4) |

- **Interpretation.** These are supporting checks already in the manuscript before the F3-b rewrite (capacity, clustering robustness, SHAP stability, ceiling decomposition); none changes a headline claim. The R² intervals including negative values support the text's statement that structure-only models are not shown to beat a constant at 95% confidence on R².
- **Caveats.** Producer scripts were not re-run at L4; provenance is file-level. If a value is challenged, re-run the producing script listed in `outputs/deposit/manifest_repo_deposit.csv`.
- **Claims:** supports C-001-family benchmark statements; no new claim.
