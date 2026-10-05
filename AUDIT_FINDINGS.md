# Selectivity paper: audit findings

Recovered from the read-only audit workflow `selectivity-overclaim-audit` (run `wf_8a5134dc-4f1`, 2026-09-22).

## Provenance and what survives

That run audited the manuscript in seven parallel passes, raised **124 findings**, and fanned them out to one verification agent each. It was **aborted before completion**: 30 of 124 verifications finished, and the harness stored only the first 401 characters of every agent result.

The consequence for this list: **the finding texts themselves are not recoverable.** What survives is the finding identifier, the section it came from, and — for the 30 verified ones — the verdict plus the opening of the verifier's reasoning, which usually states the substance. Everything below is quoted from those stored fragments and is marked where it breaks off. Nothing here has been re-derived from the paper, so treat an unverified row as a lead, not as a defect.

| | count |
|---|---|
| Findings raised | 124 |
| Verified before abort | 30 |
| Confirmed real | 17 |
| Rejected by verifier | 13 |
| Never verified | 94 |

## Coverage by area

| Area | Findings | Verified |
|---|---|---|
| Discussion / clinical (`D`) | 34 | 0 |
| Abstract & Introduction (`F`) | 31 | 19 |
| Results (`R`) | 29 | 11 |
| Cross-consistency (`C`) | 17 | 0 |
| Candidate measure (`CAND`) | 3 | 0 |
| Framing (`FRAME`) | 2 | 0 |
| Miscellaneous (`MISC`) | 2 | 0 |
| Rényi / entropy (`RENYI`) | 2 | 0 |
| Family structure (`FAM`) | 1 | 0 |
| Independence (`INDEP`) | 1 | 0 |
| Setup / methods (`SETUP`) | 1 | 0 |
| Tables (`TAB`) | 1 | 0 |

## Verified findings

Verdict is the verification agent's, not the auditor's. Reasoning is truncated at 401 characters by the harness; the break is marked with `[…]`.

### `F01` — **CONFIRMED** (high confidence)
> CONFIRMED, but the defect sits in contribution (i), not in the quoted literature sentence, and the auditor's supporting rationale is partly overstated.

(1) Both sentences appear verbatim. introduction.tex:15-18 "The field employs at least four distinct quantitative measures: the S-score,\\cite{Karaman2008} selectivity entropy,\\cite{Uitdeha… […]

### `F02` — **CONFIRMED** (high confidence)
> I reproduced all three components of the alleged overclaim, and the paper's own body text contradicts the abstract on two of them.

(1) Sentence verified verbatim. /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/abstract.tex lines 11-12: "We derive four properties that any well-formed selectivity measure should satisfy, then s… […]

### `F03` — rejected (medium confidence)
> Quote confirmed verbatim at abstract.tex lines 12-13. All three stored numbers reproduce exactly. But the allegation does not hold.

(a) Not a missing caveat. desiderata.tex attaches the caveat inside the same paragraph as the D3 claim, with the number: "D3 holds because the candidate has no free activity baseline: the boundary between ac… […]

### `F04` — **CONFIRMED** (high confidence)
> STEP 1 — Sentence confirmed verbatim. /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/abstract.tex, lines 12-13: "Finally, we propose a measure that does, recovering its full-panel ranking from about a third of a panel." The same figure recurs at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/introduction.tex:88… […]

### `F05` — rejected (high confidence)
> STEP 1 — sentence confirmed. abstract.tex lines 12-13 read: "Finally, we propose a measure that does, recovering its full-panel ranking from about a third of a panel." (The bash `cat` rendering dropped articles; the Read of the file matches the quote.)

STEP 2 — evidence checked, and the auditor's central premise does not hold.
(a) The p*… […]

### `F06` — **CONFIRMED** (high confidence)
> REAL. Sentence confirmed verbatim in /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/abstract.tex lines 5-7: "Across four profiling datasets spanning three assay technologies (68-704 compounds each), we show these definitions cluster into two families measuring different properties, and that rank instability is predictable from … […]

### `F07` — **CONFIRMED** (high confidence)
> Sentence confirmed verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/abstract.tex lines 5-7. All three pieces of the auditor's evidence reproduce, and the scope mismatch is real.

(1) Scope. methods.tex:264-266 states outright: "Associations between rank instability and binding profile features were assessed using Sp… […]

### `F08` — rejected (high confidence)
> The sentence appears verbatim in abstract.tex (lines 8-10), but the allegation does not hold up.

(1) The auditor's primary technical evidence is factually wrong. instability_by_family.py does NOT pool ranks across definitions: analyze() builds four separate per-family arrays and takes .std(axis=0) within each sweep independently (lines 76-8… […]

### `F09` — rejected (high confidence)
> Not a real finding; the sentence is defensible as written, and the auditor's central claim is contradicted by the cited paper itself.

(1) Quote verified verbatim in introduction.tex (lines 27-30): "Bosc et al.\\cite{Bosc2017} showed that these definitions can produce substantially different compound rankings on the same dataset, and the KI… […]

### `F10` — **CONFIRMED** (medium confidence)
> The sentence appears verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/introduction.tex lines 37-39: "Go/no-go decisions worth hundreds of millions of dollars in development costs are made on the basis of selectivity assessments whose sensitivity to definitional choice is unknown." It also survives unchanged from ol… […]

### `F11` — rejected (high confidence)
> Defensible as written; the auditor's central logical objection is itself mistaken.

1. SENTENCE CONFIRMED verbatim at introduction.tex:53-55.

2. THE "THEREFORE" DOES FOLLOW. The auditor overlooked that the correlations are Spearman RANK correlations, not Pearson. Verified at additional_datasets_analysis.py:82: `corr = np.array([[spearma… […]

### `F12` — **CONFIRMED** (medium confidence)
> Mixed finding: the alleged overclaim does not hold, but one alleged false fact does, and I reproduced it independently.

(1) OVERCLAIM (the auditor's headline) — NOT real. The sentence enumerates exactly what the paper delivers: "characterizing how the competing measures actually behave, stating the properties a good measure should have, a… […]

### `F13` — rejected (high confidence)
> The title is quoted correctly and the response-letter items A.1/A.3 exist as cited, but the auditor's substantive case collapses on three independently checkable points, two of which are factual errors about the document history.

(1) THE PAPER DOES NOT ASSERT NECESSITY. The word "necessary" appears zero times in abstract.tex, introduction… […]

### `F14` — rejected (high confidence)
> The auditor's factual evidence checks out (no bootstrap/permutation/replicate-noise code in any of the 20 .py files; outputs/selectivity_analysis.txt gives "Mean rank std per drug: 8.48" and outputs/klaeger_analysis.txt gives "Mean rank std: 34.30"), but the inference drawn from it does not. Three reasons.

(1) No unsupported claim is ac… […]

### `F15` — **CONFIRMED** (high confidence)
> STEP 1 — Sentence confirmed verbatim. /Users/polina/Documents/BioInfStuff/selectivity/paper/graphical_abstract_text.txt ends with exactly: "We identify three instability sources and derive four properties that no existing selectivity measure satisfies, then propose one that does."

STEP 2 — Auditor's evidence independently reproduced:
(a) … […]

### `F16` — rejected (high confidence)
> The sentence appears verbatim in introduction.tex (lines 7-9), confirmed at byte level. But it is defensible as written, and the allegation does not survive scrutiny.

(1) The numeric claim is correct. Manning et al. 2002 enumerate 518 protein kinase genes, so "more than 500 proteins" is accurate, and 478 of 518 (92%) carry the canonical e… […]

### `F17` — **CONFIRMED** (high confidence)
> CONFIRMED, with one refinement to the auditor's framing.

(1) Sentence verified verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/introduction.tex lines 4-6: "Kinase inhibitors are among the most intensely investigated compound classes in drug discovery, with over 80 FDA-approved agents and several hundred in clinica… […]

### `F18` — rejected (high confidence)
> The sentence appears as quoted (abstract.tex lines 5-6). The allegation does not survive verification on any of its three legs.

1. THE ALLEGED INTERNAL INCONSISTENCY DOES NOT EXIST. The auditor argues Metz's readout "groups it with Davis and Klaeger rather than with Anastassiadis under the paper's own two-readout-type framing." But the pa… […]

### `F19` — rejected (medium confidence)
> The sentence appears as quoted (paper/sections/abstract.tex:2-4), but two of the auditor's three evidentiary legs are wrong, and the third is a hedging preference on a motivational premise rather than an unsupported conclusion.

(1) The Bosc attribution is not in the manuscript. The auditor says related_work.tex describes Bosc et al. as "… […]

### `R-01` — **CONFIRMED** (high confidence)
> I reproduced the problem independently.

(1) The sentence appears verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/results.tex lines 71-74: "Within the distribution-based family, the number of active kinases negatively predicts entropy instability in every dataset ($r = -0.28$ to $-0.43$), and the top$_1$--top$_2$ a… […]

### `R-02` — **CONFIRMED** (high confidence)
> Sentence verified verbatim at results.tex:285-287. "both datasets" means Davis + Klaeger in this section (line 9 vs lines 12/15; Metz and Anastassiadis are called "two further datasets" at line 260). Both cited examples are Klaeger-only: XL-228 (rank_std 7.0118) and PF-03814735 (7.8012) are the two most stable rows in klaeger_selectivity_re… […]

### `R-03` — **CONFIRMED** (high confidence)
> VERIFIED as a real false-fact plus unsupported mechanism.

(1) Sentence confirmed verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/results.tex lines 288-289: "When a compound binds nearly the entire profiled kinome, all definitions rank / it as non-selective regardless of parameterization, producing low instability.… […]

### `R-04` — **CONFIRMED** (high confidence)
> Reproduced independently and confirmed. (1) The sentence appears verbatim at results.tex:272-278 under the Type 3 heading. (2) outputs/klaeger_analysis.txt lines 33-44 match the auditor's evidence exactly: the three largest entropy-ratio disagreements are AZD-8055 (183.0, n_active=0), BMS-911543 (176.0, n_active=0) and Vatalanib (167.0, n_activ… […]

### `R-05` — rejected (high confidence)
> The sentence appears verbatim at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/results.tex:342-343, and the auditor's p* numbers reproduce exactly from outputs/panel_size_analysis.txt (Klaeger: entropy 110, gini 290, s_score 110, ratio 343; Metz: entropy 80, gini 80, s_score 50, ratio 170). But the numbers do not support the a… […]

### `R-06` — rejected (high confidence)
> Sentence confirmed verbatim at results.tex:313-318. The auditor transcribed outputs/pstar_bar_sensitivity.txt correctly but drew a false inference from it. On Metz at bar 0.90 the values are entropy 80, gini 80, s_score 50, ratio 170: ratio is slowest and the second-largest p* is 80, which is Gini's value, so Gini IS at the second-slowest posi… […]

### `R-08` — **CONFIRMED** (high confidence)
> Reproduced both components independently.

FALSE FACT (solid, high confidence). "with what would be obtained on a full kinome panel" mislabels the reference ranking. panel_size_analysis.py:60 sets ref_ranks = scores_at(M, ...) on the full data matrix and :66 subsamples with rng.choice(n_kinases, ps), so the comparison target is the dataset'… […]

### `R-09` — **CONFIRMED** (medium confidence)
> Quote verified verbatim at paper/sections/results.tex:65-70. The auditor's factual evidence checks out: cross_dataset_summary.csv gives zero_active = 0/16/44/5 and the sigma pairs exactly as printed; additional_datasets_analysis.py computes only group means (rank_std[...].mean()) with no test statistic; and grep over paper/sections/*.tex and … […]

### `R-10` — rejected (high confidence)
> The sentence appears verbatim at results.tex:253-255, and the auditor's provenance evidence is accurate: no script in the repo imports or computes any rank-sum/Mann-Whitney test (all scripts import only spearmanr, rankdata, kendalltau, skew), and no stored output contains the p-value; outputs/klaeger_analysis.txt and outputs/failure_mode_illus… […]

### `R-11` — **CONFIRMED** (high confidence)
> REPRODUCED, and the problem is slightly worse than alleged.

(1) The sentence appears verbatim in /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/results.tex, lines 215-218: "The remaining three were uniformly weaker predictors than the features already included, in every family and both datasets (full numbers in the code repo… […]

### `R-12` — **CONFIRMED** (high confidence)
> The sentence appears as quoted at /Users/polina/Documents/BioInfStuff/selectivity/paper/sections/results.tex lines 137-141. I verified Table 2 (lines 194-199) cell-by-cell against outputs/instability_by_family.txt: every value matches, so the table is the right authority to test the claim against.

CLAUSE (ii) IS A REAL FALSE FACT, and worse … […]

## Unverified findings

Raised by the audit pass but never reached verification. Identifiers only.

- **Discussion / clinical**: D-01, D-02, D-03, D-04, D-05, D-06, D-07, D-08, D-09, D-10, D-11, D-12, D-13, D-14, D-15, D-16, D-17, D-18, D-19, D-20, D-21, D-22, D-23, D1-01, D1-02, D1-03, D1-04, D2-01, D2-02, D3-01, D3-02, D3-03, D3-04, D3-05
- **Abstract & Introduction**: F20, F21, F22, F23, F24, F25, F26, F27, F28, F29, F30, F31
- **Results**: R-07, R-13, R-14, R-15, R-16, R-17, R-18, R-19, R-20, R-21, R-22, R-23, R-24, R-25, R-26, R-27, R-28, R-29
- **Cross-consistency**: C1, C10, C11, C12, C13, C14, C15, C16, C17, C2, C3, C4, C5, C6, C7, C8, C9
- **Candidate measure**: CAND-01, CAND-02, CAND-03
- **Framing**: FRAME-01, FRAME-02
- **Miscellaneous**: MISC-01, MISC-02
- **Rényi / entropy**: RENYI-01, RENYI-02
- **Family structure**: FAM-01
- **Independence**: INDEP-01
- **Setup / methods**: SETUP-01
- **Tables**: TAB-01

## Regenerating

The audit script is at `~/.claude/projects/-Users-polina-Documents-BioInfStuff-meta-bio-measure-framework-paper/8443f7c3-f60e-4cad-b69b-828263d7003a/workflows/scripts/selectivity-overclaim-audit-wf_8a5134dc-4f1.js`. Re-running it reproduces the audit pass but will not reproduce these identifiers, which were assigned per run.

