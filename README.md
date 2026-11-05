# Kinase Inhibitor Selectivity: Definitional Instability and Required Properties

**Paper:** Definitional Instability in Kinase Inhibitor Selectivity: An Empirical
Characterization and Motivated Properties for a Well-Formed Measure
**Author:** Polina Vinogradova (Input Output Global, Singapore)
**Preprint:** https://doi.org/10.26434/chemrxiv.15001618/v1

---

## Part 1: Plain-language explanation

### Background: kinases and cancer drugs

Your body's cells constantly send signals — "grow," "divide," "die." Kinases are
the molecular switches that carry these signals: you have over 500 of them, and
they work by tagging other proteins with a phosphate group, switching those
proteins on or off. When a kinase malfunctions (usually from a mutation) it can
drive uncontrolled growth — cancer. A **kinase inhibitor** is a drug that blocks
a specific kinase. Imatinib (Gleevec) was the first famous example.

### The problem: kinases look almost identical

All 500+ kinases share a nearly identical region where drugs bind, so a drug
designed for one kinase very often hits several others too. Sometimes that is
harmless; sometimes it causes side effects. The property describing how
concentrated a drug's effect is on its intended target versus everything else is
called **selectivity**.

### How selectivity is measured

You test a drug against a panel of hundreds of kinases and record how tightly it
binds each one — a list of 300–500 numbers. Then you apply a formula to collapse
that list into a single selectivity score. The field uses at least four different
formulas, invented by different groups:

- **S-score** — what fraction of the panel does the drug hit above a cutoff?
- **Selectivity entropy** — borrows Shannon's information-theory entropy; low
  entropy means activity concentrated on few targets.
- **Gini coefficient** — borrows the economics inequality measure; high Gini
  means binding concentrated on few kinases.
- **Ratio** — how much more tightly does the drug bind its top target than its
  next-best one? It ignores everything else.

These are used interchangeably in the literature. Our question: does it matter
which one you use?

### What this paper found

It matters, in a structured way. The S-score, entropy, and Gini
all ask roughly the same question (how spread out is the binding?) and, with
entropy and Gini on a common log-affinity scale, agree more with one another
than with the ratio. The **ratio** asks a genuinely different question (how big is the gap to
the single nearest competitor?) and is the consistent outlier. So a drug that
hits one kinase strongly and twenty others moderately looks **selective** by
ratio but **promiscuous** by entropy — neither answer is wrong; they answer
different questions.

We confirmed this across **four datasets spanning three assay technologies**
(competition binding, chemoproteomics, and functional-activity assays), most
clearly on the two pK_d datasets. All analyses count one reference target per
human protein kinase: the released panels include mutant and domain constructs
of the same gene (Davis) and non-kinase proteins (Klaeger). We also found that
the disagreements are associated with the shape of a drug's binding profile,
and that ratio-based scores become unreliable on small panels. On small random
panels, the share of drugs that can be scored at all limits recovery of the
rankings more than the choice of formula.

### What we proposed

Three **motivated properties** for a selectivity formula:

1. A drug with no binding above the activity threshold on the profiled panel
   should not get a selectivity score (or should have it flagged) — its score
   would rest on sub-threshold values alone, which each formula treats by its
   own convention.
2. The ranking should not depend much on where, between the assay floor and the
   activity threshold, an analyst places the activity cutoff.
3. Adding a kinase the drug binds weakly (above the assay floor but below the
   activity threshold) should not make it look *more* selective.

No formula we examined satisfies property 2: moving the activity cutoff
reorders 8–13% of drug pairs on three datasets, so the cutoff must be fixed in
the analysis protocol. Property 1 is a reporting rule any formula can adopt.
Property 3 holds for simple formulas such as a count of active kinases, but not
for the S-score, the Gini, or the entropy with a low baseline. We give a
**candidate measure**, a gated, smooth-hinge entropy, that provably satisfies
property 3 while still weighting weak binding. It ranks these datasets almost
identically to the entropy, so on these data it is not an improvement on it.

### The bigger picture

This is a first step. In the 1940s Claude Shannon didn't just propose a formula
for measuring information — he showed it was the *only* one satisfying a few
reasonable requirements. The long-term goal is the analogous result for
selectivity. This paper characterizes the mess and states properties a solution
should have; finding properties that single out one measure is left for future
work.

---

## Part 2: Reproducing the results

### Requirements

Python 3.10 with:

    numpy pandas scipy matplotlib openpyxl xlrd requests

Install with conda:

    conda create -n selectivity python=3.10
    conda activate selectivity
    pip install numpy pandas scipy matplotlib openpyxl xlrd requests

On Apple Silicon, set this before running any script:

    export KMP_DUPLICATE_LIB_OK=TRUE

### Datasets

| Dataset | Readout | Size (cpd × kinase) | Source |
|---------|---------|---------------------|--------|
| Davis | pK_d (competition binding) | 68 × 433 released, 68 × 363 curated | github.com/dingyan20/Davis-Dataset-for-DTA-Prediction |
| Klaeger | pK_d (chemoproteomics) | 222 × 343 released, 222 × 219 curated | Science 2017, suppl. (doi.org/10.1126/science.aan4368) |
| Anastassiadis | % inhibition @0.5 µM (functional) | 178 × 300 released, 178 × 291 curated | Nat. Biotechnol. 2011, suppl. (doi.org/10.1038/nbt.2017) |
| Metz | pK_i | 704 × 172 | Nat. Chem. Biol. 2011, suppl. (doi.org/10.1038/nchembio.530) |

The Gao (Biochem. J. 2013) and Patricelli (Chem. Biol. 2011) datasets named by a
reviewer were not included: Gao is the same single-concentration functional
modality as Anastassiadis (already represented) and is not openly
redistributable, and Patricelli profiles only a handful of inhibitors — too few
for a ranking-stability analysis. See the manuscript Methods for details.

Raw supplementary files (downloaded programmatically):

- Anastassiadis: `https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.2017/MediaObjects/41587_2011_BFnbt2017_MOESM23_ESM.xls`
- Metz: `https://static-content.springer.com/esm/art%3A10.1038%2Fnchembio.530/MediaObjects/41589_2011_BFnchembio530_MOESM137_ESM.xls`

### Data files

| File | Description |
|------|-------------|
| davis_affinity_released.csv, davis_drugs.csv, davis_proteins.csv | Davis pK_d data as released (433 constructs) |
| davis_affinity.csv | Davis curated panel (one reference construct per human protein kinase; curate_panels.py) |
| aan4368_Table_S2.xlsx | Klaeger raw supplementary data |
| klaeger_matrix_released.csv / klaeger_matrix.csv | Klaeger drug × target pK_d matrix as released (343 targets) / curated (219 protein kinases) |
| anastassiadis.xls / anastassiadis_matrix_released.csv / anastassiadis_matrix.csv | Anastassiadis raw / processed (300 targets) / curated (291) % inhibition matrix |
| metz.xls / metz_matrix.csv | Metz raw / processed pK_i matrix (≥50-kinase coverage subset) |
| cross_dataset_summary.csv | Cross-dataset two-family + instability summary |
| selectivity_results.csv, klaeger_selectivity_results.csv | Per-compound analysis outputs |
| faers_counts.csv, clinical_safety_data.csv, selectivity_outcomes_merged.csv | Clinical-outcome analysis data (Limitations) |
| label_texts.json, extended_labels.json | Raw FDA drug label text underlying clinical_safety_data.csv |
| bosc_catds_comparison.csv | Bosc et al. window/ranking score and CATDS-proxy comparison (Discussion) |
| clinical_outcomes_correlation.csv | FAERS/discontinuation-rate vs. selectivity-rank correlations (Discussion, Limitations) |

**Clinical data provenance.** `faers_counts.csv` is pulled live from the openFDA
adverse-event API by `faers_pull.py` and is fully script-reproducible (though
report counts grow over time, so re-running will not reproduce the exact
frozen numbers in the paper). `clinical_safety_data.csv`, `label_texts.json`,
and `extended_labels.json` are **manually curated**: discontinuation rates and
grade 3/4 adverse-event rates were extracted by hand from FDA drug label text
(no script performs this extraction). `selectivity_outcomes_merged.csv` holds
the frozen FAERS counts and rates; its median-rank columns are superseded,
because `clinical_outcomes_correlation.py` takes each drug's current ranks from
`klaeger_selectivity_results.csv`. `bosc_catds_comparison.py` and `clinical_outcomes_correlation.py` are
the permanent, reproducible scripts for the two correlation analyses computed
from these files and cited in the Discussion/Limitations sections.

### Scripts

| Script | Description |
|--------|-------------|
| selectivity_analysis.py | Main Davis selectivity analysis (two-family correlations, instability) |
| klaeger_analysis.py | Klaeger analysis (reproduces the Table 1 Klaeger triangle and the zero-active instability) |
| additional_datasets_analysis.py | Anastassiadis + Metz; reproduces the four-dataset two-family result (Figure 2) and `cross_dataset_summary.csv` |
| instability_by_family.py | Instability predictors per definition family (Figure 3, Table 2) |
| binding_profiles.py | Most stable, most unstable, and most ratio-sensitive Davis profiles (Figure 4) |
| metric_fidelity_robustness.py | Operationalization-sensitivity check; original vs. log-affinity entropy/Gini (Figure 1) |
| panel_size_analysis.py | Panel-size subsampling on all four datasets (Figure 5) |
| candidate_measure.py | Candidate measure: orientation, agreement with the entropy, D2/D3 checks on Klaeger |
| faers_pull.py | Pulls FAERS adverse-event counts from the openFDA API (clinical analysis, Limitations) |
| bosc_catds_comparison.py | Reproduces Bosc et al. window/ranking score and CATDS-proxy correlations against the ratio family (Discussion, SI Table S2) |
| clinical_outcomes_correlation.py | Reproduces FAERS/discontinuation-rate vs. selectivity-rank correlations, using the frozen FAERS counts and the current Klaeger median ranks (Discussion, Limitations) |
| klaeger_filter_sensitivity.py | Klaeger high-confidence filter: entries removed, and sensitivity of Table 1 / panel-size p* to relaxing it (Methods, SI Section S6) |
| curate_panels.py | Builds the curated panels (one reference target per human protein kinase) from the released Davis, Klaeger and Anastassiadis files (Methods) |
| klaeger_protein_kinases.py | Klaeger panel composition (120 of 343 targets are not protein kinases) and the main Klaeger results on the released and curated panels (SI Table S6) |
| davis_gene_level.py | The main Davis results on the released constructs, the curated panel and the best affinity per gene (SI Table S6) |
| panel_recovery.py | Recovery from random subpanels for all measures, with coverage, pair recovery and an influence diagnostic (Results, Figure 7) |
| metz_coverage.py | Per-compound kinase-coverage distribution among the 704 retained Metz compounds (Methods) |
| dataset_descriptive_stats.py | Active-fraction, zero-active-threshold, and missing-entry percentages cited in Methods for all four datasets, plus the Davis stricter-cutoff (pKd>7) robustness check |
| d2_near_tie_sensitivity.py | The ratio's reference-choice sensitivity: Klaeger compounds with an active kinase and a top1-top2 gap < 0.05 and their R(1) vs. R(2) swing (Results, Type 2) |
| original_operationalization_robustness.py | Cross-definition correlations under the original (association-constant entropy, native %inhibition Gini) operationalizations (Methods, SI Table S3) |
| undefined_score_convention.py | States the entropy/Gini 0/0 convention for zero-active compounds and the alternative defined-configs-only rank-sigma contrast (Results, D1) |
| pstar_bar_sensitivity.py | Panel-size p* sensitivity to the r > 0.90 recoverability bar (Results, SI Table S1) |
| candidate_validation.py | Cross-dataset validation of the candidate measure: orientation, D2 (emphasis- and gate-swept), D3, smoothing-width sweep on Klaeger (Results, SI Table S4 and Section S5) |
| d3_pair_reversals.py | D2 (distributional consistency) for every measure on all four datasets: share of pairs reversing over B = [x0, tau*], exact for the S-score and on 201 grid points otherwise (Results, Table 3) |
| davis_klaeger_overlap.py | Resolves Davis/Klaeger compound overlap via PubChem CID matching (46 compounds; Discussion) |
| platform_vs_definition.py | Platform-vs-definition ranking-disagreement benchmark on the 46 shared Davis/Klaeger compounds (Discussion) |
| failure_mode_illustration.py | Figure 6: D1 (zero-active instability) and D2 (gate-swept vs. emphasis-swept instability, all over B = [5.0, 6.0]) |
| d4_entropy_vs_candidate.py | D3 (monotonicity) for the entropy vs. the candidate: the condition on T and its sharpness, the entropy counterexample, real-data frequency, and the effect of setting T = tau* - x0 (Methods proposition; Results, "Why the Smooth Hinge Matters") |
| feature_placeholder_sensitivity.py | Table 2 predictors recomputed without the 0 placeholder for undefined skewness/std/range (Results) |
| claim_checks.py | Numbers for specific claims: Type 1 target capture, Klaeger correlations without zero-active compounds and bootstrap, mean vs. median ranks, the Anastassiadis ceiling, Type 3 and stability patterns, Davis/Klaeger and Klaeger/Anastassiadis overlap, the Metz coverage dependence, the zero-active contrast per definition on every dataset, the Gini panel size with undefined compounds excluded on every dataset, active affinity range vs. instability with the number of active kinases held fixed, the weight near the threshold behind the skewness–Gini association, range vs. S-score instability on the 46 shared Davis/Klaeger compounds, where each definition ranks the most stable compounds on Davis and Klaeger, D3 under realistic panel expansion, where each definition ranks zero-active compounds, Table 2 with a Bonferroni correction, the candidate's sensitivity to the activity threshold, and its stability under measurement noise (SI Table S7 maps every result to its script and section) |
| instability_predictor_candidates.py | Five further candidate predictors of rank instability (CR1, HHI, top1 pKd, near-tie count, active-affinity CV) against Table 2's features (Results) |
| family_construction.py | Variants of the candidate: quantile floor (which coincides with the featured measure on these data), Rényi order; D1–D3 checks on Davis and Klaeger (SI Section S5, Table S5) |
| family_panelsize.py | Panel-size p* of the candidate by Rényi order on Klaeger and Metz (SI Section S5, Table S5) |
| family_alternatives.py | Two further candidate variants kept as negative results (soft gate, rank-based weights); not reported in the paper |
| script_logging.py | Shared helper imported by every script: mirrors printed output to `outputs/<script>.txt` and copies regenerated manuscript figures into `paper/figures/` and `paper/submission_figures/` |
| toc_graphic.py | Rebuilds the graphical-abstract/TOC image from the Davis correlation matrix and Klaeger panel-size curves |

Some script names and printed labels keep the property numbering of the
previous revision: D3 there is the paper's D2 (distributional consistency) and
D4 the paper's D3 (monotonicity). The paper's former D2 (bounded gap
sensitivity) was removed.

### Replication steps

**Preprocess Klaeger** (produces `klaeger_matrix_released.csv`, expected shape `(222, 343)`; then run `curate_panels.py` to write the curated panels every analysis reads):

    python3 -c "
    import pandas as pd, numpy as np
    df = pd.read_excel('aan4368_Table_S2.xlsx', sheet_name='Kinobeads')
    df = df[df['Target Classification'] == 'High confidence'].copy()
    df['pKd'] = 9 - np.log10(df['Apparent Kd'])
    matrix = df.groupby(['Drug','Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
    matrix.to_csv('klaeger_matrix_released.csv')
    print(matrix.shape)
    "

**Run the analyses:**

    python3 curate_panels.py                # curated panels (run first)
    python3 selectivity_analysis.py         # Davis
    python3 klaeger_analysis.py             # Klaeger (Table 1 values)
    python3 instability_by_family.py        # Figure 3, Table 2
    python3 additional_datasets_analysis.py # Anastassiadis + Metz (Figure 2, summary)
    python3 metric_fidelity_robustness.py   # Figure 1
    python3 panel_size_analysis.py          # Figure 5
    python3 candidate_measure.py            # candidate measure (orientation, agreement with the entropy)
    python3 panel_recovery.py               # recovery from random subpanels (Figure 7)
    python3 faers_pull.py                    # clinical data, Limitations (needs internet)
    python3 bosc_catds_comparison.py         # Bosc/CATDS comparison (Discussion)
    python3 clinical_outcomes_correlation.py # FAERS/discontinuation correlations (Discussion, Limitations)

`additional_datasets_analysis.py` reads `anastassiadis_matrix.csv` / `metz_matrix.csv`
if present, otherwise regenerates them from the raw `.xls` files.
`clinical_outcomes_correlation.py` requires `selectivity_outcomes_merged.csv`
and `clinical_safety_data.csv`, both manually curated (see Data files above)
and already committed — it does not need internet access.

### Expected key results

| Result | Value |
|--------|-------|
| Two-family: ratio vs. distribution family, all compounds | r = 0.51–0.56 (Davis), 0.09–0.66 (Klaeger), 0.48–0.55 (Anastassiadis), 0.14–0.19 (Metz) |
| Within-distribution-family correlation, compounds with an active kinase | r = 0.81–1.00 |
| Ratio vs. entropy correlation (Davis / Klaeger) | r = 0.513 / 0.344 |
| top1–top2 gap vs. ratio instability (Klaeger) | r = −0.333, p < 0.001 |
| n_active vs. entropy instability (Klaeger) | r = −0.369, p < 0.001 |
| Zero-active vs. active rank std (Klaeger) | 78.0 vs. 32.1 |
| Zero-active vs. active rank std (Anastassiadis / Metz) | 36.6 vs. 21.8 / 178.5 vs. 108.1 |
| Pairs reversing as the activity boundary moves (S-score, Davis/Klaeger/Anastassiadis) | 8.4–13.0% |
| Entropy panel-size stability threshold | 50 targets (Davis, Anastassiadis), 80 (Klaeger, Metz) |
| Ratio panel-size stability threshold | 350 of 363 (Davis), 219 of 219 (Klaeger), 260 of 291 (Anastassiadis), 170 of 172 (Metz) |
| Entropy operationalization sensitivity (literal vs. log-affinity) | r = 0.37 (Davis), −0.31 (Klaeger) |

### Repository structure

    .
    ├── README.md, RESEARCH_NOTES.md
    ├── selectivity_analysis.py, klaeger_analysis.py
    ├── additional_datasets_analysis.py, metric_fidelity_robustness.py
    ├── panel_size_analysis.py, candidate_measure.py, faers_pull.py
    ├── bosc_catds_comparison.py, clinical_outcomes_correlation.py
    ├── davis_*.csv, klaeger_matrix.csv, aan4368_Table_S2.xlsx
    ├── anastassiadis.xls, anastassiadis_matrix.csv
    ├── metz.xls, metz_matrix.csv
    ├── cross_dataset_summary.csv, selectivity_results.csv, klaeger_selectivity_results.csv
    ├── faers_counts.csv, clinical_safety_data.csv, selectivity_outcomes_merged.csv
    ├── label_texts.json, extended_labels.json (manually curated FDA label text)
    ├── bosc_catds_comparison.csv, clinical_outcomes_correlation.csv
    ├── *.png  (analysis figures)
    ├── outputs/  (printed results of every script, one .txt per script)
    └── paper/
        ├── main.tex, references.bib, main.pdf
        ├── supporting_information.tex, supporting_information.pdf
        ├── figures/  (metric_fidelity, cross_dataset_correlations, instability_by_family,
        │             binding_profiles, panel_size_stability, failure_modes,
        │             panel_recovery, toc_graphic)
        ├── submission_figures/, submission_tables/  (journal upload copies)
        └── sections/ (abstract, introduction, related_work, methods, results,
                       desiderata, discussion, conclusion)

### Citation

    Vinogradova, P. (2025). Definitional Instability in Kinase Inhibitor
    Selectivity: An Empirical Characterization Across Four Datasets and Required
    Properties for a Well-Formed Measure. ChemRxiv.
    https://doi.org/10.26434/chemrxiv.15001618/v1
