# Analysis Choices and Ranking Stability in Kinase Inhibitor Selectivity Scores

**Paper:** Analysis Choices and Ranking Stability in Kinase Inhibitor
Selectivity Scores (submitted revision; the preprint below has an earlier
title)
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

It matters, and the answer depends on analysis choices that are usually left
implicit: how affinities are weighted, which drugs are ranked, how pairs with
no detected binding are treated, which targets count as kinases, and where the
activity cutoff sits. With entropy and Gini on a common log-affinity scale, and
among drugs that bind some kinase above the activity threshold, the S-score,
entropy and Gini agree more with one another than with the **ratio**, which
asks a different question (how big is the gap to the nearest competitor?). A
drug that hits one kinase strongly and twenty others moderately looks
**selective** by ratio but **promiscuous** by entropy. The entropy as
originally published, with association-constant weights, does not group with
the S-score and Gini.

We checked this on **four datasets spanning three assay technologies**
(competition binding, chemoproteomics, and functional-activity assays). All
analyses count one reference target per human protein kinase: the released
panels include mutant and domain constructs of the same gene (Davis) and
proteins that are not kinases (Klaeger). On small random panels some drugs can
no longer be scored, and formulas that follow the share of targets with
detected binding recover the rankings from fewer targets than the ratio or the
entropy as originally published.

### What we proposed

Four **criteria** for a selectivity formula:

1. Activity eligibility: a drug with no binding above the activity threshold on
   the profiled panel should not get a selectivity score (or should have it
   flagged), since its score would rest on sub-threshold values alone.
2. Locality: moving the activity cutoff should change a drug's score only
   when the cutoff passes one of that drug's measured values.
3. Adding a kinase the drug binds weakly (above the floor but below the
   activity threshold) should not make it look *more* selective.
4. Orientation: on the same panel, binding equally strongly to more kinases
   should not make a drug look more selective, and binding one kinase should
   score as more selective than binding all of them equally.

Moving the activity cutoff reorders 8–13% of drug pairs on three datasets
for the S-score, entropy and Gini, so the cutoff must be fixed in the analysis
protocol.
Criterion 1 is a reporting rule any formula can adopt. Criterion 2 holds for
the S-score and the count of active kinases, but not for the entropy, the Gini
or the ratio, which subtract or substitute the cutoff value. Criterion 3 holds
for simple formulas such as a count of active kinases or the S-score with its
cutoff at the floor, but not for the S-score with a higher cutoff, the Gini,
or the entropy with a low baseline. We derive a **gated entropy**, a Shannon
entropy with exponential affinity weights. Among formulas that weight each
kinase by its distance above the cutoff and satisfy Shannon's grouping rule,
criteria 2 and 4 leave only these weights, and criterion 3 holds exactly when
their scale (the affinity difference over which a weight grows by a factor of
e) is at least the distance between the floor and the activity
threshold. On these data it ranks drugs similarly to the share of targets with
detected binding but recovers rankings from small panels more slowly, so it is
not an empirical improvement on simpler formulas. With twice the scale it
recovers them as fast and ranks drugs almost as that share does, so the scale
should be reported. Every formula studied
satisfies criterion 4 when pairs without detected binding get zero weight; the
active count with its sign reversed, which meets criteria 1 to 3 as the gated
entropy does, fails it.

### The bigger picture

In the 1940s Claude Shannon showed that his formula for information was the
only one satisfying a few reasonable requirements. This paper states
requirements for selectivity formulas and characterizes how the existing
formulas disagree; it does not single out one formula.

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
| Klaeger | pK_d (chemoproteomics) | 222 × 343 released, 222 × 217 curated | Science 2017, suppl. (doi.org/10.1126/science.aan4368) |
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
| target_curation_davis.csv / target_curation_klaeger.csv / target_curation_anastassiadis.csv | Every released target, whether the curated panel keeps it, and why (curate_panels.py) |
| aan4368_Table_S2.xlsx | Klaeger raw supplementary data |
| klaeger_matrix_released.csv / klaeger_matrix.csv | Klaeger drug × target pK_d matrix as released (343 targets) / curated (217 targets the source annotates as protein kinases) |
| anastassiadis.xls / anastassiadis_matrix_released.csv / anastassiadis_matrix.csv | Anastassiadis raw / processed (300 targets) / curated (291) % inhibition matrix |
| metz.xls / metz_matrix.csv | Metz raw / processed pK_i matrix (≥50-kinase coverage subset) |
| cross_dataset_summary.csv | Cross-dataset correlation and instability summary |
| selectivity_results.csv, klaeger_selectivity_results.csv | Per-compound analysis outputs |
| faers_counts.csv, selectivity_outcomes_merged.csv | Clinical-outcome analysis data (Limitations) |
| clinical_safety_data.csv, label_texts.json, extended_labels.json | Hand-extracted FDA label discontinuation rates and the label text; no longer used, because an audit found extraction errors |
| bosc_catds_comparison.csv | Bosc et al. window/ranking score and CATDS-proxy comparison (SI Section S6) |
| clinical_outcomes_correlation.csv | FAERS vs. selectivity-rank correlations (Discussion, Limitations) |

**Clinical data provenance.** `faers_counts.csv` is pulled live from the openFDA
adverse-event API by `faers_pull.py` and is fully script-reproducible (though
report counts grow over time, so re-running will not reproduce the exact
frozen numbers in the paper). `selectivity_outcomes_merged.csv` holds
the frozen FAERS counts and rates; its median-rank columns are superseded,
because `clinical_outcomes_correlation.py` takes each drug's current ranks from
`klaeger_selectivity_results.csv`. `bosc_catds_comparison.py` and `clinical_outcomes_correlation.py` are
the permanent, reproducible scripts for the two correlation analyses computed
from these files and cited in the Discussion/Limitations sections.

### Scripts

| Script | Description |
|--------|-------------|
| selectivity_analysis.py | Main Davis selectivity analysis (correlations between definitions, instability) |
| klaeger_analysis.py | Klaeger analysis (Klaeger correlations of Table 3 and the zero-active instability) |
| additional_datasets_analysis.py | Anastassiadis + Metz; correlations on four datasets, all and active compounds (SI Figure S1), and `cross_dataset_summary.csv` |
| instability_by_family.py | Profile features against instability, active compounds ranked among themselves, with the all-compound ranking as a check (SI Table S4, Figure S3) |
| binding_profiles.py | Most ratio-sensitive and most unstable Davis profiles (Figure 2) |
| metric_fidelity_robustness.py | Published vs. log-affinity entropy/Gini, with floor-coded pairs at the floor or given zero weight (Figure 1) |
| panel_size_analysis.py | Panel-size subsampling on all four datasets without the D1 gate (SI Section S5, Figure S4) |
| candidate_measure.py | Candidate measure: correlation with the number of active kinases, agreement with the entropy, gate sweep and D3 checks on Klaeger |
| faers_pull.py | Pulls FAERS adverse-event counts from the openFDA API (clinical analysis, Limitations) |
| bosc_catds_comparison.py | Reproduces Bosc et al. window/ranking score and CATDS-proxy correlations against the ratio family (SI Section S6, Table S8) |
| clinical_outcomes_correlation.py | Reproduces FAERS vs. selectivity-rank correlations, using the frozen FAERS counts and the current Klaeger median ranks (Discussion, Limitations) |
| klaeger_filter_sensitivity.py | Klaeger high-confidence filter: entries removed, and sensitivity of Table 3 and of the gated and ungated p* to relaxing it (SI Section S1) |
| curate_panels.py | Builds the curated panels (one reference target per human protein kinase) from the released Davis, Klaeger and Anastassiadis files, and the target_curation_*.csv tables (Methods) |
| klaeger_protein_kinases.py | Klaeger panel composition (122 of 343 targets removed besides four protein groups) and the main Klaeger results on the released and curated panels (SI Table S1) |
| davis_gene_level.py | The main Davis results on the released constructs, the curated panel and the best affinity per gene (SI Table S1) |
| panel_recovery.py | Recovery from random subpanels for all measures, including the S-score at the floor and the published entropy, with coverage, pair recovery, p* at other bars, the probability that a compound stays scored, an influence diagnostic, and the gated entropy at twice and four times its scale (Results, Figure 3, "Scale of the Gated Entropy and Recovery") |
| metz_coverage.py | Per-compound kinase-coverage distribution among the 704 retained Metz compounds (Methods) |
| dataset_descriptive_stats.py | Active-fraction, zero-active-threshold, and missing-entry percentages cited in Methods for all four datasets, plus the Davis check with the S-score sweep shifted to [6.5, 9.0] (SI Section S1) |
| d2_near_tie_sensitivity.py | The ratio's reference-choice sensitivity: Klaeger compounds with an active kinase and a top1-top2 gap < 0.05 and their R(1) vs. R(2) swing (Results, Parameter Sensitivity) |
| original_operationalization_robustness.py | Cross-definition correlations under the original (association-constant entropy, native %inhibition Gini) operationalizations, two floor treatments, and on Klaeger every measured affinity kept at its value (Results, Table 3, SI Table S2) |
| undefined_score_convention.py | States the entropy/Gini 0/0 convention for zero-active compounds and the alternative defined-configs-only rank-sigma contrast (Results; SI Section S3) |
| pstar_bar_sensitivity.py | Ungated panel-size p* at four correlation bars (SI Table S7) |
| candidate_validation.py | Cross-dataset validation of the candidate measure: correlation with the number of active kinases, pairs reversing when the gate or the entropy's baseline is swept, D3, scale sweep on Klaeger (Results, SI Table S5 and Section S4) |
| d3_pair_reversals.py | Boundary sensitivity for every measure on all four datasets: share of pairs reversing over B = [x0, tau*], exact for the S-score; the gated entropy's gate also at every observed value, where its score jumps; 201 grid points otherwise (Results, SI Table S3) |
| davis_klaeger_overlap.py | Resolves Davis/Klaeger compound overlap via PubChem CID matching (46 compounds; Discussion) |
| platform_vs_definition.py | Platform-vs-definition ranking-disagreement benchmark on the 46 shared Davis/Klaeger compounds and 201 shared kinases, renamed symbols mapped (Results; SI Section S6) |
| failure_mode_illustration.py | SI Figure S2: D1 (zero-active instability) and instability when the boundary moves (entropy baseline vs. gated entropy gate, over B = [5.0, 6.0]) |
| d4_entropy_vs_candidate.py | D3 (monotonicity) for the entropy, the published entropy and the gated entropy: the condition on T and its sharpness, the entropy fact behind Proposition 1, counterexamples, real-data frequency, the switch from softplus weights, the jump at the floor and the strengthening example; D4 (orientation) on reference profiles, Proposition 3 and D2 (locality) (Methods; Results, "Criteria and the Gated Entropy"; SI Section S4) |
| feature_placeholder_sensitivity.py | Profile features recomputed without the 0 placeholder for undefined skewness/std/range (SI Section S3) |
| claim_checks.py | Numbers for specific claims: target capture for compounds with no detected binding, Klaeger correlations without zero-active compounds and bootstrap, mean vs. median ranks, the Anastassiadis ceiling, broad profiles and stability patterns, Davis/Klaeger overlap, the Metz coverage dependence, the zero-active contrast per definition on every dataset, the Gini panel size with undefined compounds excluded on every dataset, active affinity range vs. instability with the number of active kinases held fixed, the weight near the threshold behind the skewness–Gini association, range vs. S-score instability on the 46 shared Davis/Klaeger compounds, where each definition ranks the most stable compounds on Davis and Klaeger, D3 under realistic panel expansion, where each definition ranks zero-active compounds, SI Table S4 with a Bonferroni correction, the gated entropy's sensitivity to the activity threshold, and its stability under measurement noise (SI Table S9 maps every result to its script and section) |
| instability_predictor_candidates.py | Five further candidate predictors of rank instability (CR1, HHI, top1 pKd, near-tie count, active-affinity CV) against the features of SI Table S4 (SI Section S3) |
| family_construction.py | Variants of the candidate: quantile floor (which coincides with the featured measure on these data), Rényi order, equal and association-constant weights; orientation, D3 and agreement on Davis and Klaeger (SI Section S4, Table S6) |
| family_panelsize.py | Panel-size p* of the candidate by Rényi order on Klaeger and Metz (SI Section S4, Table S6) |
| family_alternatives.py | Two variants of the earlier softplus candidate kept as negative results (soft gate, rank-based weights); not reported in the paper |
| script_logging.py | Shared helper imported by every script: mirrors printed output to `outputs/<script>.txt` and copies regenerated manuscript figures into `paper/figures/` and `paper/submission_figures/` |
| toc_graphic.py | Rebuilds the graphical-abstract/TOC image from the Davis correlation matrix and Davis recovery curves |

Some script names keep the property numbering of earlier revisions:
`d3_pair_reversals.py` measures boundary sensitivity (the former D3,
distributional consistency, now reported as a measurement), and
`d4_entropy_vs_candidate.py` checks the paper's D3 (monotonicity). It also
checks D4 (orientation), Proposition 3, D4 with floor-coded pairs kept at the
floor, D4 for the entropy and Gini at baselines below the floor (sections 9 to
12 of its output), and D2 (locality, section 13).

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
    python3 klaeger_analysis.py             # Klaeger (Table 3 values)
    python3 instability_by_family.py        # SI Figure S3, Table S4
    python3 additional_datasets_analysis.py # Anastassiadis + Metz (SI Figure S1, summary)
    python3 metric_fidelity_robustness.py   # Figure 1
    python3 panel_size_analysis.py          # SI Figure S4
    python3 candidate_measure.py            # candidate measure (correlation with active count, agreement with the entropy)
    python3 panel_recovery.py               # recovery from random subpanels (Figure 3)
    python3 binding_profiles.py             # Figure 2
    python3 failure_mode_illustration.py    # SI Figure S2
    python3 d3_pair_reversals.py            # boundary sensitivity (SI Table S3)
    python3 claim_checks.py                  # numbers for specific claims (run last)
    python3 faers_pull.py                    # clinical data, Limitations (needs internet)
    python3 bosc_catds_comparison.py         # Bosc/CATDS comparison (SI Section S6)
    python3 clinical_outcomes_correlation.py # FAERS correlations (Discussion, Limitations)

Every other script in the table above runs the same way; SI Table S9 maps each reported result to its script.
`additional_datasets_analysis.py` reads `anastassiadis_matrix.csv` / `metz_matrix.csv`
if present, otherwise regenerates them from the raw `.xls` files.
`clinical_outcomes_correlation.py` requires `selectivity_outcomes_merged.csv`
(see Data files above), already committed — it does not need internet access.

### Expected key results

| Result | Value |
|--------|-------|
| Two-family: ratio vs. distribution family, all compounds | r = 0.51–0.56 (Davis), 0.08–0.66 (Klaeger), 0.48–0.55 (Anastassiadis), 0.14–0.19 (Metz) |
| Within-distribution-family correlation, compounds with an active kinase | r = 0.81–0.999 |
| Ratio vs. entropy correlation (Davis / Klaeger) | r = 0.513 / 0.344 |
| top1–top2 gap vs. ratio instability (Klaeger, active compounds ranked among themselves) | r = −0.275, p < 0.001 |
| n_active vs. entropy instability (Klaeger, active compounds ranked among themselves) | r = −0.312, p < 0.001 |
| Zero-active vs. active rank std (Klaeger) | 77.9 vs. 32.2 |
| Zero-active vs. active rank std (Anastassiadis / Metz) | 36.6 vs. 21.8 / 178.5 vs. 108.1 |
| Pairs reversing as the activity boundary moves (S-score, Davis/Klaeger/Anastassiadis) | 8.4–13.0% |
| Recovery with D1 applied, targets to mean Spearman 0.90 (Davis/Klaeger/Anastassiadis/Metz) | S-score at floor, entropy at floor, Gini: 20/50/30/65; gated entropy: 50/65/30/80 (20/50/30/65 with twice the scale); ratio: >300/200/>250/>150; published entropy: 300/200/–/150; P_max: 300/200/–/>150 |
| Panel-size threshold without the D1 gate, entropy (SI Section S5) | 50 targets (Davis, Anastassiadis), 80 (Klaeger, Metz) |
| Published vs. log-affinity entropy (Davis/Klaeger) | r = 0.37/−0.31 with floor-coded pairs at the floor, 0.74/0.79 with them given zero weight |

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
    ├── faers_counts.csv, selectivity_outcomes_merged.csv
    ├── clinical_safety_data.csv, label_texts.json, extended_labels.json (not used)
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
        └── sections/ (abstract, introduction, methods, results, discussion,
                       conclusion)

### Citation

    Vinogradova, P. (2025). Definitional Instability in Kinase Inhibitor
    Selectivity: An Empirical Characterization Across Four Datasets and Required
    Properties for a Well-Formed Measure. ChemRxiv.
    https://doi.org/10.26434/chemrxiv.15001618/v1
