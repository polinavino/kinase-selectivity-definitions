# Kinase Inhibitor Selectivity: Definitional Instability and Required Properties

**Paper:** Definitional Instability in Kinase Inhibitor Selectivity: An Empirical
Characterization Across Four Datasets and Required Properties for a Well-Formed Measure
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

It matters, in a structured and predictable way. The S-score, entropy, and Gini
all ask roughly the same question (how spread out is the binding?) and mostly
agree. The **ratio** asks a genuinely different question (how big is the gap to
the single nearest competitor?) and is the consistent outlier. So a drug that
hits one kinase strongly and twenty others moderately looks **selective** by
ratio but **promiscuous** by entropy — neither answer is wrong; they answer
different questions.

We confirmed this across **four datasets spanning three assay technologies**
(competition binding, chemoproteomics, and functional-activity assays), so it is
not an artifact of one dataset or one assay. We also found that the disagreements
are predictable from the shape of a drug's binding profile, and that ratio-based
scores become unreliable on the small (50–100 kinase) panels used in practice,
while entropy-based scores stabilize above ~110 kinases.

### What we proposed

Four **required properties** that any good selectivity formula should satisfy:

1. A drug with no meaningful binding to any kinase should not get a selectivity
   score — it is just noise.
2. Small changes in which off-target you compare against should not cause large
   swings in the score.
3. The ranking should not change completely from a small change to a background
   parameter.
4. Adding a very weak accidental binding interaction should not make a drug look
   *more* selective.

No existing formula satisfies all four at once. We then propose a **candidate
measure** that does: the selectivity entropy with two fixes — a reliability gate
(property 1) and a smooth rather than hard activity cutoff (property 3). It
recovers its full-panel ranking from roughly half a typical panel, so it stays
practical.

### The bigger picture

This is a first step. In the 1940s Claude Shannon didn't just propose a formula
for measuring information — he showed it was the *only* one satisfying a few
reasonable requirements. The long-term goal is the analogous result for
selectivity. This paper characterizes the mess, states the properties a solution
should have, and offers a candidate that meets them; proving whether it is the
*unique* such measure is left for future work.

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
| Davis | pK_d (competition binding) | 68 × 433 | github.com/dingyan20/Davis-Dataset-for-DTA-Prediction |
| Klaeger | pK_d (chemoproteomics) | 222 × 343 | Science 2017, suppl. (doi.org/10.1126/science.aan4368) |
| Anastassiadis | % inhibition @0.5 µM (functional) | 178 × 300 | Nat. Biotechnol. 2011, suppl. (doi.org/10.1038/nbt.2017) |
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
| davis_affinity.csv, davis_drugs.csv, davis_proteins.csv | Davis pK_d data |
| aan4368_Table_S2.xlsx | Klaeger raw supplementary data |
| klaeger_matrix.csv | Processed Klaeger drug × kinase pK_d matrix |
| anastassiadis.xls / anastassiadis_matrix.csv | Anastassiadis raw / processed % inhibition matrix |
| metz.xls / metz_matrix.csv | Metz raw / processed pK_i matrix (≥50-kinase coverage subset) |
| cross_dataset_summary.csv | Cross-dataset two-family + instability summary |
| selectivity_results.csv, klaeger_selectivity_results.csv | Per-compound analysis outputs |
| faers_counts.csv, clinical_safety_data.csv, selectivity_outcomes_merged.csv | Clinical-outcome (SI) analysis data |

### Scripts

| Script | Description |
|--------|-------------|
| selectivity_analysis.py | Main Davis selectivity analysis (two-family correlations, instability) |
| klaeger_analysis.py | Klaeger analysis (reproduces Table 1 lower triangle and Table 2) |
| additional_datasets_analysis.py | Anastassiadis + Metz; reproduces the four-dataset two-family result (Figure 1) and `cross_dataset_summary.csv` |
| metric_fidelity_robustness.py | Operationalization-sensitivity check; original vs. log-affinity entropy/Gini (Figure S1) |
| panel_size_analysis.py | Panel-size subsampling on Klaeger (Figure 4) |
| candidate_measure.py | Candidate measure satisfying D1–D4; verifies D3/D4 and panel-size convergence (`candidate_panel_convergence.png`) |
| faers_pull.py | Pulls FAERS adverse-event counts from the openFDA API (clinical SI analysis) |

### Replication steps

**Preprocess Klaeger** (produces `klaeger_matrix.csv`, expected shape `(222, 343)`):

    python3 -c "
    import pandas as pd, numpy as np
    df = pd.read_excel('aan4368_Table_S2.xlsx', sheet_name='Kinobeads')
    df = df[df['Target Classification'] == 'High confidence'].copy()
    df['pKd'] = 9 - np.log10(df['Apparent Kd'])
    matrix = df.groupby(['Drug','Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
    matrix.to_csv('klaeger_matrix.csv')
    print(matrix.shape)
    "

**Run the analyses:**

    python3 selectivity_analysis.py         # Davis
    python3 klaeger_analysis.py             # Klaeger (Table 1 / Table 2 values)
    python3 additional_datasets_analysis.py # Anastassiadis + Metz (Figure 1, summary)
    python3 metric_fidelity_robustness.py   # Figure S1
    python3 panel_size_analysis.py          # Figure 4
    python3 candidate_measure.py            # candidate measure (Figure 5, D3/D4 checks)
    python3 faers_pull.py                    # clinical SI data (needs internet)

`additional_datasets_analysis.py` reads `anastassiadis_matrix.csv` / `metz_matrix.csv`
if present, otherwise regenerates them from the raw `.xls` files.

### Expected key results

| Result | Value |
|--------|-------|
| Two-family: ratio vs. distribution family | r = 0.34–0.45 (Davis), 0.27–0.62 (Klaeger), 0.52–0.56 (Anastassiadis), 0.14–0.19 (Metz) |
| Within-distribution-family correlation (all datasets) | r = 0.74–0.99 |
| Ratio vs. entropy correlation (Davis / Klaeger) | r = 0.343 / 0.480 |
| top1–top2 gap vs. ratio instability (Klaeger) | r = −0.341, p < 0.001 |
| n_active vs. entropy instability (Klaeger) | r = −0.455, p < 0.001 |
| Zero-active vs. active rank std (Klaeger) | ~73.8 vs. ~32.4 |
| Zero-active vs. active rank std (Anastassiadis / Metz) | 37.6 vs. 22.0 / 178.7 vs. 108.2 |
| Entropy panel-size stability threshold | ~110 kinases |
| Ratio panel-size stability threshold | > 320 kinases |
| Entropy operationalization sensitivity (literal vs. log-affinity) | r = 0.42 (Davis), −0.32 (Klaeger) |

### Repository structure

    .
    ├── README.md, RESEARCH_NOTES.md
    ├── selectivity_analysis.py, klaeger_analysis.py
    ├── additional_datasets_analysis.py, metric_fidelity_robustness.py
    ├── panel_size_analysis.py, candidate_measure.py, faers_pull.py
    ├── davis_*.csv, klaeger_matrix.csv, aan4368_Table_S2.xlsx
    ├── anastassiadis.xls, anastassiadis_matrix.csv
    ├── metz.xls, metz_matrix.csv
    ├── cross_dataset_summary.csv, *_selectivity_results.csv
    ├── faers_counts.csv, clinical_safety_data.csv, selectivity_outcomes_merged.csv
    ├── *.png  (analysis figures)
    └── paper/
        ├── main.tex, references.bib, main.pdf
        ├── figures/  (binding_profiles, instability_by_family, panel_size_stability,
        │             cross_dataset_correlations, metric_fidelity, toc_graphic)
        └── sections/ (abstract, introduction, related_work, methods, results,
                       desiderata, discussion, conclusion)

### Citation

    Vinogradova, P. (2025). Definitional Instability in Kinase Inhibitor
    Selectivity: An Empirical Characterization Across Four Datasets and Required
    Properties for a Well-Formed Measure. ChemRxiv.
    https://doi.org/10.26434/chemrxiv.15001618/v1
