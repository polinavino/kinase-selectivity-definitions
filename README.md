# Kinase Inhibitor Selectivity: Towards a Formal Definition

**Paper:** Towards a Formal Definition of Kinase Inhibitor Selectivity: Empirical Characterization of Definitional Instability and Proposed Desiderata
**Author:** Polina Vinogradova, Input Output Global
**Preprint:** ChemRxiv (DOI to be added)

---

## Part 1: Plain-language explanation

### Background: kinases and cancer drugs

Your body's cells are constantly sending signals to each other — "grow," "divide," "die," "make this protein." Kinases are the molecular switches that carry these signals. You have over 500 different kinases, and they work by physically tagging other proteins with a small chemical label (a phosphate group), which switches those proteins on or off.

When kinases malfunction — usually because of a mutation — they can send the wrong signals and cause cells to grow uncontrollably. That's cancer.

A kinase inhibitor is a drug that blocks a specific kinase from working. The idea is: if a mutant kinase is causing cancer by sending "grow" signals all the time, block that kinase and the signal stops. Imatinib (Gleevec) was the first famous example — it blocks a kinase called BCR-ABL and transformed a previously fatal blood cancer into a manageable condition.

### The problem: kinases look almost identical

All 500+ kinases share a nearly identical region where drugs bind. This means that when you design a drug to block one kinase, it very often accidentally blocks several others too. Sometimes that is fine or even helpful. But sometimes blocking the wrong kinase causes side effects — heart problems, liver damage, immune suppression.

The property that describes how concentrated a drug's effect is on its intended target versus everything else is called **selectivity**.

### How selectivity is measured

You take your drug and test it against a panel of hundreds of kinases in a lab. For each kinase you measure how tightly the drug binds (the binding affinity). You end up with a list of 300-500 numbers, one per kinase. Then you apply a mathematical formula to that list to get a single selectivity score. A high score means the drug mostly hits one kinase. A low score means it hits many kinases roughly equally.

### The problem this paper found

The field uses at least four different formulas, each invented by different research groups:

- **S-score** counts how many kinases the drug hits above a cutoff threshold. Simple but sensitive to where you draw the line.
- **Selectivity entropy** borrows a formula from information theory — the same mathematics used to measure uncertainty in a phone signal — and applies it to the binding profile. Low entropy means activity concentrated on few targets.
- **Gini coefficient** borrows a formula from economics originally designed to measure income inequality. Applied here, it measures how unequally a drug distributes its binding across kinases.
- **Ratio** simply asks: how much more tightly does the drug bind its main target compared to its next-best target?

These four formulas are used interchangeably in the scientific literature. Our question was: does it matter which formula you use?

It matters a lot — but in a structured, predictable way.

The S-score, entropy, and Gini are all asking roughly the same question: how spread out is the drug's binding activity across the whole kinase panel? They mostly agree on rankings. The ratio definition asks a completely different question: how much better is the drug at hitting its top target compared to its second-best target? It ignores everything else entirely.

This means a drug that hits one kinase very strongly and twenty others moderately will look highly selective by ratio (big gap between first and second) but moderately promiscuous by entropy (activity spread across many targets). Neither answer is wrong — they are answering different questions.

We also found that if you only test a drug against 50 kinases instead of 300, ratio-based scores become essentially meaningless, while entropy-based scores remain reasonably stable above about 110 kinases.

### What we proposed

We proposed four properties that any good selectivity formula should have — called desiderata:

1. A drug with no meaningful binding to any kinase should not get a selectivity score — it is just noise
2. Small changes in which off-target you compare against should not cause large swings in the score
3. The ranking should not change completely just because you changed a background parameter slightly
4. Adding a very weak accidental binding interaction should not make a drug look more selective

No existing formula satisfies all four properties simultaneously.

### The bigger picture

This paper is the first step toward something more ambitious. In the 1940s, Claude Shannon did not just propose a formula for measuring information — he proved mathematically that his formula was the only one satisfying a small set of reasonable requirements. We want to do something similar for selectivity: find the formula (or prove no perfect formula exists) by starting from first principles. That is a harder problem left for future work. This paper characterizes the mess before cleaning it up.

---

## Part 2: Reproducing the results

### Requirements

Python 3.10 with the following packages:

    numpy pandas scipy matplotlib openpyxl requests

Install with conda:

    conda create -n selectivity python=3.10
    conda activate selectivity
    pip install numpy pandas scipy matplotlib openpyxl requests

On Apple Silicon, set this environment variable before running any script:

    export KMP_DUPLICATE_LIB_OK=TRUE

### Data files

| File | Description | Source |
|------|-------------|--------|
| davis_affinity.csv | Drug-kinase binding affinities (pKd), long format | github.com/dingyan20/Davis-Dataset-for-DTA-Prediction |
| davis_drugs.csv | Drug names and SMILES | Same |
| davis_proteins.csv | Kinase names and sequences | Same |
| aan4368_Table_S2.xlsx | Klaeger et al. 2017 raw data | science.org/doi/10.1126/science.aan4368 (supplementary) |
| klaeger_matrix.csv | Processed Klaeger drug-kinase matrix | Generated by preprocessing step |
| faers_counts.csv | FAERS adverse event counts | Generated by faers_pull.py |
| clinical_safety_data.csv | FDA label discontinuation rates | Manually extracted from FDA labels |
| selectivity_results.csv | Davis analysis output | Generated by selectivity_analysis.py |
| klaeger_selectivity_results.csv | Klaeger analysis output | Generated by klaeger_analysis.py |
| selectivity_outcomes_merged.csv | Selectivity scores merged with FAERS | Generated by analysis |

### Scripts

| Script | Description |
|--------|-------------|
| selectivity_analysis.py | Main Davis selectivity analysis |
| klaeger_analysis.py | Klaeger selectivity analysis |
| faers_pull.py | Pulls adverse event counts from openFDA API |
| panel_size_analysis.py | Subsampling analysis of panel size dependence |

### Replication steps

**Step 1: preprocess the Klaeger data**

Run the following to convert the raw Excel file to a drug x kinase pKd matrix:

    python3 -c "
    import pandas as pd, numpy as np
    df = pd.read_excel('aan4368_Table_S2.xlsx', sheet_name='Kinobeads')
    df = df[df['Target Classification'] == 'High confidence'].copy()
    df['pKd'] = 9 - np.log10(df['Apparent Kd'])
    matrix = df.groupby(['Drug','Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
    matrix.to_csv('klaeger_matrix.csv')
    print(matrix.shape)
    "

Expected output: (222, 343)

**Step 2: Davis selectivity analysis**

    python3 selectivity_analysis.py

Produces: selectivity_results.csv, selectivity_analysis.png, binding_profiles.png, instability_by_family.png

**Step 3: Klaeger selectivity analysis**

    python3 klaeger_analysis.py

Produces: klaeger_selectivity_results.csv and printed correlation/instability tables

**Step 4: FAERS adverse event data**

Requires internet access. Queries the openFDA API.

    python3 faers_pull.py

Produces: faers_counts.csv

**Step 5: panel size subsampling analysis**

    python3 panel_size_analysis.py

Produces: panel_size_stability.png and printed stability table

### Expected key results

| Result | Value |
|--------|-------|
| Ratio vs entropy correlation (Davis) | r = 0.343 |
| Ratio vs entropy correlation (Klaeger) | r = 0.480 |
| top1_top2_gap vs ratio instability (Klaeger) | r = -0.293, p < 0.001 |
| n_active vs entropy instability (Klaeger) | r = -0.474, p < 0.001 |
| Zero-active drug rank std | ~74 |
| Active drug rank std | ~32 |
| Entropy stability threshold | ~110 kinases |
| Ratio stability threshold | >320 kinases |

### Repository structure

    .
    +-- README.md
    +-- RESEARCH_NOTES.md
    +-- selectivity_analysis.py
    +-- klaeger_analysis.py
    +-- faers_pull.py
    +-- panel_size_analysis.py
    +-- davis_affinity.csv
    +-- davis_drugs.csv
    +-- davis_proteins.csv
    +-- aan4368_Table_S2.xlsx
    +-- klaeger_matrix.csv
    +-- klaeger_selectivity_results.csv
    +-- selectivity_results.csv
    +-- faers_counts.csv
    +-- clinical_safety_data.csv
    +-- selectivity_outcomes_merged.csv
    +-- *.png
    +-- paper/
        +-- main.tex
        +-- references.bib
        +-- sections/

### Citation

    Vinogradova, P. (2025). Towards a Formal Definition of Kinase Inhibitor
    Selectivity: Empirical Characterization of Definitional Instability and
    Proposed Desiderata. ChemRxiv. DOI: [to be added]
