# Kinase Inhibitor Selectivity: Towards a Formal Definition

## Project status
Under major revision for Molecular Informatics (resubmission). Empirical
analysis complete on four datasets. Manuscript revised in response to reviewer
comments (language/terminology, additional datasets, equation derivations).

## Research question
Existing kinase inhibitor selectivity definitions (S-score, entropy, Gini, ratio)
are used interchangeably but measure different things. This paper:
1. Systematically characterizes where and why they disagree
2. Identifies three mechanistic sources of instability
3. Proposes desiderata toward a formal axiomatic definition
4. Explores correlation with clinical outcomes (FAERS + RCT data)

## Key findings so far
- Ratio definition clusters separately from S-score/entropy/Gini in all four
  datasets (cross-family r=0.14-0.62 vs within-distribution-family r=0.74-0.99)
- Three instability sources identified:
  * Type 1: Zero-active drugs (no pKd>6 binding) — definitional noise, rank_std ~74
  * Type 2: Near-tied top targets (small top1_top2_gap) — ratio-specific instability (Klaeger r=-0.341***)
  * Type 3: High n_active with large gap — entropy-ratio disagreement
- n_active / active_range / active_std predict instability within entropy/Gini/S-score
  family (Klaeger n_active~entropy r=-0.455***)
- Replicated across four datasets / three assay technologies: Davis (68, pKd),
  Klaeger (222, pKd), Anastassiadis (178, % inhibition), Metz (704, pKi)
- Operationalization check: literal vs log-affinity entropy/Gini diverge
  (entropy r=0.42 Davis, -0.32 Klaeger) — Figure S1
- Clinical (FAERS + FDA-label) analysis: null, confounded — demoted to SI

## Datasets
- Davis 2011: davis_affinity.csv, davis_drugs.csv, davis_proteins.csv
- Klaeger 2017: aan4368_Table_S2.xlsx -> klaeger_matrix.csv
- Anastassiadis 2011: anastassiadis.xls -> anastassiadis_matrix.csv (% inhibition)
- Metz 2011: metz.xls -> metz_matrix.csv (pKi, >=50-kinase subset)
- NOT used: Gao 2013 (same modality as Anastassiadis, not redistributable),
  Patricelli 2011 (too few inhibitors for ranking analysis)
- FAERS / clinical: faers_counts.csv, selectivity_outcomes_merged.csv

## Key scripts
- selectivity_analysis.py  — main analysis on Davis
- instability_by_family decomposition — inline in session (needs extraction to script)
- Klaeger analysis — inline in session (needs extraction to script)
- faers_counts pull — inline in session

## Paper structure agreed
1. Abstract
2. Introduction (formal definition framing — Arrow, Shannon/Khinchin as precedents)
3. Related Work
   3.1 Existing selectivity metrics and comparisons
       - Karaman 2008 (S-score), Uitdehaag 2011 (entropy), Graczyk 2010 (partition index)
       - Bosc et al 2017 (WS/RS metrics, closest prior work — same datasets, fixed params)
       - Bajorath/Klaeger reconciliation 2018 (consistency finding, different approach)
       - KInhibition 2018 (noted conflict, practical not theoretical response)
   3.2 Axiomatic approaches to ranking/measurement
       - Shannon/Khinchin axioms (entropy uniqueness — the model for what we want)
       - Arrow impossibility (ranking axioms, possible impossibility result)
       - Truth discovery axiomatization (direct structural precedent, Pavillidis 2022)
4. Methods
   4.1 Datasets (Davis, Klaeger, Anastassiadis, Metz)
   4.2 Selectivity definitions and parameterization (original -> implemented -> meaning)
   4.3 Rank stability analysis methodology
5. Results
   5.1 Definitions cluster into two families
   5.2 Two-family structure replicates across assay technologies (4 datasets, Fig 1)
   5.3 Instability is structured and predictable
   5.4 Three-way taxonomy of instability sources
   5.5 Panel size dependence
6. Required Properties for a well-formed selectivity measure (formerly "Desiderata")
   D1: Reliability threshold (motivated by zero-active finding)
   D2: Bounded gap sensitivity (motivated by ratio instability)
   D3: Distributional consistency (motivated by entropy/Gini param sensitivity)
   D4: Monotonicity under weak off-target addition
7. Discussion (clinical/FAERS null result now in Limitations + SI)
8. Conclusion

## Target venue
Primary: Molecular Informatics (Wiley) -- under major revision
Preprint: ChemRxiv (https://doi.org/10.26434/chemrxiv.15001618/v1)

## TODO
- [ ] Extract Klaeger analysis to standalone script
- [ ] Get RCT grade 3/4 adverse event rates (DailyMed or manual)
- [ ] Correlate RCT data with selectivity scores
- [ ] Draft Section 1 (Introduction)
- [ ] Draft Section 3 (Related Work) — notes above are complete
- [ ] Draft Section 6 (Desiderata) — this is the novel conceptual contribution

## Clinical outcomes analysis (completed)
- FDA label extraction: 33/35 drugs with discontinuation rates
- FAERS: 46 drugs, all outcome correlations null (r<0.05, all p>0.14)
- RCT discontinuation rates: also null (r=-0.047 to +0.002, all p>0.14)
- Grade 3/4 rates: only 5 drugs with data — insufficient for analysis
- Conclusion: confounded by indication severity, prescription volume, in vitro/in vivo gap
- Framing: null result reported honestly, identifies rigorous future validation approach
- Files: faers_counts.csv, clinical_safety_data.csv, selectivity_outcomes_merged.csv

## TODO (updated)
- [x] Extract Klaeger analysis to standalone script (partial — inline scripts exist)
- [x] Get RCT grade 3/4 adverse event rates
- [x] Correlate RCT data with selectivity scores
- [ ] Draft Section 3 (Related Work)
- [ ] Draft Section 6 (Desiderata)
- [ ] Draft Section 2 (Introduction)
- [ ] Draft remaining sections
