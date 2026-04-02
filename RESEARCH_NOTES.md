# Kinase Inhibitor Selectivity: Towards a Formal Definition

## Project status
Work in progress. Empirical analysis complete on two datasets.
Paper draft in progress.

## Research question
Existing kinase inhibitor selectivity definitions (S-score, entropy, Gini, ratio)
are used interchangeably but measure different things. This paper:
1. Systematically characterizes where and why they disagree
2. Identifies three mechanistic sources of instability
3. Proposes desiderata toward a formal axiomatic definition
4. Explores correlation with clinical outcomes (FAERS + RCT data)

## Key findings so far
- Ratio definition clusters separately from S-score/entropy/Gini (r=0.27-0.48 vs r=0.75-0.91)
- Three instability sources identified:
  * Type 1: Zero-active drugs (no pKd>6 binding) — definitional noise, rank_std ~74
  * Type 2: Near-tied top targets (small top1_top2_gap) — ratio-specific instability (r=-0.293***)
  * Type 3: High n_active with large gap — entropy-ratio disagreement
- active_range and active_std predict instability within entropy/Gini/S-score family
- Replicated on both Davis (68 drugs) and Klaeger (222 drugs) datasets
- FAERS analysis: null result (confounded by indication severity — noted as limitation)
- RCT adverse event data extraction: in progress (DailyMed/FDA labels)

## Datasets
- Davis 2011: davis_affinity.csv, davis_drugs.csv, davis_proteins.csv
- Klaeger 2017: aan4368_Table_S2.xlsx -> klaeger_matrix.csv
- FAERS counts: faers_counts.csv
- Merged outcomes: selectivity_outcomes_merged.csv

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
   4.1 Datasets (Davis, Klaeger)
   4.2 Selectivity definitions and parameterization
   4.3 Rank stability analysis methodology
   4.4 FAERS and RCT outcome data
5. Results
   5.1 Definitions cluster into two families
   5.2 Instability is structured and predictable
   5.3 Three-way taxonomy of instability sources
   5.4 Clinical outcome correlations (FAERS null + RCT attempt)
6. Desiderata for a well-formed selectivity measure
   D1: Reliability threshold (motivated by zero-active finding)
   D2: Bounded gap sensitivity (motivated by ratio instability)
   D3: Distributional consistency (motivated by entropy/Gini param sensitivity)
   D4: Monotonicity under weak off-target addition
7. Discussion
8. Conclusion

## Target venue
Primary: Journal of Chemical Information and Modeling (JCIM)
Preprint: ChemRxiv first

## TODO
- [ ] Extract Klaeger analysis to standalone script
- [ ] Get RCT grade 3/4 adverse event rates (DailyMed or manual)
- [ ] Correlate RCT data with selectivity scores
- [ ] Draft Section 1 (Introduction)
- [ ] Draft Section 3 (Related Work) — notes above are complete
- [ ] Draft Section 6 (Desiderata) — this is the novel conceptual contribution
