# Kinase Inhibitor Selectivity: Research Notes

## Project status
Second revision for Molecular Informatics (manuscript 7756455), responding to
Reviewer 1's round-2 review (`paper/reviewer1_round2_review.txt`; response in
`paper/response_to_reviewers_round2.txt`). Analysis complete on four datasets;
every number in the paper is printed by a script into `outputs/` (see README).
Preprint: ChemRxiv (https://doi.org/10.26434/chemrxiv.15001618/v1).

## Research question
Existing kinase inhibitor selectivity definitions (S-score, entropy, Gini,
ratio) are used interchangeably. The paper:
1. Characterizes where and why they disagree, on four datasets curated to one
   reference target per human protein kinase (curate_panels.py)
2. Identifies three sources of rank instability
3. States three properties (D1 reliability threshold, D2 distributional
   consistency, D3 monotonicity under weak off-target addition) and evaluates
   the four definitions, two comparators and a candidate against them
4. Evaluates recovery from random subpanels separately (panel_recovery.py)

## Key findings
- Two families, with entropy and Gini on a common baseline-subtracted
  log-affinity scale: among compounds with an active kinase, S-score/entropy/
  Gini agree (r = 0.81-0.999) and the ratio agrees less with each. With all
  compounds, on Klaeger the Gini sits closer to the ratio (0.66 vs 0.51; not
  significant); on Anastassiadis the split comes mainly from the 45 zero-active
  compounds. Under the published association-constant entropy the grouping
  does not hold on any affinity dataset (SI Table S3).
- Instability sources (Klaeger): Type 1, zero-active compounds (22 of 222; all
  20 most unstable; pooled rank sigma 78.0 vs 32.1, from the definitions'
  conventions for such profiles); among the 20 most unstable active
  compounds, 6 Type 2 (near-tied top targets) and 8 Type 3 (>= 10 active).
- Replicated relations: ratio-gap on all four datasets (Anastassiadis only
  after excluding 8 ceiling-saturated compounds); entropy-n_active on all four
  (r = -0.30 to -0.43; Davis p = 0.01, not significant after Bonferroni).
- Properties: no measure satisfies D2 (the S-score's threshold reverses
  8.4-13.0% of pairs on three datasets; entropy/Gini 7.9-12.5%; candidate gate
  9.2-14.2%; ratio floor fallback up to 13.5%), so the activity boundary must
  be fixed by protocol. D3 holds for the active count, the entropy with
  beta = tau*, the ratio and the candidate; not for the S-score, Gini or the
  entropy with beta < tau*. The former D2 (bounded gap sensitivity) was
  dropped: different reference off-targets answer different questions.
- Candidate (gated smooth-hinge entropy): D1 and D3 (proved for additions in
  (x0, tau*]) while weighting binding below tau*; ranks almost identically to
  the entropy (Spearman 0.994-0.998); gate discontinuity costs little under
  measurement noise (Spearman >= 0.975).
- Recovery: weighting weak binding (entropy beta = x0, Gini, candidate)
  reaches Spearman 0.9 on co-scored compounds from fewer targets on Davis,
  Klaeger, Anastassiadis (20/50/30 vs 30/100/50 for the S-score/count); reversed
  on Metz. Coverage limits pair recovery; S-score and count identical.
- Released panels (SI S6): Davis 26 of 68 compounds have same-gene top two
  constructs; Klaeger 120 of 343 targets are not protein kinases.
- Operationalization check: literal vs log-affinity entropy r = 0.37 (Davis),
  -0.31 (Klaeger) (Figure 1).

## Datasets
- Davis 2011: davis_affinity_released.csv (68 x 433) -> davis_affinity.csv (curated, 363)
- Klaeger 2017: aan4368_Table_S2.xlsx -> klaeger_matrix_released.csv (222 x 343) -> klaeger_matrix.csv (curated, 219)
- Anastassiadis 2011: anastassiadis.xls -> anastassiadis_matrix_released.csv (178 x 300) -> anastassiadis_matrix.csv (curated, 291)
- Metz 2011: metz.xls -> metz_matrix.csv (704 x 172, pKi, >= 50-kinase subset)
- Not used: Gao 2013 (same modality as Anastassiadis, not redistributable),
  Patricelli 2011 (too few inhibitors for a ranking analysis)

## Clinical outcomes analysis (Limitations)
- Proportions of FAERS reports that were serious, fatal or involved
  hospitalization (n = 46) and FDA-label discontinuation rates (n = 33) vs median
  selectivity rank: all |r| <= 0.20, all p > 0.17 (clinical_outcomes_correlation.py)
- Only nominal signal: raw FAERS report count vs ratio rank, r = 0.32, p = 0.03,
  likely prescription volume
- Grade 3/4 rates: too few drugs with data for analysis
- Confounded by indication and prescription volume; reported as a null result
