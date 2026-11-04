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
1. Characterizes where and why they disagree, on four datasets
2. Identifies three sources of rank instability
3. Motivates four properties (D1-D4) for a well-formed measure
4. Shows with a candidate measure that D1-D4 can be met together

## Key findings
- Two families, with entropy and Gini on a common baseline-subtracted
  log-affinity scale: S-score/entropy/Gini agree (within-family r = 0.72-1.00),
  the ratio agrees with none (cross-family r = 0.14-0.62). On Anastassiadis the
  split comes mainly from the 44 zero-active compounds (without them the ratio
  correlates 0.79-0.84 with the family). Under the published
  association-constant entropy the grouping does not hold (SI Table S3).
- Instability sources (Klaeger): Type 1, zero-active compounds (16 of 222;
  pooled rank sigma 74.0 vs 31.2, driven by the definitions disagreeing about
  them: S-score and entropy rank them near the top, the ratio in the
  bottom third); Type 2, near-tied top targets (gap vs ratio instability
  r = -0.341); Type 3, broad profiles, often with a dominant target.
- Replicated relations: ratio-gap on all four datasets (Anastassiadis only
  after excluding 8 ceiling-saturated compounds); entropy-n_active on all four
  (r = -0.28 to -0.42).
- Properties: no existing definition satisfies all four; D1 is the one all four
  violate. Candidate (gated, smooth-hinge effective-target number) satisfies
  all four, D3 only with its activity boundary fixed at the assay floor; D4 is
  proved for T >= 0.58 (tau* - x0).
- Candidate and entropy rank these datasets almost identically
  (Spearman 0.994-0.999); panel-size p* = 50 (Davis, Anastassiadis), 80
  (Klaeger, Metz) for both.
- Operationalization check: literal vs log-affinity entropy r = 0.42 (Davis),
  -0.34 (Klaeger) (Figure 1).

## Datasets
- Davis 2011: davis_affinity.csv, davis_drugs.csv, davis_proteins.csv (68 x 433)
- Klaeger 2017: aan4368_Table_S2.xlsx -> klaeger_matrix.csv (222 x 343)
- Anastassiadis 2011: anastassiadis.xls -> anastassiadis_matrix.csv (178 x 300, % inhibition)
- Metz 2011: metz.xls -> metz_matrix.csv (704 x 172, pKi, >= 50-kinase subset)
- Not used: Gao 2013 (same modality as Anastassiadis, not redistributable),
  Patricelli 2011 (too few inhibitors for a ranking analysis)

## Clinical outcomes analysis (Limitations)
- FAERS rates (n = 46) and FDA-label discontinuation rates (n = 33) vs median
  selectivity rank: all |r| < 0.17, all p > 0.27 (clinical_outcomes_correlation.py)
- Only nominal signal: raw FAERS report count vs ratio rank, r = 0.31, p = 0.035,
  likely prescription volume
- Grade 3/4 rates: too few drugs with data for analysis
- Confounded by indication and prescription volume; reported as a null result
