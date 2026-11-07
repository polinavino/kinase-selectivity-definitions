# Kinase Inhibitor Selectivity: Research Notes

## Project status
Second revision for Molecular Informatics (manuscript 7756455), responding to
Reviewer 1's round-2 review (`paper/reviewer1_round2_review.txt`; response in
`paper/response_to_reviewers_round2.txt`). Analysis complete on four datasets;
every number in the paper is printed by a script into `outputs/` (see README).
Preprint: ChemRxiv (https://doi.org/10.26434/chemrxiv.15001618/v1).

## Paper structure (after consolidation)
Introduction (merged with Background); Methods (datasets and target selection,
score implementations, ranking comparisons, criteria and the gated entropy,
recovery design); Results (analysis choices and agreement; parameter
sensitivity and profile associations; criteria and the gated entropy;
recovery from small panels); Discussion; Conclusion. Main text about 9,000
words; SI Sections S1-S8 follow the same
order.

## Research question
Existing kinase inhibitor selectivity definitions (S-score, entropy, Gini,
ratio) are used interchangeably. The paper:
1. Characterizes where and why they disagree, on four datasets curated to one
   reference target per human protein kinase (curate_panels.py)
2. Relates rank instability to features of the binding profile
3. States four criteria (D1 activity eligibility, D2 distributional
   consistency, D3 monotonicity under weak-target addition, D4 orientation on
   profiles with equally strong targets) and evaluates
   the four definitions, three comparators and a candidate against them
4. Evaluates recovery from random subpanels separately (panel_recovery.py)

## Key findings
- Agreement depends on analysis choices: which compounds are included, how
  the entropy weights affinities, and how floor-coded pairs are treated
  (Results, first section; Table 3). With entropy and Gini on a common baseline-subtracted
  log-affinity scale, among compounds with an active kinase, S-score/entropy/
  Gini agree (r = 0.81-0.999) and the ratio agrees less with each. With all
  compounds, on Klaeger the Gini sits closer to the ratio (0.66 vs 0.51; not
  significant); on Anastassiadis the split comes mainly from the 45 zero-active
  compounds; on Metz rankings largely follow coverage.
- Published (association-constant) entropy correlates most with the ratio on
  Davis, Klaeger and Metz under both treatments of floor-coded pairs. Against
  the log-affinity entropy: r = 0.37/-0.31 (Davis/Klaeger) with floor-coded
  pairs at the floor, 0.74/0.79 with them given zero weight (Figure 1).
- Compounds with no active kinase (Klaeger, 22 of 222): all 20 most unstable
  when all compounds are ranked together; aggregated rank sigma 77.9 vs 32.2,
  from the definitions' conventions for such profiles. Among the 20 most
  unstable active compounds (ranked among themselves), 5 have a top1-top2 gap
  below 0.1, 5 have >= 10 active kinases and 10 neither.
- Ranking population: profile analyses of active compounds rank them among
  themselves (D1). Ranked among all 222 Klaeger compounds, skewness was the
  Gini's strongest correlate (+0.26); ranked among active compounds it is
  +0.07 (p = 0.31), and n_active correlates with Gini instability on both pKd
  datasets.
- Relations: ratio-gap on Davis, Klaeger, Metz (-0.44, -0.28, -0.26), not on
  Anastassiadis (0.00; -0.15 without saturated compounds); entropy-n_active on
  all four (r = -0.30 to -0.42; Davis p = 0.01, not significant after
  Bonferroni).
- Properties: no boundary parameter tested satisfies D2 (rule: within eps on every
  dataset and, for the ratio, every k). S-score threshold reverses 8.4-13.0%
  of pairs on three datasets; entropy/Gini 7.9-12.4%; candidate gate
  9.2-14.4%; ratio clipping level b in R_b(k) = x(1) - max(x(k+1), b) up to
  13.7%. D3 holds for the active count, the S-score at tau = x0, the entropy
  with beta = tau*, the ratio, P_max and the gated entropy (not the
  published entropy); not for the S-score with
  tau > x0, the Gini or the entropy with beta < tau*. The former D2 (bounded
  gap sensitivity) was dropped. D4 holds for every score studied when
  floor-coded pairs have zero weight; the entropy and Gini fail it at
  baselines below the floor (part of the Anastassiadis sweep); the reversed
  active count (which meets D1-D3 as the gated entropy does), the gated
  entropy without its gate and the published entropy with floor-coded pairs
  kept at the floor fail it (P_max keeps D4 under either floor treatment). Proposition 3: turning a floor-coded pair
  into a target no stronger than the weakest weighted one never raises -H,
  P_max or the Gini. Proposition 2: for fixed weights, Shannon's conditions
  make the gated entropy unique.
- Candidate (gated smooth-hinge entropy): D1 and D3 (proved) while weighting
  binding below tau* continuously; ranks almost identically to the entropy
  (0.994-0.998) and to the S-score at x0 (0.998-0.9997). Noise check with the
  D1 gate reapplied: Spearman >= 0.978.
- Recovery (D1 applied, Figure 3): measures that follow the share of targets with
  detected binding (S-score at x0, entropy beta = x0, Gini, candidate) reach
  Spearman 0.9 at 20/50/30/65 targets (Davis/Klaeger/Anastassiadis/Metz);
  S-score at tau*/count 30/80/50/50; published entropy 300/200/-/150; partition
  index P_max 300/200/-/>150; ratio more than the largest subpanel tested
  except on Klaeger (200). At p*, 62-90% of subpanels exceed 0.90. Coverage
  limits pair recovery.
- Released panels (SI S1, Table S1): Davis 26 of 68 compounds have same-gene top two
  constructs; Klaeger curated to the 217 targets the source annotates as
  protein kinases (BCR, TRRAP removed: not annotated, possibly co-enriched).

## Datasets
- Davis 2011: davis_affinity_released.csv (68 x 433) -> davis_affinity.csv (curated, 363)
- Klaeger 2017: aan4368_Table_S2.xlsx -> klaeger_matrix_released.csv (222 x 343) -> klaeger_matrix.csv (curated, 217; values below the 5.0 floor set to it)
- Anastassiadis 2011: anastassiadis.xls -> anastassiadis_matrix_released.csv (178 x 300) -> anastassiadis_matrix.csv (curated, 291)
- Metz 2011: metz.xls -> metz_matrix.csv (704 x 172, pKi, >= 50-kinase subset)
- Not used: Gao 2013 (same modality as Anastassiadis, not redistributable),
  Patricelli 2011 (too few inhibitors for a ranking analysis)

## Clinical outcomes analysis (Limitations)
- Proportions of FAERS reports that were serious, fatal or involved
  hospitalization (n = 46) vs median selectivity rank: all |r| <= 0.20, all
  p > 0.19 (clinical_outcomes_correlation.py)
- FDA-label discontinuation rates were dropped: an audit against the stored
  label text found dose-reduction rates recorded as discontinuation
  (ceritinib 45%, pazopanib 36%), values not in the stored text (ibrutinib,
  lapatinib, sorafenib) and a combination-therapy rate (axitinib)
- Only nominal signal: raw FAERS report count vs ratio rank, r = 0.30, p = 0.04,
  likely prescription volume
- Confounded by indication and prescription volume; reported as a null result
