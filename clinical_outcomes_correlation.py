"""
Reproduces the FAERS / FDA-label discontinuation-rate correlation analysis
referenced in the Discussion and Limitations sections: does in vitro
selectivity (median rank under each of the four definition families) predict
clinical adverse-event outcomes for FDA-approved Klaeger compounds?

Merges the two frozen clinical data files committed to this repository:
  * selectivity_outcomes_merged.csv -- per-drug FAERS report counts and rates
    (serious_rate, death_rate, hosp_rate) plus each drug's median selectivity
    rank under every definition family, for n=46 FDA-approved Klaeger
    compounds with FAERS records. FAERS report counts grow over time (this is
    not what faers_pull.py's live API pull reproduces); this file is the
    frozen snapshot the paper's numbers are computed from.
  * clinical_safety_data.csv -- discontinuation_rate and grade34_rate manually
    extracted from FDA drug label text (label_texts.json, extended_labels.json)
    for n=33/35 of those drugs with a reported discontinuation rate.

No script produces these two files themselves (they were built/merged by
hand from FAERS pulls and manual label review); this script is the
permanent, reproducible record of the correlation analysis computed from them.
"""
import script_logging; script_logging.capture(__file__)
import os
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))

RANK_COLS = ['median_rank_s_score', 'median_rank_entropy', 'median_rank_gini', 'median_rank_ratio']
RATE_COLS = ['serious_rate', 'death_rate', 'hosp_rate', 'discontinuation_rate']


def analyze():
    so = pd.read_csv(f"{HERE}/selectivity_outcomes_merged.csv")
    cs = pd.read_csv(f"{HERE}/clinical_safety_data.csv")
    merged = so.merge(cs[['drug', 'discontinuation_rate', 'grade34_rate']], on='drug', how='left')

    rows = []
    for rate in RATE_COLS:
        for rank in RANK_COLS:
            sub = merged[[rate, rank]].dropna()
            r, p = spearmanr(sub[rate], sub[rank])
            rows.append(dict(rate=rate, rank=rank, n=len(sub), r=round(r, 3), p=round(p, 4)))
    df = pd.DataFrame(rows)

    r_tot, p_tot = spearmanr(merged['total_reports'], merged['median_rank_ratio'])
    return df, merged, r_tot, p_tot


if __name__ == "__main__":
    df, merged, r_tot, p_tot = analyze()
    df.to_csv(f"{HERE}/clinical_outcomes_correlation.csv", index=False)
    print(df.to_string(index=False))
    print(f"\nmax |r| across all rate x rank pairs = {df['r'].abs().max():.3f}, "
          f"min p = {df['p'].min():.4f}")
    print(f"n (FAERS rates) = {merged[RATE_COLS[0]].notna().sum()}, "
          f"n (discontinuation_rate) = {merged['discontinuation_rate'].notna().sum()}")
    print(f"\ntotal_reports vs median_rank_ratio: r = {r_tot:.3f}, p = {p_tot:.4f}")
