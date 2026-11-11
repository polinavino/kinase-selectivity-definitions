"""
Reproduces the FAERS correlation analysis referenced in the Limitations
(SI Section S7): does in vitro selectivity (median rank under each of the four
definition families) predict clinical adverse-event outcomes for FDA-approved
Klaeger compounds?

Uses the frozen FAERS snapshot selectivity_outcomes_merged.csv: per-drug FAERS
report counts and rates (serious_rate, death_rate, hosp_rate) for n=46
FDA-approved Klaeger compounds with FAERS records. FAERS report counts grow
over time, so faers_pull.py's live API pull does not reproduce them exactly;
this file is the snapshot the paper's numbers are computed from.

FDA-label discontinuation rates (clinical_safety_data.csv) are no longer used:
an audit against the stored label text found dose-reduction rates recorded as
discontinuation rates (ceritinib, pazopanib), values not found in the stored
text (ibrutinib, lapatinib, sorafenib) and a combination-therapy rate
(axitinib).
"""
import script_logging; script_logging.capture(__file__)
import os
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))

RANK_COLS = ['median_rank_s_score', 'median_rank_entropy', 'median_rank_gini', 'median_rank_ratio']
RATE_COLS = ['serious_rate', 'death_rate', 'hosp_rate']


def analyze():
    so = pd.read_csv(f"{HERE}/selectivity_outcomes_merged.csv")
    # The FAERS counts in the merged file are a frozen snapshot, but its rank
    # columns are not: take each drug's median ranks from the current Klaeger
    # analysis (klaeger_analysis.py) so the correlations track the paper.
    ranks = pd.read_csv(f"{HERE}/klaeger_selectivity_results.csv")[['drug', 'rank_std'] + RANK_COLS]
    so = so.drop(columns=RANK_COLS + ['rank_std', 'n_active'], errors='ignore').merge(ranks, on='drug', how='left')
    panel = pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0)
    so['n_active'] = so['drug'].map((panel > 6.0).sum(axis=1))
    assert so[RANK_COLS].notna().all().all(), "a FAERS drug is missing from klaeger_selectivity_results.csv"
    merged = so

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
    # Refresh only the derived selectivity fields; retain the frozen report
    # totals and rates rather than querying a different live FAERS snapshot.
    merged.to_csv(f"{HERE}/selectivity_outcomes_merged.csv", index=False)
    df.to_csv(f"{HERE}/clinical_outcomes_correlation.csv", index=False)
    print(df.to_string(index=False))
    print(f"\nmax |r| across all rate x rank pairs = {df['r'].abs().max():.3f}, "
          f"min p = {df['p'].min():.4f}")
    print(f"n (FAERS rates) = {merged[RATE_COLS[0]].notna().sum()}")
    print(f"\ntotal_reports vs median_rank_ratio: r = {r_tot:.4f}, p = {p_tot:.4f}")
