"""Read-only checks supporting the November 11, 2026 manuscript audit.

Run from any directory. Prints JSON; does not execute analysis modules or write
their outputs. Function definitions are loaded through AST to avoid their
top-level logging, data writes and figure generation.
"""
import ast
import json
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
FAM = ["s_score", "entropy", "gini", "ratio"]


def functions(filename):
    tree = ast.parse((ROOT / filename).read_text())
    scope = dict(np=np, pd=pd, rankdata=rankdata, spearmanr=spearmanr, FAM=FAM)
    module = ast.Module(
        body=[node for node in tree.body if isinstance(node, ast.FunctionDef)],
        type_ignores=[],
    )
    exec(compile(module, filename, "exec"), scope)
    return scope


def csv(filename, **kwargs):
    return pd.read_csv(ROOT / filename, **kwargs)


def main():
    out = {}
    d = csv("davis_affinity.csv").pivot(
        index="Drug_Index", columns="Protein_Index", values="Affinity")
    k = csv("klaeger_matrix.csv", index_col=0)
    a = csv("anastassiadis_matrix.csv", index_col=0)
    m = csv("metz_matrix.csv", index_col=0)
    raw_k = pd.read_excel(ROOT / "aan4368_Table_S2.xlsx", sheet_name="Kinobeads")
    raw_k["pKd"] = 9 - np.log10(raw_k["Apparent Kd"])
    high = raw_k[raw_k["Target Classification"] == "High confidence"]
    rebuilt_k = high.groupby(["Drug", "Gene Name"])["pKd"].mean().unstack(fill_value=5)
    rebuilt_k = rebuilt_k.reindex(index=k.index, columns=k.columns).clip(lower=5)
    raw_a = pd.read_excel(ROOT / "anastassiadis.xls", sheet_name="Sheet1", header=None)
    inhibition = np.nan_to_num(np.clip(
        100 - raw_a.iloc[3:, 1:].apply(pd.to_numeric, errors="coerce").values,
        0, 100), nan=0).T
    rebuilt_a = pd.DataFrame(inhibition, index=raw_a.iloc[1, 1:].tolist(),
                             columns=raw_a.iloc[3:, 0].tolist())
    rebuilt_a = rebuilt_a.loc[a.index, a.columns]
    raw_m = pd.read_excel(ROOT / "metz.xls", sheet_name=0).set_index("Cmpd_ID")
    rebuilt_m = raw_m.loc[m.index, m.columns].apply(pd.to_numeric, errors="coerce").fillna(4)
    out["raw_matrix_max_absolute_errors"] = {
        "Klaeger": float(np.max(np.abs(rebuilt_k.values - k.values))),
        "Anastassiadis": float(np.max(np.abs(rebuilt_a.values - a.values))),
        "Metz": float(np.max(np.abs(rebuilt_m.values - m.values))),
    }

    docx = ROOT / "paper/molecular_informatics_submission.docx"
    with ZipFile(docx) as archive:
        xml = ET.fromstring(archive.read("word/document.xml"))
    text = " ".join(node.text or "" for node in xml.iter(
        "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t"))
    out["upload_word_document"] = {
        "path": str(docx.relative_to(ROOT)),
        "contains_old_title": "Definitional Instability in Kinase Inhibitor Selectivity" in text,
        "contains_softplus": "softplus" in text,
        "contains_old_D2": "bounded gap sensitivity" in text,
        "contains_new_title": "Analysis Choices and Ranking Stability" in text,
    }

    # The published cross-dataset diagnostic uses an aggregated sigma for
    # every definition. Repeat with each definition's own parameter sigma.
    defs = functions("claim_checks.py")
    sweep = lambda matrix: defs["sweeps"](
        matrix, np.arange(5.5, 8.25, .25), np.arange(5, 6.75, .25), 5)
    rd, rk = sweep(d.values), sweep(k.values)
    overlap = csv("davis_klaeger_overlap.csv")
    di = d.index.get_indexer(overlap.Davis_Drug_Index)
    ki = k.index.get_indexer(overlap.Klaeger_Name)
    agg_d = np.vstack(list(rd.values())).std(0)
    agg_k = np.vstack(list(rk.values())).std(0)
    associations = {}
    for family in FAM:
        displacement = np.abs(rankdata(np.median(rd[family], 0)[di]) -
                              rankdata(np.median(rk[family], 0)[ki]))
        row = {}
        for label, instability in [
            ("Davis_aggregate", agg_d[di]), ("Klaeger_aggregate", agg_k[ki]),
            ("Davis_within_definition", rd[family].std(0)[di]),
            ("Klaeger_within_definition", rk[family].std(0)[ki]),
        ]:
            result = spearmanr(displacement, instability)
            row[label] = dict(r=float(result.statistic), p=float(result.pvalue))
        associations[family] = row
    out["instability_vs_cross_dataset_displacement"] = associations

    # Filter sensitivity on the intersection versus each filter's eligible set.
    relaxed = raw_k.dropna(subset=["pKd"]).groupby(
        ["Drug", "Gene Name"])["pKd"].mean().unstack(fill_value=5)
    relaxed = relaxed.reindex(index=k.index, columns=k.columns).clip(lower=5)
    hk, hr = (k.values > 6).any(1), (relaxed.values > 6).any(1)
    cells = defs["corr_cells"]
    ch, cr = cells(sweep(k.values[hk])), cells(sweep(relaxed.values[hr]))
    ci, cj = cells(sweep(k.values[hk & hr])), cells(sweep(relaxed.values[hk & hr]))
    out["filter_sensitivity"] = dict(
        high_confidence_eligible=int(hk.sum()), relaxed_eligible=int(hr.sum()),
        newly_eligible=k.index[hr & ~hk].tolist(),
        max_change_intersection=float(max(abs(ci[key] - cj[key]) for key in ci)),
        max_change_own_eligible_sets=float(max(abs(ch[key] - cr[key]) for key in ch)),
    )

    # Illustrate the uncertainty hidden by the clinical null wording.
    rates = csv("selectivity_outcomes_merged.csv")[["drug", "hosp_rate"]]
    rates = rates.merge(csv("klaeger_selectivity_results.csv"), on="drug")
    x, y = rates.hosp_rate.values, rates.median_rank_gini.values
    rng = np.random.default_rng(42)
    bootstrap = []
    for _ in range(3000):
        idx = rng.integers(0, len(x), len(x))
        bootstrap.append(spearmanr(x[idx], y[idx]).statistic)
    result = spearmanr(x, y)
    out["clinical_hospitalization_vs_gini"] = dict(
        n=len(x), r=float(result.statistic), p=float(result.pvalue),
        percentile_95_CI=np.quantile(bootstrap, [.025, .975]).tolist(),
        bootstrap_draws=3000, seed=42,
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
