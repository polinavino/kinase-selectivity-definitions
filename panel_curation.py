"""Shared, deterministic target selection for the Anastassiadis panel."""
from pathlib import Path

import numpy as np
import pandas as pd


def anastassiadis_gene(name):
    gene = name.split("/")[0].split(" cyclin")[0].strip()
    return {"LYN B": "LYN", "PKCb1": "PKCb", "PKCb2": "PKCb",
            "PKG1a": "PKG1", "PKG1b": "PKG1"}.get(gene, gene)


def curate_anastassiadis(frame):
    """Retain the first-listed complex or isoform of each kinase gene."""
    first = {}
    for column in frame.columns:
        first.setdefault(anastassiadis_gene(column), column)
    columns = [column for column in frame.columns
               if first[anastassiadis_gene(column)] == column]
    return frame.loc[:, columns].copy()


def load_anastassiadis_panel(directory):
    """Load a curated panel, regenerating missing files without bypassing curation.

    Validate target identities and compound order even for an existing cache.
    A released panel accidentally stored under the curated name fails loudly.
    """
    directory = Path(directory)
    released_path = directory / "anastassiadis_matrix_released.csv"
    curated_path = directory / "anastassiadis_matrix.csv"
    if released_path.exists():
        released = pd.read_csv(released_path, index_col=0)
    else:
        raw = pd.read_excel(directory / "anastassiadis.xls", sheet_name="Sheet1", header=None)
        residual = raw.iloc[3:, 1:].apply(pd.to_numeric, errors="coerce").values
        values = np.nan_to_num(np.clip(100 - residual, 0, 100), nan=0.0).T
        released = pd.DataFrame(values, index=raw.iloc[1, 1:].tolist(),
                                columns=raw.iloc[3:, 0].tolist())
        released.to_csv(released_path)
    expected = curate_anastassiadis(released)
    if not curated_path.exists():
        expected.to_csv(curated_path)
        return expected
    cached = pd.read_csv(curated_path, index_col=0)
    if list(cached.columns) != list(expected.columns) or list(cached.index) != list(expected.index):
        raise ValueError("Anastassiadis cache does not match the curated target identities "
                         "and compound order; regenerate it with curate_panels.py.")
    if not np.allclose(cached.values, expected.values, rtol=0, atol=1e-12, equal_nan=True):
        raise ValueError("Anastassiadis cache differs from the released data; "
                         "regenerate it with curate_panels.py.")
    return cached
