"""Regression checks for raw regeneration and stale Anastassiadis caches."""
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from panel_curation import load_anastassiadis_panel

ROOT = Path(__file__).resolve().parents[1]


class AnastassiadisPanelTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.expected = pd.read_csv(ROOT / "anastassiadis_matrix.csv", index_col=0)

    def tearDown(self):
        self.temp.cleanup()

    def check_panel(self, actual):
        self.assertEqual(actual.shape, (178, 291))
        self.assertEqual(list(actual.columns), list(self.expected.columns))
        self.assertEqual(list(actual.index), list(self.expected.index))
        np.testing.assert_allclose(actual.values, self.expected.values, rtol=0, atol=1e-12)

    def test_regeneration_from_released_panel_keeps_curation(self):
        shutil.copy2(ROOT / "anastassiadis_matrix_released.csv", self.directory)
        self.check_panel(load_anastassiadis_panel(self.directory))
        self.check_panel(pd.read_csv(self.directory / "anastassiadis_matrix.csv", index_col=0))

    def test_regeneration_from_raw_spreadsheet_keeps_curation(self):
        shutil.copy2(ROOT / "anastassiadis.xls", self.directory)
        self.check_panel(load_anastassiadis_panel(self.directory))
        released = pd.read_csv(self.directory / "anastassiadis_matrix_released.csv", index_col=0)
        self.assertEqual(released.shape, (178, 300))

    def test_released_panel_under_curated_filename_is_rejected(self):
        shutil.copy2(ROOT / "anastassiadis_matrix_released.csv", self.directory)
        shutil.copy2(ROOT / "anastassiadis_matrix_released.csv",
                     self.directory / "anastassiadis_matrix.csv")
        with self.assertRaisesRegex(ValueError, "target identities"):
            load_anastassiadis_panel(self.directory)

    def test_same_size_cache_with_wrong_target_is_rejected(self):
        shutil.copy2(ROOT / "anastassiadis_matrix_released.csv", self.directory)
        wrong = self.expected.rename(columns={self.expected.columns[0]: "wrong-target"})
        wrong.to_csv(self.directory / "anastassiadis_matrix.csv")
        with self.assertRaisesRegex(ValueError, "target identities"):
            load_anastassiadis_panel(self.directory)


if __name__ == "__main__":
    unittest.main()
