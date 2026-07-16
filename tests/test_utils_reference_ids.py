import unittest

import pandas as pd

from fl_ai_reagent_stocker.utils import (
    format_brace_entry,
    format_pmcid_display,
    format_pmid_display,
    format_pmid_list_display,
    format_reference_id_columns,
    format_reference_id_display,
    parse_stock_candidate_label,
)


class TestReferenceIdFormatters(unittest.TestCase):
    def test_format_pmid_display_adds_prefix(self):
        self.assertEqual(format_pmid_display("12345678"), "PMID12345678")
        self.assertEqual(format_pmid_display("PMID12345678"), "PMID12345678")
        self.assertEqual(format_pmid_display(""), "-")
        self.assertEqual(format_pmid_display("-"), "-")

    def test_format_pmcid_display_adds_prefix(self):
        self.assertEqual(format_pmcid_display("123456"), "PMCID123456")
        self.assertEqual(format_pmcid_display("PMC123456"), "PMCID123456")
        self.assertEqual(format_pmcid_display("PMCID123456"), "PMCID123456")

    def test_format_reference_id_display_prefers_pmid(self):
        self.assertEqual(
            format_reference_id_display("123", "456"),
            "PMID123",
        )
        self.assertEqual(
            format_reference_id_display("", "456"),
            "PMCID456",
        )

    def test_format_pmid_list_display(self):
        self.assertEqual(
            format_pmid_list_display("123; 456; 123"),
            "PMID123; PMID456",
        )
        self.assertEqual(format_pmid_list_display("-"), "-")

    def test_parse_stock_candidate_label(self):
        self.assertEqual(
            parse_stock_candidate_label("(7126, BDSC)"),
            ("7126", "Bloomington"),
        )
        self.assertEqual(parse_stock_candidate_label("(8123)"), ("8123", ""))

    def test_format_brace_entry(self):
        self.assertEqual(
            format_brace_entry("tim-GAL4", "BDSC", "PMID123", "sleep defect"),
            "{tim-GAL4, BDSC, PMID123, sleep defect}",
        )

    def test_format_reference_id_columns(self):
        df = pd.DataFrame(
            {
                "PMID": ["123; 456", "-"],
                "PMCID": ["PMC789", ""],
                "FBal PMID": ["999", "-"],
            }
        )
        formatted = format_reference_id_columns(df)
        self.assertEqual(formatted["PMID"].iloc[0], "PMID123; PMID456")
        self.assertEqual(formatted["PMCID"].iloc[0], "PMCID789")
        self.assertEqual(formatted["FBal PMID"].iloc[0], "PMID999")


if __name__ == "__main__":
    unittest.main()
