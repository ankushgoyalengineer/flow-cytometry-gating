import pytest
import numpy as np
import pandas as pd
from flowtool.gating import pct, apply_gating, assign_labels, CUTOFFS


def test_pct_empty_parent_is_nan():
    assert np.isnan(pct([], []))


def test_pct_equals_right_percentage():
    """pct uses the lengths of its arguments, so [2, 2] within [4, 4, 4, 4] is 2 of 4 = 50%."""
    assert pct([2, 2], [4, 4, 4, 4]) == 50


@pytest.fixture
def temp_df():
    # Values chosen so each row falls in exactly one group with the default CUTOFFS (700, 550, 450, 500)
    return pd.DataFrame({
        "FL1-H": [800, 100, 100, 800, 800, 300],
        "FL2-H": [600, 800, 300, 800, 300, 100],
        "FL3-H": [900, 100, 100, 100, 100, 100],
        "FL4-H": [500, 800, 800, 100, 100, 300]
    })


def test_apply_gating_result(temp_df):
    test_summary = apply_gating(temp_df, **CUTOFFS)

    assert test_summary["live"] == 5
    assert test_summary["events"] - test_summary["live"] == 1
    assert test_summary["cd45_neg"] == 1
    assert test_summary["donor"] == 2
    assert test_summary["recipient"] == 2
    assert test_summary["mac1_pos_recipient"] == 1
    assert test_summary["mac1_neg_recipient"] == 1
    assert test_summary["mac1_pos_donor"] == 1
    assert test_summary["mac1_neg_donor"] == 1


def test_values_at_cutoff_go_to_correct_side(temp_df):
    # Move one value in each of four rows to exactly its cut-off
    temp_df.loc[5, "FL1-H"] = CUTOFFS["cd45_1_cutoff"]
    temp_df.loc[4, "FL2-H"] = CUTOFFS["mac1_cutoff"]
    temp_df.loc[0, "FL3-H"] = CUTOFFS["pi_cutoff"]
    temp_df.loc[2, "FL4-H"] = CUTOFFS["cd45_2_cutoff"]

    test_summary_2 = apply_gating(temp_df, **CUTOFFS)

    # A value exactly at a cut-off counts as above it (>=)
    assert test_summary_2["live"] == 5
    assert test_summary_2["recipient"] == 2
    assert test_summary_2["donor"] == 3
    assert test_summary_2["mac1_pos_donor"] == 2


def test_apply_gating_assign_labels_agree(temp_df):
    test_summary = apply_gating(temp_df, **CUTOFFS)
    labelled = assign_labels(temp_df, **CUTOFFS)
    counts = pd.Series(labelled).value_counts()

    assert counts["dead"] == test_summary["events"] - test_summary["live"]
    assert counts["recipient Mac1+"] == test_summary["mac1_pos_recipient"]
    assert counts["recipient Mac1-"] == test_summary["mac1_neg_recipient"]
    assert counts["donor Mac1+"] == test_summary["mac1_pos_donor"]
    assert counts["donor Mac1-"] == test_summary["mac1_neg_donor"]
    assert counts["CD45-negative"] == test_summary["cd45_neg"]
    assert "unassigned" not in counts