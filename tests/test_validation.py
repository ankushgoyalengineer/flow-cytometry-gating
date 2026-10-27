import pandas as pd
from flowtool.validation import score_crosstab, validate_against_labels
import pytest
from pathlib import Path
from flowtool.fcs_io import inspect_fcs_file
from flowtool.gating import apply_gating, CUTOFFS

DATA = Path("data/FlowCAP-I HSCT")
LABELS = Path("data/hsct_ks_labels.csv")
needs_data = pytest.mark.skipif(not DATA.exists(), reason="dataset not downloaded")
needs_labels = pytest.mark.skipif(not LABELS.exists(), reason="dataset not available")

SUBJECT_TO_FILE = {5: "005.fcs", 6: "006.fcs", 9: "009.fcs", 12: "012.fcs"}

@pytest.fixture
def temp_df():
    # This function creates and returns a fresh DataFrame every time it is called
    return pd.DataFrame({0: [0, 2], 1: [8, 2], 2: [0, 10]}, 
    index=["my group A", "my group B"])

def test_confirm_scores_are_correct(temp_df):
    scores = score_crosstab(temp_df).set_index("expert_class")

    assert 0 not in scores.index
    assert len(scores) == 2
    assert scores.loc[1, "my_group"] == "my group A"
    assert scores.loc[1, "precision"] == pytest.approx(1.0, abs=0.001)
    assert scores.loc[1, "recall"] == pytest.approx(0.8, abs=0.001)
    assert scores.loc[1, "f_measure"] == pytest.approx(0.889, abs=0.001)

    assert scores.loc[2, "my_group"] == "my group B"
    assert scores.loc[2, "precision"] == pytest.approx(0.714, abs=0.001)
    assert scores.loc[2, "recall"] == pytest.approx(1.0, abs=0.001)
    assert scores.loc[2, "f_measure"] == pytest.approx(0.833, abs=0.001)

def test_confirm_small_classes_flagged(temp_df):
    temp_df[3]=[3,0]
    scores = score_crosstab(temp_df)
    target_row = scores[scores["expert_class"] == 3]
    assert target_row["note"].iloc[0] == "too few cells to score"

@needs_data
def test_file001_counts():
    
    _, file001_df = inspect_fcs_file(DATA/"001.fcs")
    summary_counts = apply_gating(file001_df, **CUTOFFS)

    assert summary_counts["live"]==9854
    assert summary_counts["recipient"]==9461
    assert summary_counts["donor"]==4
    assert summary_counts["cd45_neg"]==389

@needs_data
@needs_labels
def test_validation_reproduces_result():
    validated_df = validate_against_labels(LABELS, DATA, SUBJECT_TO_FILE, CUTOFFS )
    validated_df = validated_df[validated_df["note"]==""]

    assert min(validated_df["f_measure"]) == pytest.approx(0.957, abs=0.001)
    assert (validated_df["f_measure"].median()) == pytest.approx(0.985, abs=0.001)



