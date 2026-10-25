import numpy as np
import pandas as pd

from flowtool.fcs_io import inspect_fcs_file
from flowtool.gating import assign_labels

def score_crosstab(tab, min_events=10):
    """For each expert class (except the catch-all class 0), find the group of mine
    that overlaps it most, then calculate precision, recall and F-measure."""
    rows = []
    for expert_class in tab.columns:
        if expert_class == 0:                      # class 0 = experts' "everything else"
            continue
        column = tab[expert_class]                 # how my groups are spread across this expert class
        expert_total = column.sum()                # all cells the experts put in this class

        best_group = column.idxmax()               # my group with the most cells in this class
        both = column[best_group]                  # cells we both put in this population
        my_total = tab.loc[best_group].sum()       # all cells I put in my matched group

        precision = both / my_total                # of my cells, share the experts agree with
        recall = both / expert_total               # of the experts' cells, share I found
        f_measure = 2 * precision * recall / (precision + recall)

        rows.append({
            "expert_class": expert_class,
            "my_group": best_group,
            "expert_cells": int(expert_total),
            "precision": precision,
            "recall": recall,
            "f_measure": f_measure,
            "note": "too few cells to score" if expert_total < min_events else "",
        })
    return pd.DataFrame(rows)

KS_TO_FL = {"FITC.CD45.1": "FL1-H", "PE.Ly65Mac1": "FL2-H",
            "PI.LiveDead": "FL3-H", "APC.CD45.2": "FL4-H"}

def validate_against_labels(labels_csv, folder, subject_to_file, cutoffs):
    """Attach the ks expert labels to each matching FCS file, label every event with
    my gates, and score the agreement file by file."""
    fl_cols = list(KS_TO_FL.values())
    author_df = pd.read_csv(labels_csv).rename(columns=KS_TO_FL)

    all_scores = []
    for subject, file_name in subject_to_file.items():
        _, df = inspect_fcs_file(f"{folder}/{file_name}")
        subject_rows = author_df[author_df["subject"] == subject].reset_index(drop=True)

        assert len(subject_rows) == len(df), f"{file_name}: different number of events"
        assert np.allclose(subject_rows[fl_cols].to_numpy(), df[fl_cols].to_numpy()), \
            f"{file_name}: values don't match"

        df["expert_label"] = subject_rows["label"].to_numpy()
        df["my_label"] = assign_labels(df, **cutoffs)

        tab = pd.crosstab(df["my_label"], df["expert_label"])
        scores = score_crosstab(tab)
        scores.insert(0, "file", file_name)
        all_scores.append(scores)

    return pd.concat(all_scores, ignore_index=True)