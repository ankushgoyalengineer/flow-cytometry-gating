import pandas as pd

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