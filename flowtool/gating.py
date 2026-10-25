import numpy as np
import pathlib
import pandas as pd
from flowtool.fcs_io import inspect_fcs_file


CUTOFFS = {"pi_cutoff": 700, "cd45_2_cutoff": 550, "cd45_1_cutoff": 450, "mac1_cutoff": 500}

def pct(part, parent):
    """Percentage of part within parent; NaN if the parent is empty."""
    if len(parent) == 0:
        return np.nan
    return 100 * len(part) / len(parent)


def apply_gating(df_sample, pi_cutoff, cd45_2_cutoff, cd45_1_cutoff, mac1_cutoff):
    """Gate one sample and return counts and percentages for every population."""
    # Live / dead
    live = df_sample[df_sample["FL3-H"] < pi_cutoff]
    live_percent = pct(live, df_sample)

    # Recipient: CD45.2-high
    recipient = live[live["FL4-H"] >= cd45_2_cutoff]
    recipient_percent = pct(recipient, live)

    # CD45.2-low, split by CD45.1 into donor and CD45-negative
    cd45_low = live[live["FL4-H"] < cd45_2_cutoff]
    donor = cd45_low[cd45_low["FL1-H"] >= cd45_1_cutoff]
    donor_percent = pct(donor, live)
    cd45_neg = cd45_low[cd45_low["FL1-H"] < cd45_1_cutoff]
    cd45_neg_percent = pct(cd45_neg, live)
    assert len(donor) + len(cd45_neg) == len(cd45_low)
    assert len(recipient) + len(donor) + len(cd45_neg) == len(live)

    # Mac1 within recipient
    mac1_pos_recipient = recipient[recipient["FL2-H"] >= mac1_cutoff]
    mac1_neg_recipient = recipient[recipient["FL2-H"] < mac1_cutoff]
    mac1_pos_recipient_percent = pct(mac1_pos_recipient, recipient)
    assert len(mac1_pos_recipient) + len(mac1_neg_recipient) == len(recipient)

    # Mac1 within donor
    mac1_pos_donor = donor[donor["FL2-H"] >= mac1_cutoff]
    mac1_neg_donor = donor[donor["FL2-H"] < mac1_cutoff]
    mac1_pos_donor_percent = pct(mac1_pos_donor, donor)
    assert len(mac1_pos_donor) + len(mac1_neg_donor) == len(donor)

    # Donor share of all myeloid cells (donor myeloid chimerism)
    all_mac1_pos = len(mac1_pos_recipient) + len(mac1_pos_donor)
    donor_myeloid_chimerism = np.nan if all_mac1_pos == 0 else 100 * len(mac1_pos_donor) / all_mac1_pos

    summary = {"events": len(df_sample),
               "live": len(live),
               "live_pct": live_percent,
               "recipient": len(recipient),
               "recipient_pct": recipient_percent,
               "donor": len(donor),
               "donor_pct": donor_percent,
               "cd45_neg": len(cd45_neg),
               "cd45_neg_pct": cd45_neg_percent,
               "mac1_pos_recipient": len(mac1_pos_recipient),
               "mac1_neg_recipient": len(mac1_neg_recipient),
               "mac1_pos_recipient_pct": mac1_pos_recipient_percent,
               "mac1_pos_donor": len(mac1_pos_donor),
               "mac1_neg_donor": len(mac1_neg_donor),
               "mac1_pos_donor_pct": mac1_pos_donor_percent,
               "donor_myeloid_chimerism_pct": donor_myeloid_chimerism,
               "pi_cutoff": pi_cutoff,
               "cd45_2_cutoff": cd45_2_cutoff,
               "cd45_1_cutoff": cd45_1_cutoff,
               "mac1_cutoff": mac1_cutoff
               }
    return summary


def assign_labels(df_sample, pi_cutoff, cd45_2_cutoff, cd45_1_cutoff, mac1_cutoff):
    """Give every event the name of the group it falls into, using the same gates as apply_gating."""
    dead = df_sample["FL3-H"] >= pi_cutoff
    recipient = ~dead & (df_sample["FL4-H"] >= cd45_2_cutoff)
    donor = ~dead & (df_sample["FL4-H"] < cd45_2_cutoff) & (df_sample["FL1-H"] >= cd45_1_cutoff)
    cd45_neg = ~dead & (df_sample["FL4-H"] < cd45_2_cutoff) & (df_sample["FL1-H"] < cd45_1_cutoff)
    mac1_pos = df_sample["FL2-H"] >= mac1_cutoff

    conditions = [dead,
                  recipient & mac1_pos, recipient & ~mac1_pos,
                  donor & mac1_pos, donor & ~mac1_pos,
                  cd45_neg]
    names = ["dead", "recipient Mac1+", "recipient Mac1-",
             "donor Mac1+", "donor Mac1-", "CD45-negative"]
    return np.select(conditions, names, default="unassigned")


def create_table(folder="data/FlowCAP-I HSCT", out_file="results/hsct_gating_summary.csv", cutoffs=CUTOFFS):
    """Gate every FCS file in a folder with the same cut-offs and save one row per file."""
    sorted_files = sorted(pathlib.Path(folder).glob("*.fcs"))
    all_summaries = []
    for file in sorted_files:
        _, temp_df = inspect_fcs_file(file)
        temp_summary = apply_gating(temp_df, **cutoffs)
        temp_summary["file"] = file.name
        all_summaries.append(temp_summary)

    table_df = pd.DataFrame(all_summaries)
    table_df = table_df[["file"] + [c for c in table_df.columns if c != "file"]]   # file name as first column
    table_df.to_csv(out_file, index=False)
    return table_df