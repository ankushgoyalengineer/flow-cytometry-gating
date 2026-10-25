import argparse
import time
from pathlib import Path

from flowtool.gating import CUTOFFS, create_table
from flowtool.validation import validate_against_labels

LIKELY_CONTROLS = ["027.fcs", "028.fcs", "029.fcs", "030.fcs"]
LABELLED_SAMPLES = {5: "005.fcs", 6: "006.fcs", 9: "009.fcs", 12: "012.fcs"}


def main():
    parser = argparse.ArgumentParser(description="Gate FCS files and summarise the populations.")
    parser.add_argument("--folder", required=True, help="Folder containing the FCS files")
    parser.add_argument("--out", default="results", help="Folder for the output tables")
    parser.add_argument("--labels", help="Optional: ks expert-label CSV, to validate the gating")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(exist_ok=True)

    start = time.perf_counter()

    # 1. Gate every file
    table = create_table(folder=args.folder, out_file=out / "hsct_gating_summary.csv", cutoffs=CUTOFFS)
    table["likely_control"] = table["file"].isin(LIKELY_CONTROLS)
    table.to_csv(out / "hsct_gating_summary.csv", index=False)
    samples = table[~table["likely_control"]]
    print(f"Gated {len(table)} files ({len(samples)} samples, {len(table) - len(samples)} likely controls)")
    print(samples.filter(like="_pct").median().round(2).to_string())

    # 2. Optional validation against expert labels
    if args.labels:
        scores = validate_against_labels(args.labels, args.folder, LABELLED_SAMPLES, CUTOFFS)
        scores.to_csv(out / "hsct_validation_scores.csv", index=False)
        valid = scores[scores["note"] == ""]
        print(f"\nF-measure across {len(valid)} populations: "
              f"min {valid['f_measure'].min():.3f}, max {valid['f_measure'].max():.3f}, "
              f"median {valid['f_measure'].median():.3f}")

    print(f"\nFinished in {time.perf_counter() - start:.1f} seconds")


if __name__ == "__main__":
    main()