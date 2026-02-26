import argparse
import csv
import matplotlib.pyplot as plt

def read_rows(fp):
    rows = []
    with open(fp, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append({
                "group": row["group"],
                "sample": row["sample"],
                "hap1_A": float(row["hap1_A"]),
                "hap2_A": float(row["hap2_A"]),
                "flip1": int(row["flip1"]),
                "flip2": int(row["flip2"]),
            })
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", default="data/real_eval.tsv")
    ap.add_argument("--out_scatter", default="results/real_eval_scatter.png")
    ap.add_argument("--out_box", default="results/real_eval_box.png")
    ap.add_argument("--title", default="Chr21 real-data validation (A=AFR)")
    args = ap.parse_args()

    rows = read_rows(args.infile)
    afr = [r for r in rows if r["group"] == "AFR"]
    eur = [r for r in rows if r["group"] == "EUR"]

    # --- Scatter: hap1_A vs hap2_A ---
    plt.figure()
    plt.scatter([r["hap1_A"] for r in afr], [r["hap2_A"] for r in afr], label="AFR-like")
    plt.scatter([r["hap1_A"] for r in eur], [r["hap2_A"] for r in eur], label="EUR-like")
    plt.xlabel("hap1 A fraction")
    plt.ylabel("hap2 A fraction")
    plt.title(args.title + " (hap1 vs hap2)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out_scatter, dpi=200)
    plt.close()
    print(f"Wrote: {args.out_scatter}")

    afr_mean = [(r["hap1_A"] + r["hap2_A"]) / 2 for r in afr]
    eur_mean = [(r["hap1_A"] + r["hap2_A"]) / 2 for r in eur]

    plt.figure()
    plt.boxplot([afr_mean, eur_mean], labels=["AFR-like", "EUR-like"])
    plt.ylabel("mean A fraction (avg of hap1,hap2)")
    plt.title(args.title + " (mean A distribution)")
    plt.tight_layout()
    plt.savefig(args.out_box, dpi=200)
    plt.close()
    print(f"Wrote: {args.out_box}")

if __name__ == "__main__":
    main()