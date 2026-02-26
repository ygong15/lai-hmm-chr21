import csv
import sys
import subprocess
from pathlib import Path

VCF = "data/1000G_chr21_pruned.vcf.gz"
EMISS = "data/emissions.tsv"
GENO = "data/genotype.tsv"
PRED = "data/pred.tsv"
OUTFILE = "data/real_eval.tsv"

AFR = [
    "HG02983","HG03054","NA20126","HG03559","HG03268",
    "HG03280","HG03130","HG02570","NA19834","NA19320"
]

EUR = [
    "HG01784","HG00099","NA07000","HG00106","NA20585",
    "HG00311","HG00250","HG00178","HG00736","NA12716"
]

PY = sys.executable


def run_cmd(args):
    """Run a command as a list (no shell), exit on failure."""
    r = subprocess.run(args)
    if r.returncode != 0:
        print("Command failed:", " ".join(args))
        sys.exit(r.returncode)


def compute_metrics(pred_fp):
    n = 0
    a1 = a2 = 0
    f1 = f2 = 0
    p1 = p2 = None

    with open(pred_fp, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            s1 = row["POP1"]
            s2 = row["POP2"]

            a1 += (s1 == "A")
            a2 += (s2 == "A")

            if p1 is not None and s1 != p1:
                f1 += 1
            if p2 is not None and s2 != p2:
                f2 += 1

            p1, p2 = s1, s2
            n += 1

    if n == 0:
        raise RuntimeError(f"No rows found in {pred_fp}")
    return a1 / n, a2 / n, f1, f2


def run_sample(group, sample):
    print(f"Running {group} {sample}")

    run_cmd([PY, "scripts/make_inputs.py", VCF, sample, "AFR_AF", "EUR_AF", EMISS, GENO])
    run_cmd([PY, "main.py"])

    return compute_metrics(PRED)


def main():
    Path("data").mkdir(exist_ok=True)

    with open(OUTFILE, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["group", "sample", "hap1_A", "hap2_A", "flip1", "flip2"])

        for s in AFR:
            h1, h2, f1, f2 = run_sample("AFR", s)
            writer.writerow(["AFR", s, f"{h1:.6f}", f"{h2:.6f}", f1, f2])

        for s in EUR:
            h1, h2, f1, f2 = run_sample("EUR", s)
            writer.writerow(["EUR", s, f"{h1:.6f}", f"{h2:.6f}", f1, f2])

    summary = {"AFR": [0.0, 0.0, 0], "EUR": [0.0, 0.0, 0]}  # sum_h1, sum_h2, count
    with open(OUTFILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            g = row["group"]
            summary[g][0] += float(row["hap1_A"])
            summary[g][1] += float(row["hap2_A"])
            summary[g][2] += 1

    print("\nGroup means:")
    for g in ("AFR", "EUR"):
        h1 = summary[g][0] / summary[g][2]
        h2 = summary[g][1] / summary[g][2]
        print(f"{g}: mean hap1_A={h1:.4f}, mean hap2_A={h2:.4f} (n={summary[g][2]})")

    print(f"\nWrote {OUTFILE}")


if __name__ == "__main__":
    main()