import argparse
import csv
import math
import random

EPS = 1e-6

def clamp(p: float) -> float:
    return max(EPS, min(1 - EPS, p))

def flip(bit: str) -> str:
    return "1" if bit == "0" else "0"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emissions", default="data/emissions.tsv")
    ap.add_argument("--out_genotype", default="data/sim_genotype.tsv")
    ap.add_argument("--out_truth", default="data/sim_truth.tsv")
    ap.add_argument("--rho", type=float, default=1e-8,
                    help="per-base recombination rate for r_i = 1-exp(-rho*Δbp)")
    ap.add_argument("--qA_hap1", type=float, default=0.5,
                    help="P(hap1 starts in ancestry A)")
    ap.add_argument("--qA_hap2", type=float, default=0.5,
                    help="P(hap2 starts in ancestry A)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flip_rate", type=float, default=0.0,
                    help="per-haplotype allele flip probability (genotyping error)")
    args = ap.parse_args()

    if not (0.0 <= args.flip_rate <= 1.0):
        raise SystemExit("--flip_rate must be in [0, 1].")

    random.seed(args.seed)

    # Read emissions rows
    rows = []
    with open(args.emissions, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            chrom = row["CHR"]
            pos = int(row["POS"])
            pA = clamp(float(row["pA_alt"]))
            pB = clamp(float(row["pB_alt"]))
            rows.append((chrom, pos, pA, pB))

    if len(rows) < 2:
        raise SystemExit("Not enough rows in emissions.tsv")

    # Initialize ancestries
    anc1 = "A" if random.random() < args.qA_hap1 else "B"
    anc2 = "A" if random.random() < args.qA_hap2 else "B"

    with open(args.out_genotype, "w", newline="") as gout, open(args.out_truth, "w", newline="") as tout:
        gout.write("CHR\tPOS\thap1\thap2\n")
        tout.write("CHR\tPOS\tTRUE1\tTRUE2\n")

        prev_pos = rows[0][1]

        for i, (chrom, pos, pA, pB) in enumerate(rows):
            if i > 0:
                delta = pos - prev_pos
                if delta < 0:
                    raise SystemExit(f"Positions not increasing: prev={prev_pos}, cur={pos}")

                # per-step recomb probability per haplotype
                r_i = 1.0 - math.exp(-args.rho * delta)
                r_i = max(0.0, min(1.0, r_i)) 

                # each hap switches independently
                if random.random() < r_i:
                    anc1 = "A" if anc1 == "B" else "B"
                if random.random() < r_i:
                    anc2 = "A" if anc2 == "B" else "B"
                prev_pos = pos

            # Sample alleles (0=REF, 1=ALT) based on ancestry-specific ALT freq
            p_alt1 = pA if anc1 == "A" else pB
            p_alt2 = pA if anc2 == "A" else pB

            hap1 = "1" if random.random() < p_alt1 else "0"
            hap2 = "1" if random.random() < p_alt2 else "0"

            # flip alleles independently
            if args.flip_rate > 0.0:
                if random.random() < args.flip_rate:
                    hap1 = flip(hap1)
                if random.random() < args.flip_rate:
                    hap2 = flip(hap2)

            gout.write(f"{chrom}\t{pos}\t{hap1}\t{hap2}\n")
            tout.write(f"{chrom}\t{pos}\t{anc1}\t{anc2}\n")

    print("Wrote:")
    print(" ", args.out_genotype)
    print(" ", args.out_truth)

if __name__ == "__main__":
    main()