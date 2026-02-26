#!/usr/bin/env python3
import argparse
import csv
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lai_hmm.model import emission_logprob


def load_positions_and_calls(emissions_fp: str, genotype_fp: str):
    rows = []
    with open(emissions_fp, newline="") as efile, open(genotype_fp, newline="") as gfile:
        e_reader = csv.DictReader(efile, delimiter="\t")
        g_reader = csv.DictReader(gfile, delimiter="\t")

        for e_row, g_row in zip(e_reader, g_reader):
            if e_row["POS"] != g_row["POS"] or e_row["CHR"] != g_row["CHR"]:
                raise ValueError(
                    f"Inputs misaligned: emissions={e_row['CHR']}:{e_row['POS']} "
                    f"genotype={g_row['CHR']}:{g_row['POS']}"
                )

            chrom = e_row["CHR"]
            pos = int(e_row["POS"])
            pA = float(e_row["pA_alt"])
            pB = float(e_row["pB_alt"])
            hap1 = g_row["hap1"]  # "0" or "1"
            hap2 = g_row["hap2"]  # "0" or "1"
            rows.append((chrom, pos, hap1, hap2, pA, pB))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emissions", default="data/emissions.tsv")
    ap.add_argument("--genotype", default="data/sim_genotype.tsv")
    ap.add_argument("--output", default="data/sim_pred_baseline.tsv")
    args = ap.parse_args()

    rows = load_positions_and_calls(args.emissions, args.genotype)
    states = ["AA", "AB", "BA", "BB"]

    with open(args.output, "w", newline="") as out:
        out.write("CHR\tPOS\tPOP1\tPOP2\n")
        for chrom, pos, hap1, hap2, pA, pB in rows:
            # pick best state independently at each SNP
            best_state = max(
                states,
                key=lambda s: emission_logprob(s, hap1, hap2, pA, pB)
            )
            out.write(f"{chrom}\t{pos}\t{best_state[0]}\t{best_state[1]}\n")

    print("Wrote:", args.output)


if __name__ == "__main__":
    main()