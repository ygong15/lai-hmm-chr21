import argparse
import csv

from lai_hmm.model import emission_logprob
from lai_hmm.viterbi import viterbi


def load_inputs(emissions_fp: str, genotype_fp: str):
    emissions = []
    positions = []

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
            pos_bp = int(e_row["POS"])
            positions.append((chrom, pos_bp))

            pA = float(e_row["pA_alt"])
            pB = float(e_row["pB_alt"])
            hap1 = g_row["hap1"]
            hap2 = g_row["hap2"]

            state_probs = {}
            for state in ["AA", "AB", "BA", "BB"]:
                state_probs[state] = emission_logprob(state, hap1, hap2, pA, pB)

            emissions.append(state_probs)

    return positions, emissions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=1e-8)
    ap.add_argument("--emissions", default="data/emissions.tsv")
    ap.add_argument("--genotype", default="data/genotype.tsv")
    ap.add_argument("--output", default="data/pred.tsv")
    args = ap.parse_args()

    positions, emissions = load_inputs(args.emissions, args.genotype)
    positions_bp = [pos for (_chr, pos) in positions]

    path = viterbi(emissions, positions_bp, args.rho)

    with open(args.output, "w", newline="") as out:
        out.write("CHR\tPOS\tPOP1\tPOP2\n")
        for (chrom, pos_bp), state in zip(positions, path):
            out.write(f"{chrom}\t{pos_bp}\t{state[0]}\t{state[1]}\n")


if __name__ == "__main__":
    main()