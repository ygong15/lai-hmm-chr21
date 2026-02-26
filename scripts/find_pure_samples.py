import gzip
import math
import random

VCF = "data/1000G_chr21_pruned.vcf.gz"
POP_A = "AFR_AF"   # A = AFR
POP_B = "EUR_AF"   # B = EUR

EPS = 1e-6
MAX_SITES = 2000

def clamp(p):
    return max(EPS, min(1 - EPS, p))

def parse_info(s):
    d = {}
    for x in s.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d

def logprob_gt(gt, p_alt):
    a,b = gt.split("|")
    p = clamp(p_alt)
    def lp(h):
        return math.log10(p if h=="1" else (1-p))
    return lp(a) + lp(b)

def main():
    with gzip.open(VCF, "rt") as f:
        samples = None
        scores = None  # per sample: sum(logP_B - logP_A)

        n_used = 0

        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                scores = [0.0] * len(samples)
                continue

            if samples is None:
                continue

            parts = line.rstrip("\n").split("\t")
            ref, alt = parts[3], parts[4]
            if len(ref)!=1 or len(alt)!=1 or "," in alt:
                continue

            info = parse_info(parts[7])
            if POP_A not in info or POP_B not in info:
                continue
            if info[POP_A] == "." or info[POP_B] == ".":
                continue

            pA = float(info[POP_A])
            pB = float(info[POP_B])

            n_used += 1
            if n_used > MAX_SITES:
                break

            gts = parts[9:]
            for i, gt in enumerate(gts):
                if "|" not in gt:
                    continue
                # score positive = more EUR-like (B), negative = more AFR-like (A)
                scores[i] += (logprob_gt(gt, pB) - logprob_gt(gt, pA))

    ranked = sorted(zip(samples, scores), key=lambda x: x[1])

    print(f"Used sites: {n_used}")
    print("\nMost AFR-like (lowest score):")
    for s, sc in ranked[:10]:
        print(s, sc)

    print("\nMost EUR-like (highest score):")
    for s, sc in ranked[-10:][::-1]:
        print(s, sc)

if __name__ == "__main__":
    main()