import gzip
import sys

def parse_info(info_str: str) -> dict:
    d = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d

def main(vcf_gz, sample, popA_key, popB_key, emissions_out, genotype_out):
    with gzip.open(vcf_gz, "rt") as f, \
         open(emissions_out, "w") as eout, \
         open(genotype_out, "w") as gout:

        sample_idx = None

        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                if sample not in samples:
                    raise SystemExit(f"Sample {sample} not found in VCF. Example: {samples[0]}")
                sample_idx = 9 + samples.index(sample)

                eout.write("CHR\tPOS\tREF\tALT\tpA_alt\tpB_alt\n")
                gout.write("CHR\tPOS\thap1\thap2\n")
                continue

            # data line
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parse_info(parts[7])

            # keep simple SNPs only
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                continue
            if "," in alt:
                continue 

            if popA_key not in info or popB_key not in info:
                continue

            try:
                pA = float(info[popA_key].split(",")[0])
                pB = float(info[popB_key].split(",")[0])
            except ValueError:
                continue

            # genotype for target sample
            gt = parts[sample_idx]
            if gt == "." or gt == "./.":
                continue

            # expect phased 0|1
            if "|" not in gt:
                continue
            a, b = gt.split("|", 1)
            if a not in ("0", "1") or b not in ("0", "1"):
                continue

            eout.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{pA}\t{pB}\n")
            gout.write(f"{chrom}\t{pos}\t{a}\t{b}\n")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Usage:\n"
            "  make_inputs.py <vcf.gz> <sample> <popA_key> <popB_key> <emissions.tsv> <genotype.tsv>\n\n"
            "Example (AFR vs EUR):\n"
            "  make_inputs.py 1000G_chr21_pruned.vcf.gz HG00100 AFR_AF EUR_AF emissions.tsv genotype.tsv",
            file=sys.stderr
        )
        sys.exit(2)

    vcf_gz, sample, popA_key, popB_key, emissions_out, genotype_out = sys.argv[1:]
    main(vcf_gz, sample, popA_key, popB_key, emissions_out, genotype_out)