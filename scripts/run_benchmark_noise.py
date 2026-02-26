import subprocess
import pandas as pd
import sys

PYTHON = sys.executable

flip_rates = [0.0, 0.005, 0.01]
rhos = [5e-8, 1e-7, 5e-7]
seeds = range(1, 11)

rows = []

for flip in flip_rates:
    for rho in rhos:
        for seed in seeds:

            # simulate
            subprocess.run([
                PYTHON, "scripts/simulate_admixed.py",
                "--emissions", "data/emissions.tsv",
                "--out_genotype", "data/sim_genotype.tsv",
                "--out_truth", "data/sim_truth.tsv",
                "--rho", str(rho),
                "--seed", str(seed),
                "--flip_rate", str(flip)
            ], check=True, stdout=subprocess.DEVNULL)

            # run HMM
            subprocess.run([
                PYTHON, "main.py",
                "--rho", str(rho),
                "--emissions", "data/emissions.tsv",
                "--genotype", "data/sim_genotype.tsv",
                "--output", "data/sim_pred.tsv"
            ], check=True, stdout=subprocess.DEVNULL)

            # evaluate
            out = subprocess.check_output([
                PYTHON, "scripts/benchmark_metrics.py",
                "--truth", "data/sim_truth.tsv",
                "--pred", "data/sim_pred.tsv"
            ]).decode().strip().splitlines()

            acc = float(out[0].split()[1])
            true_sw = int(out[1].split()[-1])
            pred_sw = int(out[2].split()[-1])
            skew = float(out[3].split()[-1]) if true_sw > 0 else None

            rows.append({
                "flip_rate": flip,
                "rho": rho,
                "seed": seed,
                "accuracy": acc,
                "true_sw": true_sw,
                "pred_sw": pred_sw,
                "skew": skew
            })

df = pd.DataFrame(rows)
df.to_csv("results/bench_noise.tsv", sep="\t", index=False)

print("Saved results to results/bench_noise.tsv")
print()
print(df.groupby(["flip_rate","rho"])[["accuracy","skew"]].mean().round(6))