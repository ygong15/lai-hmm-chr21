import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/bench_noise.tsv", sep="\t")

# Remove runs with true_sw == 0 (skew undefined)
df = df[df["true_sw"] > 0]

# Compute mean/std
summary = (
    df.groupby(["flip_rate", "rho"])
      .agg({"accuracy": ["mean", "std"],
            "skew": ["mean", "std"]})
      .reset_index()
)

summary.columns = [
    "flip_rate", "rho",
    "acc_mean", "acc_std",
    "skew_mean", "skew_std"
]

# Plot 1: Accuracy vs rho
plt.figure(figsize=(7,5))

for flip in sorted(summary["flip_rate"].unique()):
    sub = summary[summary["flip_rate"] == flip]
    plt.errorbar(
        sub["rho"],
        sub["acc_mean"],
        yerr=sub["acc_std"],
        marker="o",
        label=f"flip_rate={flip}"
    )

plt.xscale("log")
plt.xlabel("Recombination rate (rho)")
plt.ylabel("Accuracy")
plt.title("HMM Accuracy vs Recombination Rate")
plt.legend()
plt.tight_layout()
plt.savefig("results/accuracy_vs_rho.png", dpi=300)

# Plot 2: Skew vs rho
plt.figure(figsize=(7,5))

for flip in sorted(summary["flip_rate"].unique()):
    sub = summary[summary["flip_rate"] == flip]
    plt.errorbar(
        sub["rho"],
        sub["skew_mean"],
        yerr=sub["skew_std"],
        marker="o",
        label=f"flip_rate={flip}"
    )

plt.xscale("log")
plt.xlabel("Recombination rate (rho)")
plt.ylabel("Switch Skew (pred/true)")
plt.title("Switch Skew vs Recombination Rate")
plt.legend()
plt.tight_layout()
plt.savefig("results/skew_vs_rho.png", dpi=300)
