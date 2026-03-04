#!/usr/bin/env bash
set -euo pipefail

# Detect python interpreter
PYTHON=$(command -v python3 || command -v python)
if [[ -z "$PYTHON" ]]; then
  echo "ERROR: python not found."
  exit 1
fi

mkdir -p results

echo "=== [1/2] Preprocess: Precheck ==="
if [[ ! -f "data/1000G_chr21_pruned.vcf.gz" ]]; then
  echo "ERROR: Missing data/1000G_chr21_pruned.vcf.gz"
  exit 1
fi

echo "=== [2/2] Preprocess: Build emissions.tsv (if missing) ==="
SAMPLE="${SAMPLE:-HG00099}"
POPA_KEY="${POPA_KEY:-AFR_AF}"
POPB_KEY="${POPB_KEY:-EUR_AF}"

if [[ ! -f "data/emissions.tsv" ]]; then
  echo "Generating data/emissions.tsv using sample=$SAMPLE ..."
  $PYTHON scripts/make_inputs.py \
    data/1000G_chr21_pruned.vcf.gz \
    "$SAMPLE" \
    "$POPA_KEY" \
    "$POPB_KEY" \
    data/emissions.tsv \
    data/genotype.tsv
else
  echo "Found data/emissions.tsv (skipping make_inputs)"
fi

echo "Sanity check: emissions.tsv rows"
$PYTHON - <<'PY'
import pandas as pd
df = pd.read_csv("data/emissions.tsv", sep="\t")
print("rows:", len(df), "cols:", list(df.columns))
print(df.head(2).to_string(index=False))
PY

echo "=== [1/4] Test 1: Simulation accuracy (ground truth) ==="
$PYTHON scripts/simulate_admixed.py \
  --emissions data/emissions.tsv \
  --out_genotype data/sim_genotype.tsv \
  --out_truth data/sim_truth.tsv \
  --rho 1e-7 \
  --seed 1 \
  --flip_rate 0.0

$PYTHON main.py \
  --rho 1e-7 \
  --emissions data/emissions.tsv \
  --genotype data/sim_genotype.tsv \
  --output data/sim_pred.tsv

$PYTHON scripts/benchmark_metrics.py \
  --truth data/sim_truth.tsv \
  --pred data/sim_pred.tsv

echo "Plotting truth vs pred..."
$PYTHON scripts/plot_truth_vs_pred.py >/dev/null || true

echo "=== [2/4] Test 2: Robustness sweep (rho x flip_rate x seeds) ==="
$PYTHON scripts/run_benchmark_noise.py
$PYTHON scripts/plot_benchmark_noise.py >/dev/null || true

echo "=== [3/4] Test 3: Baseline comparison (no transitions) ==="
$PYTHON scripts/baseline_independent.py \
  --emissions data/emissions.tsv \
  --genotype data/sim_genotype.tsv \
  --output data/sim_pred_baseline.tsv

echo "--- HMM metrics ---"
$PYTHON scripts/benchmark_metrics.py \
  --truth data/sim_truth.tsv \
  --pred data/sim_pred.tsv

echo "--- Baseline metrics ---"
$PYTHON scripts/benchmark_metrics.py \
  --truth data/sim_truth.tsv \
  --pred data/sim_pred_baseline.tsv

echo "=== [4/4] Test 4: Real-data plausibility (AFR vs EUR top-10) ==="
$PYTHON scripts/find_pure_samples.py
$PYTHON scripts/run_real_eval.py
$PYTHON scripts/plot_real_eval.py >/dev/null || true

echo ""
echo "All tests completed."
echo "Outputs:"
echo "  data/sim_genotype.tsv"
echo "  data/sim_truth.tsv"
echo "  data/sim_pred.tsv"
echo "  data/sim_pred_baseline.tsv"
echo "  results/bench_noise.tsv"
echo "  results/*.png (plotting scripts saved figures)"
echo "  data/real_eval.tsv (from real data eval)"
