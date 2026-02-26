# Local Ancestry Inference with HMM (chr21)

A simplified Hidden Markov Model (HMM) implementation for **Local Ancestry Inference (LAI)** on chromosome 21.

This tool:
- Builds emission probabilities from population allele frequencies  
- Uses recombination-aware transitions  
- Runs Viterbi decoding to infer ancestry blocks  
- Benchmarks accuracy, switch detection, and robustness to noise  

The repository includes simulation, benchmarking, and real-data validation tools.

---

# Quick Start (Run Everything)

To reproduce all tests and visualizations:

```bash
bash run_all_tests.sh
```

This will:

1. Generate required inputs  
2. Run simulated benchmarks  
3. Compare to a baseline model  
4. Perform robustness analysis (noise + recombination rate)  
5. Run real-data sanity checks  
6. Produce evaluation plots  

All results will be written to:

```
data/
results/
```

---

# Installation

## Requirements

- Python 3.9+
- numpy
- pandas
- matplotlib

---

# Repository Structure

```
main.py                  - HMM CLI runner
model.py                 - Transition + emission logic
viterbi.py               - Viterbi decoding

scripts/
  make_inputs.py         - Generate emissions + genotype
  simulate_admixed.py    - Create simulated admixed genome
  benchmark_metrics.py   - Accuracy + switch metrics
  baseline_independent.py- Independent-site baseline
  run_benchmark_noise.py - Grid benchmark (rho + noise)
  run_real_eval.py       - Evaluate on real “pure” samples
  plot_*.py              - Visualization scripts

run_all_tests.sh         - Run full pipeline
```

---

# Tests Included

We built several tests to evaluate different aspects of the method.

---

## 1. Simulated Admixed Genome (Core Accuracy Test)

**Purpose:**  
Measure whether the HMM correctly recovers ancestry blocks when the true ancestry is known.

**What it does:**
- Simulates an admixed genome  
- Runs HMM inference  
- Compares prediction to ground truth  
- Computes:
  - Accuracy  
  - True vs predicted switch counts  
  - Ancestry skew  

**Why this matters:**  
This directly measures whether the HMM works as intended under ideal conditions.

Output files:

```
data/sim_pred.tsv
data/sim_truth.tsv
results/truth_vs_pred.png
```

---

## 2. Baseline Comparison (Independent Sites)

**Purpose:**  
Show that modeling transitions improves performance.

**What it does:**
- Runs a naive per-site classifier  
- Ignores recombination  
- Compares to HMM results  

**Why this matters:**  
Demonstrates that ancestry blocks are not independent — transitions matter.

---

## 3. Robustness Benchmark (Noise + Recombination Sweep)

**Purpose:**  
Test how stable the method is under:
- Different recombination rates (rho)  
- Different genotype noise levels  

**What it does:**
- Runs many simulations across:
  - rho grid  
  - flip_rate grid  
  - multiple random seeds  
- Aggregates results  

Outputs:

```
results/bench_noise.tsv
results/accuracy_vs_rho.png
results/skew_vs_rho.png
```

**Why this matters:**  
Evaluates model robustness and stability.

---

## 4. Real Data Sanity Check (Pure AFR / EUR Samples)

**Purpose:**  
Verify behavior on real 1000G samples assumed to be ancestry-pure.

**What it does:**
- Runs inference on known AFR and EUR samples  
- Measures fraction of predicted ancestry  
- Evaluates switch frequency  

Outputs:

```
data/real_eval.tsv
data/real_eval_scatter.png
data/real_eval_box.png
```

**Why this matters:**  
Confirms the model does not artificially hallucinate switches or mixed ancestry on pure samples.

---

# Small Test Dataset

The repository includes a small pruned chromosome 21 dataset:

```
data/1000G_chr21_pruned.vcf.gz
data/1000G_chr21_pruned.vcf.gz.tbi
```

This dataset is sufficient to:
- Generate emissions  
- Simulate admixed genomes  
- Run all benchmarks  

No external downloads are required.

---

# Model Overview

The HMM models ancestry as hidden states:

States:
- POP1  
- POP2  

Transition probabilities:
- Controlled by recombination rate (rho)  
- Encourage long ancestry blocks with occasional switches  

Emission probabilities:
- Based on population-specific allele frequencies  
- Represent likelihood of genotype under each ancestry  

Decoding:
- Viterbi algorithm for most probable ancestry path    

---

# Troubleshooting

If you see:

```
FileNotFoundError: emissions.tsv
```

Make sure `make_inputs.py` has been run (this is handled automatically in `run_all_tests.sh`).