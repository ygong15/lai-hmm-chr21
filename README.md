# Local Ancestry Inference with HMM (chr21)

In this project, we built a simplified Hidden Markov Model (HMM) implementation for **Local Ancestry Inference (LAI)** on chromosome 21.

This tool does the following 4 things:
- Builds emission probabilities from population allele frequencies  
- Uses recombination-aware transitions  
- Runs Viterbi decoding to infer ancestry blocks  
- Benchmarks accuracy, switch detection, and robustness to noise  

The repository includes simulation, benchmarking, and real-data validation tools.

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

# Installation

## Requirements

- Python 3.9+
- numpy
- pandas
- matplotlib

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

# Troubleshooting

If you see:

```
FileNotFoundError: emissions.tsv
```

Make sure `make_inputs.py` has been run (this is handled automatically in `run_all_tests.sh`).

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

# Tests Included

We built 4 tests to evaluate different aspects of the method.

---

## 1. Simulated Admixed Genome (Core Accuracy Test)

**Purpose:**  
Measure whether the HMM correctly recovers ancestry blocks when the true ancestry is known

**What it does:**
- Simulates an admixed genome  
- Runs HMM inference  
- Compares prediction to ground truth  
- Computes:
  - Accuracy  
  - True vs predicted switch counts  
  - Ancestry skew  

Output files:

```
data/sim_pred.tsv        # predicted ancestry from HMM
data/sim_truth.tsv       # ground truth ancestry from simulation
results/truth_vs_pred.png  # visualization of prediction vs truth
```

---

## 2. Baseline Comparison (Independent Sites)

**Purpose:**  
Show that modeling ancestry transitions improves prediction compared to independent per-site classification

**What it does:**
- Runs a naive per-site classifier  
- Ignores recombination  
- Compares to HMM results  

Output files:

```
data/sim_pred.tsv           # HMM predictions
data/sim_pred_baseline.tsv  # baseline per-site predictions
```
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
- Aggregates performance metrics across runs

Outputs:

```
results/bench_noise.tsv        # aggregated benchmark results
results/accuracy_vs_rho.png    # accuracy vs recombination rate
results/skew_vs_rho.png        # ancestry skew vs recombination rate
```

---

## 4. Real Data Sanity Check (Pure AFR / EUR Samples)

**Purpose:**  
Verify behavior on real 1000G samples that are expected to have mostly single ancestry to confirm that the model does not artificially hallucinate switches or mixed ancestry on pure samples

**What it does:**
- Runs inference on known AFR and EUR samples  
- Computes predicted ancestry proportions and switch counts

Outputs:

```
data/real_eval.tsv           # summary statistics for real samples
data/real_eval_scatter.png   # ancestry proportion visualization
data/real_eval_box.png       # switch frequency comparison
```
