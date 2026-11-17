# Benchmarking Bayesian phylogenetic inference from nucleotide sequence data: code and data

Code and per-run results for the article

> Chung Y. Benchmarking Bayesian phylogenetic inference from nucleotide sequence data. *Genes & Genomics* (submitted).

The article compares Markov chain Monte Carlo, run with BEAST X, and annealed sequential Monte Carlo, run with the
reference implementation of Wang et al. (2020), on simulated alignments and on the M1809 alignment of 59 mammalian
mitochondrial sequences.

Two things can be reproduced from this repository:

1. **The figures, tables and quoted numbers of the article**, from the per-run result tables in `pipeline/tables/`.
   This needs only Python and takes about a minute. See [Reproducing the figures and tables](#reproducing-the-figures-and-tables).
2. **The runs themselves**, with the scripts in `pipeline/`. This needs BEAST X, BEAGLE, Seq-Gen and the annealed SMC
   program, and several thousand CPU-hours. See [Re-running the experiments](#re-running-the-experiments).

## Contents

| Path | Contents |
|---|---|
| `pipeline/` | The scripts that generated the data and ran and scored both programs, as they were run on the server |
| `pipeline/tables/` | One row per run: accuracy, CPU time, wall-clock time (the input of everything in `analysis/`) |
| `analysis/` | Scripts that turn `pipeline/tables/` into the figures, tables and numbers of the article |
| `output/figures/` | Figs. 1–5 and S1–S4 as produced by `analysis/make_figures.py` (PDF and 600-dpi PNG) |
| `output/tables/` | Tables S1–S9 and the summary tables as produced by `analysis/` |

## Reproducing the figures and tables

Requirements: Python 3 with NumPy, SciPy and Matplotlib (`pip install -r requirements.txt`).

```bash
cd analysis
python3 make_summaries.py        # means, SDs, paired differences  -> output/tables/
python3 make_figures.py          # Figs. 1-5, S1-S4                -> output/figures/
python3 make_tables.py           # Tables S1-S9                    -> output/tables/
python3 check_paper_numbers.py   # recomputes the numbers quoted in the text; prints "0 failed"
```

| Item in the article | Script | Input in `pipeline/tables/` |
|---|---|---|
| Fig. 1, Fig. 2, Fig. S1 (error and posterior support against N) | `make_figures.py` | `equal_n_runs.tsv`, `alpha_rerun.csv` |
| Fig. 3, Fig. 4 (cost; error against cost), Figs. S2–S3 | `make_figures.py` | `alpha_rerun.csv` |
| Fig. 5 (parallel scaling) | `make_figures.py` | `scaling_runs.tsv` |
| Fig. S4 (M1809) | `make_figures.py` | `emp_converge_EMP_M1809.json`, `emp_converge_EMPA5_M1809.json` |
| Table S1 (N = 100) | `make_tables.py` | `equal_n_runs.tsv`, `alpha_rerun.csv` |
| Table S2 (N = 2,000) | `make_tables.py` | `equal_n_runs.tsv`, `wall_mean_sd.csv` |
| Table S3 (four values of α) | `make_tables.py` | `alpha_pilot.csv`, `alpha_rerun.csv` |
| Table S4 (M1809) | `make_tables.py` | `emp_converge_EMP_M1809.json`, `emp_converge_EMPA5_M1809.json` |
| Table S5 (paired differences) | `make_tables.py` | `equal_n_runs.tsv`, `alpha_rerun.csv` |
| Table S6 (error and support across N) | `make_tables.py` | `equal_n_runs.tsv`, `alpha_rerun.csv` |
| Table S7 (CPU time, wall-clock time, parallel scaling) | `make_tables.py` | `alpha_rerun.csv`, `scaling_runs.tsv` |
| Table S8 (rate variation among lineages) | `make_tables.py` | `equal_n_runs.tsv`, `alpha_rerun.csv` |
| Table S9 (chains of 20,000,000 states) | `make_tables.py` | `conv_check.csv` |
| Numbers quoted in the abstract, Results and Online Resource 1 (Results S1) | `check_paper_numbers.py`, `make_summaries.py` | all of the above and `conv_check.csv` |

Tables 1 and 2 of the article (the simulation design and the practical guidance) contain no computed values.
The figures of the article are lettered in the Inter typeface; if it is not installed, Matplotlib's default typeface
is used and only the lettering differs.

### The result tables

| File | Rows | Contents |
|---|---|---|
| `equal_n_runs.tsv` | 3,000 | BEAST X and annealed SMC at α = 0.99: six conditions × 50 replicates × five values of N. Normalized Robinson–Foulds distance of the majority-rule consensus (`nrf_cons`) and of the MCC tree (`nrf_mcc`) from the true tree, posterior probability of the true topology (`p_true`), CPU seconds (`cpu_s`) |
| `alpha_rerun.csv` | 420 | Annealed SMC at α = 0.99999: N = 100 on replicates 1–50, N = 200–2,000 on replicates 1–5, with CPU seconds, wall-clock seconds at 32 threads and the number of annealing steps. Columns prefixed `bx_` and `a099_` give the BEAST X run and the α = 0.99 run on the same alignment and N; their wall-clock columns are filled for replicates 1–5 |
| `alpha_pilot.csv` | 69 | Annealed SMC at α = 0.999, 0.9999 and 0.99999 on a subset of alignments (Table S3) |
| `wall_mean_sd.csv` | 120 | Mean and SD of wall-clock and CPU time over replicates 1–5; subset `b` is annealed SMC at 32 threads and BEAST X at 8 BEAGLE threads |
| `scaling_runs.tsv` | 18 | Thread ladders: annealed SMC (N = 200) on one S4 and one S3 alignment at 1–64 and 1–32 threads, BEAST X (2,000,000 states) on the S4 alignment at 1–16 BEAGLE threads |
| `conv_check.csv` | 15 | BEAST X chains of 20,000,000 states on replicates 1–5 of S3, S4 and R2, beside the N = 2,000 and N = 100 runs on the same alignments |
| `emp_converge_EMP_M1809.json` (`.csv`) | 5 rungs | M1809, BEAST X and annealed SMC at α = 0.99 at N = 1,000–50,000: agreement between seeds within each program and between the programs |
| `emp_converge_EMPA5_M1809.json` (`.csv`) | 1 rung | The same at N = 1,000 with annealed SMC at α = 0.99999 |

Conditions: S1 20 taxa × 1,000 bp; S2 20 × 5,000; S3 50 × 1,000; S4 50 × 5,000 (strict clock); R1 20 × 1,000 and
R2 50 × 1,000 (branch rates multiplied by independent mean-one lognormal factors, σ = 0.3).

## Re-running the experiments

### Software

The runs of the article used Ubuntu 20.04 on a 112-core server with BEAST X v10.5.0, BEAGLE v4.0.1 (CPU, SSE),
Seq-Gen 1.3.4, OpenJDK 17, Python 3.8 with DendroPy 4.6.4 and Biopython 1.83, and the annealed SMC program built from
<https://github.com/liangliangwangsfu/annealedSMC> (commit `bfab06d`) with the two patches below.

Seq-Gen must be version 1.3.4: version 1.3.5 writes FASTA by default, which `simulate.py` does not read.

Set the two paths at the top of `pipeline/env.sh` (`BEAST_HOME`, `ASMC_HOME`), then start every session with

```bash
cd pipeline && source env.sh
```

### Patches to the annealed SMC program

```bash
python3 patch_asmc.py --asmc-home $ASMC_HOME           # writes every final particle and its weight to disk
python3 patch_asmc_threads.py --asmc-home $ASMC_HOME   # makes the one-time initialization of the move list thread-safe
bash rebuild_asmc.sh
bash verify_thread_patch.sh     # patched runs at 1, 8 and 16 threads are byte-identical; runs that aborted unpatched now finish
bash verify_thread_patch2.sh    # single-thread output is byte-identical with and without the patch; unpatched 8-thread runs are not reproducible
```

`metrics.py` is the summarization routine applied to the output of both programs; `python3 test_metrics.py` runs its
validation cases.

### Steps

| Step | Command | Result |
|---|---|---|
| Simulate the alignments (300; seeds are a function of condition and replicate) | `bash gen_data.sh`, then `REP_START=11 REPS=50 bash gen_data.sh` | `data/<dataset>/` |
| Annealed SMC at α = 0.99, 50 replicates; thread ladders | `python3 sweep.py --phases scal,asmc` | `runs_asmc/`, `runs_beastx/`, run log |
| BEAST X, 50 replicates; wall-clock runs on replicates 1–5 | `python3 equal_n_sweep.py`, `python3 score_equal_n.py`, `python3 wall_mean_sd.py` | `tables/equal_n_runs.tsv`, `tables/wall_mean_sd.csv` |
| Annealed SMC at α = 0.99999 | `python3 alpha_rerun.py --reps 1-5`, `python3 alpha_rerun.py --reps 1-50 --ns 100`, `python3 alpha_rerun.py --score-only --reps 1-50` | `tables/alpha_rerun.csv` |
| Annealed SMC at intermediate α | `python3 alpha_pilot.py`, `python3 alpha_pilot.py --score-only` | `tables/alpha_pilot.csv` |
| BEAST X chains of 20,000,000 states | `python3 conv_rerun.py`, `python3 conv_rerun.py --score-only` | `tables/conv_check.csv` |
| M1809 | `bash get_empirical.sh`, `python3 emp_converge.py`, `python3 emp_rung_score.py --dataset EMP_M1809`, `python3 emp_alpha.py`, `python3 emp_alpha.py --score` | `tables/emp_converge_*.json` |

Each driver script is resumable (finished runs are skipped) and describes its options and its scheduling rules in its
header; `--dry-run` lists the jobs without starting them. `run_beastx.sh` and `run_asmc.sh` run one analysis each and
are called by the drivers. Both programs are run under HKY with κ = 2 and equal base frequencies fixed, the model the
alignments were simulated under (`seq-gen -mHKY -t1`: a transition/transversion ratio of 1 is κ = 2 at equal
frequencies).

Notes:

- `sweep.py` also defines the phases `beastx` and `extra`, which belong to an earlier design (BEAST X at fixed chain
  lengths) and are not used in the article; `equal_n_sweep.py` likewise defines a wall-clock subset `a` (both programs
  at 8 threads) that the article does not report.
- The simulated alignments, the M1809 alignment and the raw output of the runs (posterior tree samples) are not stored
  here. The alignments are regenerated by `gen_data.sh`; M1809 is downloaded by `get_empirical.sh` from the repository
  of the annealed SMC authors.
- `pipeline/tables/alpha_rerun.csv` holds the 420 runs at α = 0.99999 that the article reports.
- In the usage comments of the scripts the working directory is written as `pipeline`; otherwise the scripts are
  unchanged from the versions that produced the results.

## Reference

Wang L, Wang S, Bouchard-Côté A (2020) An annealed sequential Monte Carlo method for Bayesian phylogenetics.
Syst Biol 69:155–183. https://doi.org/10.1093/sysbio/syz028
