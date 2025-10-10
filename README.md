# VDDP: Verifiable Distributed Differential Privacy

Welcome to the official implementation of [VDDP](https://arxiv.org/abs/2504.21752).

## Reproduction Instructions

Please follow the steps below carefully to reproduce the experimental results.

---

### Step 0: Fork (and star) this repository

This repository is designed to automatically generate raw experimental logs.
If you do **not** fork the repository, the code will attempt to push your logs to the original repository, which may fail or interfere with existing data.
To avoid this, please **fork** the repository before running any experiments. (A star is appreciated, too!)

---

### Step 1: Set up the environment

To create and activate the Conda environment, run:

```bash
conda env create --name VDDP_REPRODUCTION_ENV --file environment.yml
conda activate VDDP_REPRODUCTION_ENV
```

---

### Step 2: Run the reproduction scripts

After granting execution permissions using:

```bash
chmod +x ./scripts/*.sh
```

execute the following scripts **from the root directory of the repository** (the same directory containing this `README.md`):

```bash
./scripts/vddlm.sh
./scripts/vddlm-anal.sh
./scripts/vrr.sh
./scripts/vddgm.sh
```

If your server uses **SLURM** for job management, you may alternatively submit the scripts as batch jobs:

```bash
sbatch scripts/vddlm.sh
sbatch scripts/vddlm-anal.sh
sbatch scripts/vrr.sh
sbatch scripts/vddgm.sh
```

---

### Step 3: Generate the plots

To (re)generate all plots, run:

```bash
./scripts/plot.sh
```

---

### Step 4: Submit your reproduction report

Once you have completed your reproduction, please [open a new issue using the *Reproduction Report* template](https://github.com/jvhs0706/VDDP/issues/new?template=reproduction-report.md).
Proper attribution will be provided for all submitted reports.
