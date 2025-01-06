# VDDP: Verifiable Distributed Differential Privacy

Welcome to the official implementation of VDDP by anonymous authors.

## Prerequisites

To set up the conda environment named `VDP_ENV`, run:

```bash
conda create --name VDP_ENV --file req.txt
```

Note: The environment is named "VDP" instead of "VDDP" to accommodate generic verifiable DP mechanisms, not just VDDP.

## Reproduction

To reproduce the results of VDDLM (Figure 5 and Table 2), execute `scripts/vddlm-anal.sh` and `scripts/vddlm.sh`. We use [VDBM (BC23)](https://github.com/abiswas3/Verifiable-Differential-Privacy) as the baseline.

To reproduce the results of VRR (Table 3 and Figure 7), execute `scripts/vrr.sh`. We use [KCY21]() as the baseline.