# Optimization-Based-Stochastic-ILC
MATLAB code for "Monotonic Convergence of Stochastic Iterative Learning Control: An Optimization-Based Method"
This repository contains the official MATLAB implementation for the simulation results presented in the paper:

> **"Monotonic Convergence of Stochastic Iterative Learning Control: An Optimization-Based Method"** > *Wenjin Lv, Deyuan Meng, and Jingyao Zhang* > (Currently under review in *SIAM Journal on Control and Optimization (SICON)*)

## 📌 Overview

This code reproduces the comparative simulation (Section 5 of the manuscript), demonstrating the tracking and learning performance of the proposed optimization-based stochastic ILC algorithm. 

To provide a rigorous comparative analysis, the simulation includes:
1. **Proposed Algorithm 3.1**: An optimization-based stochastic ILC that strictly guarantees step-by-step monotonic variance reduction.
2. **SA-Shen [33]**: A stochastic approximation (SA)-based algorithm.
3. **SA-Cheng [10]**: An accelerated SA-based learning control scheme.

As analyzed in the manuscript, while the SA-based methods guarantee asymptotic convergence, they  suffer from non-monotonic transient growth due to the lack of dynamic variance optimization. The proposed algorithm perfectly resolves this inherent trade-off.

##  Requirements

* MATLAB (R2017b or newer is recommended).
* No external toolboxes are required. The code runs on basic MATLAB functions.

##  How to Run

1. Clone or download this repository to your local machine.
2. Open MATLAB and navigate to the downloaded folder.
3. Run the main script:
   ```matlab
   main_simulation.m
