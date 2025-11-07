# Pattern-Coupled Sparse Bayesian Learning with Fixed-Point Iterations – MATLAB Implementation

This repository contains MATLAB implementations supporting the paper:

**Didem Doğan and Geert Leus**,  
*“Pattern Coupled Sparse Bayesian Learning with Fixed Point Iterations for DOA and Amplitude Estimation,”*  
in *Proc. 2023 57th Asilomar Conference on Signals, Systems, and Computers*, IEEE, 2023.  
[https://ieeexplore.ieee.org/abstract/document/10477000](https://ieeexplore.ieee.org/abstract/document/10477000)

The code reproduces the numerical experiments reported in the paper, which introduces a **fast fixed-point variant of Pattern-Coupled Sparse Bayesian Learning (FP-PCSBL)** for **block-sparse signal recovery** and **DOA estimation**.

---

## 📘 Overview

This project implements **Pattern-Coupled Sparse Bayesian Learning (PCSBL)** and its **proposed fast fixed-point version (FP-PCSBL)**, comparing them to standard **SBL**-based algorithms.

Five algorithms are benchmarked:
- **EM-SBL** – Expectation–Maximization Sparse Bayesian Learning  
- **EM-PCSBL** – Pattern-Coupled SBL with EM updates  
- **EM-CSBL** – Correlated SBL (baseline variant)  
- **FP-SBL** – Fixed-Point SBL (non-coupled)  
- **FP-PCSBL** – *Proposed Algorithm*: Pattern-Coupled SBL using **fixed-point iterations**

The **FP-PCSBL** method achieves similar accuracy to EM-based PCSBL but with **significantly faster convergence**, eliminating the need for hyperparameter approximations.

---

## 🧩 Repository Structure

| Category | Files | Description |
|-----------|-------|-------------|
| **Main experiments** | `NMSE_Asilomar.m`, `NMSE_Asilomar_snapshot.m`, `NMSE_asilomarK.m`, `NMSE_asilomar_SNR.m` | Generate NMSE plots under different experimental conditions (m/n, snapshots, K, SNR). |
| **Core solvers** | `SBL.m`, `SBL_PC.m`, `MPCSBL.m`, `MPCSBL_alternative.m` | Implement FP SBL, FP PCSBL, EM SBL, EM PCSBL, and EM CSBL
| **Helpers & utilities** |`SBLSet.m`| Support functions for structure setup and hyperparameter updates. |

---

## 🧠 Methodology Summary

Each simulation follows these steps:
1. **Generate block-sparse signals** with random block locations and amplitudes.  
2. **Form the DOA steering matrix** for a uniform linear array (ULA).  
3. **Add complex Gaussian noise** at varying SNR levels.  
4. **Recover signals** using the five algorithms (EM and FP variants).  
5. **Compute NMSE** between the estimated and true magnitudes, averaged over Monte Carlo trials.

---

## 📊 Experiments

| Script | Variable Swept | Plot |
|---------|----------------|------|
| `NMSE_Asilomar.m` | **N/L ratio** | NMSE vs. number of sensors |
| `NMSE_Asilomar_snapshot.m` | **Snapshots (M)** | NMSE vs. number of snapshots |
| `NMSE_asilomarK.m` | **K (sparsity level)** | NMSE vs. K |
| `NMSE_asilomar_SNR.m` | **SNR (dB)** | NMSE vs. SNR |

Each figure shows five curves: **EM-SBL**, **EM-PCSBL**, **EM-CSBL**, **FP-SBL**, and **FP-PCSBL (proposed)**.

---

## ⚙️ Default Parameters

- **Speed of sound:** 1500 m/s  
- **Frequency:** 200 Hz → **λ = c / f**  
- **Sensor spacing:** λ / 2  
- **Angular grid:** −90° to 90°, 0.5° step  
- **Snapshots:** 1–50  
- **SNR range:** −10 dB to 30 dB  
- **Blocks:** L = 3  
- **Monte Carlo trials:** typically 100  

---

## 🧾 Reference

> D. Doğan and G. Leus,  
> *“Pattern Coupled Sparse Bayesian Learning with Fixed Point Iterations for DOA and Amplitude Estimation,”*  
> in *Proc. 2023 57th Asilomar Conference on Signals, Systems, and Computers*, IEEE, 2023.  
> DOI: [10.1109/ASILOMAR59739.2023.10477000](https://ieeexplore.ieee.org/abstract/document/10477000)

---

## 👤 Author

**Didem Doğan Başkaya**  
- LinkedIn: [in/didemdoganbaskaya](https://www.linkedin.com/in/didemdoganbaskaya)  
- GitHub: [Didemld](https://github.com/Didemld)
