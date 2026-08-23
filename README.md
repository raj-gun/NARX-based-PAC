# NARX-PAC: A Dynamical Systems and System Identification Framework for Phase-Amplitude Coupling Analysis  
*An open-source MATLAB implementation for detecting and characterising phase-amplitude coupling using nonlinear system identification and polynomial NARX models.*

Authors: [Rajintha Gunawardena](https://www.researchgate.net/profile/Shenal-Gunawardena?ev=hdr_xprf)<sup>1</sup>, [Fei He](https://github.com/feihelab)<sup>1</sup>  
1. Centre for Computational Science and Mathematical Modelling, Coventry University, Coventry CV1 5FB, UK.

[![MATLAB](https://img.shields.io/badge/MATLAB-Code-blue)](https://www.mathworks.com/products/matlab.html)  
[![arXiv](https://img.shields.io/badge/arXiv-2603.08866-b31b1b.svg)](https://arxiv.org/abs/2603.08866)

---

## Overview 📖  
**NARX-PAC** is a MATLAB implementation of a dynamical-systems and nonlinear-system-identification framework for phase-amplitude coupling (PAC) analysis. The method represents PAC as quadratic phase coupling between a slow oscillation, a fast oscillation, and the two immediate intermodulation components. A two-input, single-output, second-order input-only NARX model is identified to obtain a canonical approximation of the underlying PAC dynamics.

Unlike methods based only on filtered phase and amplitude-envelope variations, NARX-PAC identifies a generative coupling model. Noise-free simulations of the identified model are then used to estimate modulation strength, preferred phase, and PAC type. The framework also includes procedures for rejecting harmonic-related spurious PAC and a discriminator $\mathcal{D}$ map for distinguishing genuine PAC from intermodulation-related spurious PAC.

In the experiments reported in the accompanying paper, the method produced sharper and more frequency-specific coupling localisation than the benchmark filtering-based methods. It remained robust at a signal-to-noise ratio (SNR) of 2, reasonably robust at an SNR of 1, and performed well with analysis windows as short as 3-5 seconds.

### Features  
- **Canonical NARX approximation**: Identifies the minimal dynamical structure associated with PAC using a second-order input-only polynomial NARX model.  
- **NARX-PAC modulation index**: Estimates modulation strength from the immediate-sideband magnitudes relative to the high-frequency component and supports interpretation of PAC type.  
- **Preferred-phase estimation**: Uses noise-free simulations of the identified model to determine the low-frequency phase at which the high-frequency amplitude reaches its maximum.  
- **Harmonic-related spurious PAC rejection**: Applies an instantaneous-frequency criterion derived from the identified canonical dynamics.  
- **Discriminator $\mathcal{D}$ map**: Supports the distinction between genuine PAC and intermodulation-related spurious PAC.  

---

## Getting Started 🚀  

### Prerequisites  
- MATLAB.  
- Required MATLAB Toolboxes:  
  - **Signal Processing Toolbox**.  
  - **Parallel Computing Toolbox**.

The standalone input-only system-identification code required by NARX-PAC is included in the [`NonSysID-i`](NonSysID-i/) folder.

### Installation  
1. Clone the repository:  
   ```bash
   git clone https://github.com/raj-gun/NARX-based-PAC.git
   ```
   or manually download the repository.

2. Update the paths at the start of the MATLAB example script that you want to run:
   ```matlab
   addpath('\<path-to>\NonSysID-i\');
   addpath('\<path-to>\NARX_PAC\');
   addpath('\<path-to>\NARX_PAC\Utils\');
   ```

### Basic use
The principal grid-search function is [`pac_miso_Cmdg_mod_21`](NARX_PAC/pac_miso_Cmdg_mod_21.m). A basic call has the following form:

```matlab
addpath('\<path-to>\NonSysID-i\');
addpath('\<path-to>\NARX_PAC\');
addpath('\<path-to>\NARX_PAC\Utils\');

fL_vals = 4:1:10;
fH_vals = 30:1:100;
filt_typ = {'sbp','sbp'}; 
frq_bndw_LF = 1;
frq_bndw_HF = 0.5;
RCT = 3;

[Comods, diff_comod, phs_data_mat, fL_grd, fH_grd, ...
    All_freq_comb, mod_trm_clstr_ERR, All_freq_comb_ARX_1, ...
    All_freq_comb_ARX_2, narx_pac_modls_1, narx_pac_modls_2] = ...
    pac_miso_Cmdg_mod_21(signal, fL_vals, fH_vals, Fs, RCT, ...
    filt_typ, frq_bndw_LF, frq_bndw_HF);

% Discriminator D map
D_map = diff_comod;
D_map(D_map > 0) = 1;
D_map(D_map < 0) = -1;

figure;
imagesc(fL_vals, fH_vals, Comods{1});
axis xy;
colorbar;
xlabel('Low frequency (Hz)');
ylabel('High frequency (Hz)');
```

`Comods{1}` contains the thresholded NARX-PAC comodulogram, `Comods{2}` contains the raw comodulogram, and `diff_comod` contains the normalised high-frequency-magnitude minus modulation-strength map used to construct the discriminator $\mathcal{D}$ map. The signed discriminator map is stored as `D_map` in the basic-use example. The example scripts show how to apply [`IF_harmonic_test`](NARX_PAC/IF_harmonic_test.m) and [`SpuCup_intrmd_2`](NARX_PAC/SpuCup_intrmd_2.m) for post-processing.

### Repository structure
- [`NARX_PAC`](NARX_PAC/) contains the principal NARX-PAC functions and utilities.  
- [`NonSysID-i`](NonSysID-i/) contains the standalone input-only system-identification routines used to identify a NARX model which estimates the canonical approximation of a phase-amplitude coupling.  
- [`NARX-PAC paper`](NARX-PAC%20paper/) contains the experiment scripts, saved data, and plotting files associated with the paper.

### Examples
- A synthetic 7 Hz and 63 Hz PAC example is provided in [`NARX-PAC paper/7_PAC_63/Fig 13`](NARX-PAC%20paper/7_PAC_63/Fig%2013/).  
- The non-stationary 6-7 Hz and 55-65 Hz example is provided in [`NARX-PAC paper/6-7_PAC_55-65`](NARX-PAC%20paper/6-7_PAC_55-65/).  
- The example with a 9-10 Hz slow oscillation coupled to two fast oscillations is provided in [`NARX-PAC paper/9-10_PAC_35-40_n_70-80`](NARX-PAC%20paper/9-10_PAC_35-40_n_70-80/).  
- Harmonic-, transient-, and spike-related spurious-PAC examples are provided in [`NARX-PAC paper/Spurious_PAC_harmonics_and_spikes`](NARX-PAC%20paper/Spurious_PAC_harmonics_and_spikes/).  
- The rat hippocampal LFP examples are provided in [`NARX-PAC paper/LFP`](NARX-PAC%20paper/LFP/).

Scripts beginning with `PAC_OthrMthds_` reproduce analyses using the benchmark PAC methods. These scripts require the corresponding external method implementations; update their `Methods\Matlab_Code` path before use.

## Paper

If you use this code for academic purposes, kindly reference our paper as follows:

**A Dynamical Systems and System Identification Framework for Phase-Amplitude Coupling Analysis**

Rajintha Gunawardena and Fei He

arXiv: [2603.08866](https://arxiv.org/abs/2603.08866)  
DOI: [10.48550/arXiv.2603.08866](https://doi.org/10.48550/arXiv.2603.08866)

```bibtex
@misc{Gunawardena2026,
  title = {A Dynamical Systems and System Identification Framework for Phase-Amplitude Coupling Analysis},
  author = {Gunawardena, Rajintha and He, Fei},
  year = {2026},
  eprint = {2603.08866},
  archivePrefix = {arXiv},
  primaryClass = {q-bio.NC},
  doi = {10.48550/arXiv.2603.08866},
  url = {https://arxiv.org/abs/2603.08866}
}
```
