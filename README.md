<center><img src="AxionyxLogo.png" alt="Logo" width="560" align = "center"/></center>

# AxioNyx 1.1 Initial Release
April 2026

AxioNyx 1.1 is a fork of AxioNyx 1 maintained by Dr. Yourong Wang (FWPhys.com), building on the original codebase with the following enhancements:

- Restored functionality and improved stability
  - Joint Axion+Baryon+CDM (ABC) simulation support.
  - Updated API compatibility with AMReX Version 2025 and newer.

- Overall modernization and streamlining of IO and simulation loops.
  - Brand new problem initialiser with ability to read in a binary data file (from e.g. PyUltraLight).
  - Overhauled runlog writer that respects user-specified precision requirements.
  - Recalibrated pipeline to start Mixed FDM+CDM cosmology from MUSIC-generated files.

- Cosmetic and Quality of Life Improvements
  - New splash print with version information.
  - Improved output file headers for improved compatibility with ROCKSTAR.
  - Active migration from SWFFT to native AMReX FFT, removing the MPI rank count constraint present in version 1.0.

This initial release (Axionyx 1.1) accompanies Wang (2026), arXiv:2604.06038. The repository ships with one example configuration (`./Exec/MDM_NoHydro`) corresponding to the mixed FDM+CDM, dark-matter-only setup used in that paper.

# Getting Started

clone this repository to a folder of your choosing, and make sure `FFTW` is installed and findable.

```bash
cd Exec/MDM_NoHydro
make -j
```

To perform a simulation on your local machine,
```
mpirun -np 8 ./Nyx3D.xxx.EX inputs
```

# AxioNyx 1.1: Instructions on Cosmological Mixed-DM Simulations

## Cosmological Parameters

Cosmology enters an AxioNyx simulation at two levels: the initial conditions 
and the runtime inputs file. You should make sure they are consistent for each run. 

Key parameters in `inputs.uldm` include 
`nyx.comoving_h`, `nyx.comoving_OmM`, and `nyx.ratio_fdm` (the axion fraction 
of total matter, including baryons.

> **Unit note:** AxioNyx and Nyx natively work in comoving Mpc. The CAMB 
> family and most of the literature work in comoving Mpc/h. The box size in 
> `geometry.prob_hi` must reflect this conversion — a 5 Mpc/h box at h=0.675 
> should be entered as 7.407 comoving Mpc.

## Generating Initial Conditions

Mixed-DM initial conditions require separate displacement fields for the axion 
(FDM) and CDM components. The pipeline is:

AxionCAMB → transfer functions → Python converter → MUSIC → AxioNyx ICs

### Step 1: AxionCAMB

Generate transfer functions using [AxionCAMB](https://github.com/dgrin1/axionCAMB). 
An example parameter file is provided in `Exec/MDM_NoHydro/camb_inputs.ini`. 
Set `transfer_redshift(1)` to match your intended simulation start redshift.

### Step 2: Convert AxionCAMB output for MUSIC

AxionCAMB transfer function output is not directly readable by MUSIC. A Python 
conversion script is required. This repository ships a fork of the converter 
originally developed by [Alex Laguë](https://github.com/alexlague/AxionsIC).

A Jupyter Notebook based on this code is available in `./Exec/MDM_NoHydro/Scripts/axionCAMB_to_MUSIC.ipynb`. Please note that while `camb` has a Python wrapper and is used, `AxionCAMB` has to be built from source.

Follow the instructions in the notebook and convert your AxionCAMB transfer output before proceeding.

### Step 3: MUSIC

Download and build [MUSIC](https://bitbucket.org/ohahn/music/src/master/) with the Nyx (Boxlib) output plugin enabled. Generate ICs separately for the axion and CDM components. 

The example `inputs.uldm` expects two output directories:

```
cosmo.FDMinitDirName = /path/to/ICs/A/
cosmo.initDirName    = /path/to/ICs/C/
```

### Step 4: Inputs file

Set `nyx.initial_z` to match `transfer_redshift` in your CAMB parameter file. 
The axion mass `nyx.m_tt` is in units of 10⁻²² eV (so `m_tt = 0.1` corresponds 
to m = 10⁻²³ eV). Make sure this is consistent across the board as well.


# About AxioNyx 1.0
axionyx is a modification of the [Nyx] (https://github.com/AMReX-Astro/Nyx.git) code for dealing with Fuzzy Dark Matter. The AxioNyx paper https://arxiv.org/abs/2007.08256 was produced using the mixed_dm branch.

# License
axionyx is released under the LBL's modified BSD license, see the [license.txt](license.txt) file for details.
