A collection of scripts used for the paper 'Charm baryons at a variety of temperatures'

A key idea behind all these scripts is using a '[toml](https://github.com/toml-lang/toml)' file as input. This file specifies aims to specify any changeable features of the analysis as well as what to do the analysis on.

Example input files (for example only, they reference data which is not included herein) are in the 'exampleToml' folder. Example output files for most scripts are also in the 'exampleToml/output' folder.

The scripts themselves are in the 'scripts' folder

Scripts are run using i.e. python scriptname.py tomlfile.toml

---

# Figure Details
* Figure 1
  * Plots the correlators
  * plotG.py
* Figures 2-4
  * Does the parity-doubling R-ratio analysis
  * gvarParityRatio.py
  * plotInflectionPoints.py
* Figures 5-8
  * Does fits to the correlators and model averages
  * simpleFit.py
* Figure 9
  * Plots the spectrum of masses against experiment
  * plotSpectrum.py
* Figures 10-13
  * Plots the spectrum of masses as a function of temperature
  * simplePlotEMP.py
* Figure 14
  * Considers the Gell-Mann-Okubo relation
  * GellMannOkubo.py

---

# Conda Notes

## Install Environment
```conda env create -f environment.yml```
## Activate/Use
```conda activate charm```
## Update (w. new packages)

1. Edit `environment.yml`
2. Deactivate conda environment with `conda deactivate`
3. Update conda environment with `conda env update -f=environment.yml`
