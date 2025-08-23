# Opacities

separate Python scripts to precompute band-mean opacity coefficients based on an Athena++ input file (Section~\ref{sec:emissivity_opacity}) and to create and format the required input files for each RADMC-3D model presented in Section~\ref{sec:hydrostatic_comparisons};

## `athena_opac.py`

Along with a [RADMC-3D formatted]() monochromatic opacity file (e.g., [`dustkappa_dsharp.inp`](#dustkappa_dsharp.inp)), this script should be placed within an Athena++ model or run directory (e.g., [`models/dsharp_abs-sca/nfreq/64`](https://github.com/sabaronett/irrad_disk/tree/main/athena/models/dsharp_abs-sca/nfreq/64)).
Upon execution (e.g., via `$ python athena_opac.py`), it precomputes band-mean opacity coefficients based on the corresponding Athena++ input file (e.g., [`athinput.dsharp`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/64/athinput.dsharp)) and outputs the following files as required by our problem generator [`irrad_disk.cpp`](https://github.com/sabaronett/irrad_disk/blob/main/athena/src/pgen/irrad_disk.cpp):


### `freq_table.txt`
asdf

### `kappa_pf_table.txt`
asdf

### `kappa_rf_table.txt`
asdf

### `kappa_sf_table.txt`
asdf

### `temp_table.txt`
asdf


## `dustkappa_dsharp.inp`

The tabulated monochromatic dust opacities we used in Baronett, Jiang et al. (2025, Section 3.2).
See file header for details.
