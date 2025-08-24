# Opacities


## `athena_opac.py`

Along with a [RADMC-3D formatted](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#the-dustkappa-inp-files) monochromatic opacity file (e.g., [`dustkappa_dsharp.inp`](#dustkappa_dsharp.inp)), this script should be placed within an Athena++ model or run directory (e.g., [`models/dsharp_abs-sca/nfreq/64`](https://github.com/sabaronett/irrad_disk/tree/main/athena/models/dsharp_abs-sca/nfreq/64)).
Upon execution (e.g., via `$ python athena_opac.py`), it precomputes band-mean opacity coefficients based on the corresponding Athena++ input file (e.g., [`athinput.dsharp`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/64/athinput.dsharp)) and outputs the following files as required by our problem generator [`irrad_disk.cpp`](https://github.com/sabaronett/irrad_disk/blob/main/athena/src/pgen/irrad_disk.cpp):


### [`kappa_pf_table.txt`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/3/kappa_pf_table.txt)

[$N_f$](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L66) columns, with [`<problem>/n_temperature`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L86) rows, of $\kappa_{\mathrm{p},f\mathrm{(,e)}}$ [Baronett, Jiang et al. 2025, equation (8)].


### [`kappa_rf_table.txt`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/3/kappa_rf_table.txt)

[$N_f$](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L66) columns, with [`<problem>/n_temperature`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L86) rows, of $\kappa_{\mathrm{r},f}$ [Baronett, Jiang et al. 2025, equation (9)].


### [`kappa_sf_table.txt`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/3/kappa_sf_table.txt)

[$N_f$](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L66) columns, with [`<problem>/n_temperature`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L86) rows, of $\kappa_{\mathrm{s},f}$ [Baronett, Jiang et al. 2025, equation (10)].


### [`temp_table.txt`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/3/temp_table.txt)

[`<problem>/n_temperature`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L86) rows of temperatures in Kelvin, logarithmically spaced between [`<problem>/temperature_min`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L87) and [`<problem>/temperature_max`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L88).
[`irrad_disk.cpp`](https://github.com/sabaronett/irrad_disk/blob/main/athena/src/pgen/irrad_disk.cpp) linearly interpolates the opacities above based on the local thermal $T_\mathrm{therm}$ or color $T_\mathrm{c}$ temperature [Baronett, Jiang et al. 2025, Section 2.4].


## `dustkappa_dsharp.inp`

The tabulated monochromatic dust opacities we used in Baronett, Jiang et al. (2025, Section 3.2).
See file header for details.
