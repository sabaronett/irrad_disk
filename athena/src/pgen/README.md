# `irrad_disk.cpp`

## [`freq_table.txt`](https://github.com/sabaronett/irrad_disk/blob/main/athena/models/dsharp_abs-sca/nfreq/3/freq_table.txt)

If [`<problem>/frequency_table = 1`](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L89), this list of [$N_f$](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/models/dsharp_abs-sca/nfreq/3/athinput.dsharp#L66)$ - 1$ custom frequency band interfaces will be read by `irrad_disk.cpp` ([L258-L282](https://github.com/sabaronett/irrad_disk/blob/f7ee62d3f50dcabcd4281b8726955f0f91c207b6/athena/src/pgen/irrad_disk.cpp#L258-L282)).
All values must be negative or positivie which will be interpreted in code (i.e., $k_\mathrm{B}T_0/h$) or cgs (i.e., Hz) units, respectively.