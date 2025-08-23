# Athena++ Source Files

This directory follows the [Athena++ directory structure](https://github.com/PrincetonUniversity/athena/tree/main/src).

- [`nr_radiation/integrators`](https://github.com/sabaronett/irrad_disk/tree/main/athena/src/nr_radiation/integrators) and [`task_list`](https://github.com/sabaronett/irrad_disk/tree/main/athena/src/task_list) contain modifications needed to prevent the disk from dynamically evolving while allowing only the radiation field to evolve toward thermal equilibrium (Baronett, Jiang et al. 2025, Section 3.3).
- [`pgen`](https://github.com/sabaronett/irrad_disk/tree/main/athena/src/pgen) contains our problem generator `irrad_disk.cpp` (which may also be included in a future public version of [Athena++](https://github.com/PrincetonUniversity/athena)), described in Baronett, Jiang et al. (2025, Sections 2.2, 2.4, and 3.1)

Please see subdirectories for details.
