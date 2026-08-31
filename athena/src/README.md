# Athena++ Source Files

This directory follows the [Athena++ directory structure](https://github.com/PrincetonUniversity/athena/tree/main/src).

- [`nr_radiation/integrators`](nr_radiation/integrators/README.md) and [`task_list`](task_list/README.md) describe the modifications to files within those directories needed to prevent the disk from hydrodynamically evolving while allowing only the radiation field to evolve (hydrostatically) toward thermal equilibrium (Baronett, Jiang et al. ([in revision](https://ui.adsabs.harvard.edu/abs/2026arXiv260608859B/abstract)), Section 3.3).
- [`pgen`](pgen/README.md) contains our problem generator `irrad_disk.cpp` (which may also be included in a future public version of [Athena++](https://github.com/PrincetonUniversity/athena)), described in Baronett, Jiang et al. ([in revision](https://ui.adsabs.harvard.edu/abs/2026arXiv260608859B/abstract), Sections 2.2, 2.4, and 3.1)

Please see subdirectories for details.
