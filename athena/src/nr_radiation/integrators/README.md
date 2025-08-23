# `rad_source.cpp`

For consistency and to improve the performance of the iterative implicit solver ([Jiang 2021](https://ui.adsabs.harvard.edu/abs/2021ApJS..253...49J/abstract), $\S$ 3.2.2), lines similar to the following must be commented out to remove radiation pressure as a source of momentum for the fluid:

```cpp
...
void RadIntegrator::AddSourceTerms(MeshBlock *pmb, AthenaArray<Real> &u) {
  ...
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // u(IM1,k,j,i) += rad_source(1,k,j,i);
        // u(IM2,k,j,i) += rad_source(2,k,j,i);
        // u(IM3,k,j,i) += rad_source(3,k,j,i);
        ...
```
