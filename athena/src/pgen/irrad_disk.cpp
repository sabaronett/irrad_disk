//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file irrad_disk.cpp
//! \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//! spherical polar coordinates. Initial conditions are in vertical hydrostatic eqm.
//! The inner radial boundary admits stellar irradiation only along radial rays.
//! (Adapted from disk.cpp.)
//! 
//! Author: Stanley A. Baronett
//! Created: 2024-10-21
//! Updated: 2025-05-29
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // exp, log10, pow, sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

namespace {
// functions to construct and determine the disk profile and opacities
void GetCylCoord(Coordinates *pco, Real &rad, Real &phi, Real &z, int i, int j, int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
int BinarySearchIncreasing(AthenaArray<Real> &arr, int low, int high, const Real target);
Real LinearInterpolation(Real x, Real x0, Real x1, Real y0, Real y1);
void GetOpacities(const Real temp, const int ifr, Real &kappa_af, Real &kappa_pf);
void GetScatteringOpacities(const Real temp, const int ifr, Real &kappa_pf,
                            Real &kappa_rf, Real &kappa_sf);
void SetFrequencies(NRRadiation *prad);

// input file parameters which are useful to make global to this file
Real x1min;                                     // <mesh>
Real gamma_gas, dfloor;                         // <hydro>
Real Omega0;                                    // <orbital_advection>
Real t_unit, density_unit, length_unit;         // <radiation>
Real gm0, r0, rho0, p0_over_r0, pslope, dslope; // <problem> (disk)
Real r_star, t_star, kappa_a, kappa_s;          // <problem> (irradiation)

// for frequency dependent opacities
bool scattering;
int nfreq, ntemp, user_freq;
Real dlog10T;
AthenaArray<Real> freq_table;
AthenaArray<Real> temp_table;
AthenaArray<Real> kappa_pf_table;
AthenaArray<Real> kappa_rf_table;
AthenaArray<Real> kappa_sf_table;
} // namespace

// User-defined boundary conditions and radiation-related functions for disk simulations
void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void RadInnerX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void RadOuterX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem", "gm0", 1.0);
  r0 = pin->GetOrAddReal("problem", "r0", 1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem", "rho0");
  dslope = pin->GetOrAddReal("problem", "dslope", 0.0);

  // Get parameters of initial pressure and cooling
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem", "p0_over_r0", 0.0025);
    pslope = pin->GetOrAddReal("problem", "pslope", 0.0);
    gamma_gas = pin->GetReal("hydro", "gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro", "iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro", "dfloor", (1024*(float_min)));

  // Get parameters of orbital advection
  Omega0 = pin->GetOrAddReal("orbital_advection", "Omega0", 0.0);

  // Get stellar radiation parameters
  x1min = pin->GetOrAddReal("mesh", "x1min", 1.0);
  nfreq = pin->GetOrAddInteger("radiation", "n_frequency", 1);
  t_unit = pin->GetOrAddReal("radiation", "T_unit", 0.0);
  density_unit = pin->GetOrAddReal("radiation", "density_unit", 0.0);
  length_unit = pin->GetOrAddReal("radiation", "length_unit", 0.0);
  kappa_a = pin->GetOrAddReal("problem", "kappa_a", 0.0);
  kappa_s = pin->GetOrAddReal("problem", "kappa_s", 0.0);
  r_star = pin->GetOrAddReal("problem", "r_star", 0.001);
  t_star = pin->GetOrAddReal("problem", "t_star", 1.0);
  ntemp = pin->GetOrAddInteger("problem", "n_temperature", 0);
  user_freq = pin->GetOrAddInteger("problem", "frequency_table", 0);

  // Convert constant opacities as needed
  std::stringstream msg;
  if (kappa_a < 0) {
    kappa_a *= -1;
  } else if (density_unit == 0.0) {
    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
        << "`density_unit` input parameter needed to convert `kappa_a` to code units";
    ATHENA_ERROR(msg);

    return;
  } else if (length_unit == 0.0) {
    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
        << "`length_unit` input parameter needed to convert `kappa_a` to code units";
    ATHENA_ERROR(msg);

    return;
  } else {
    kappa_a *= density_unit*length_unit;
  }
  if (kappa_s < 0) {
    kappa_s *= -1;
  } else if (density_unit == 0.0) {
    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
        << "`density_unit` input parameter needed to convert `kappa_s` to code units";
    ATHENA_ERROR(msg);

    return;
  } else if (length_unit == 0.0) {
    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
        << "`length_unit` input parameter needed to convert `kappa_s` to code units";
    ATHENA_ERROR(msg);

    return;
  } else {
    kappa_s *= density_unit*length_unit;
  }

  // Prepare frequency- and temperature-dependent opacity tables
  if (nfreq > 1) {
    int ftemp_lines = 0;
    std::string line;
    FILE *ftemp_table = fopen("./temp_table.txt", "r");
    FILE *fkappa_pf_table = fopen("./kappa_pf_table.txt", "r");
    FILE *fkappa_rf_table = fopen("./kappa_rf_table.txt", "r");
    FILE *fkappa_sf_table = fopen("./kappa_sf_table.txt", "r");
    temp_table.NewAthenaArray(ntemp);
    kappa_pf_table.NewAthenaArray(ntemp, nfreq);
    kappa_rf_table.NewAthenaArray(ntemp, nfreq);
    kappa_sf_table.NewAthenaArray(ntemp, nfreq);

    // Open input table files
    if (ftemp_table == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
          << "Unable to open temp_table.txt for frequency-dependent opacities";
      ATHENA_ERROR(msg);

      return;
    }
    if (fkappa_pf_table == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
          << "Unable to open kappa_pf_table.txt for frequency-dependent opacities";
      ATHENA_ERROR(msg);

      return;
    }
    if (fkappa_rf_table == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
          << "Unable to open kappa_rf_table.txt for frequency-dependent opacities";
      ATHENA_ERROR(msg);

      return;
    }
    if (fkappa_sf_table == nullptr) {
      scattering = false;
      // std::cout << "### Scattering OFF" << std::endl;
    } else {
      scattering = true;
      // std::cout << "### Scattering ON" << std::endl;
    }

    // Check table sizes against input parameters
    std::ifstream temp_ifstream("temp_table.txt");
    while (std::getline(temp_ifstream, line))
      ++ftemp_lines;
    if (ftemp_lines != ntemp) {
      msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
          << "temp_table.txt size does not match `n_temperature` input parameter";
      ATHENA_ERROR(msg);

      return;
    } else {
      temp_ifstream.close();
    }
    if (user_freq == 1) {
      std::ifstream freq_ifstream("freq_table.txt");

      if (!freq_ifstream.is_open()) {
      msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
          << "Unable to open freq_table.txt for user-defined frequency groups";
      ATHENA_ERROR(msg);

      return;
      } else {
        int ffreq_lines = 0;

        while (std::getline(freq_ifstream, line))
          ++ffreq_lines;
        if (ffreq_lines != nfreq-1) {
          msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
              << "freq_table.txt size inconsistent with `n_frequency` input parameter";
          ATHENA_ERROR(msg);

          return;
        } else {
          freq_ifstream.close();
        }
      }
    }

    // Read file values into array tables
    for (int i=0; i<ntemp; ++i) {
      fscanf(ftemp_table, "%le", &(temp_table(i)));
      for (int j=0; j<nfreq; ++j) {
        fscanf(fkappa_pf_table, "%le", &(kappa_pf_table(i, j)));
        fscanf(fkappa_rf_table, "%le", &(kappa_rf_table(i, j)));
        if (scattering) { fscanf(fkappa_sf_table, "%le", &(kappa_sf_table(i, j))); }
      }
    }
    fclose(ftemp_table);
    fclose(fkappa_pf_table);
    fclose(fkappa_rf_table);
    fclose(fkappa_sf_table);
    dlog10T = std::log10(temp_table(1)) - std::log10(temp_table(0));
    if (user_freq == 1) {
      Real k_b = 1.3807e-16;  // Boltzmann constant [erg/K]
      Real h = 6.626196e-27;  // Planck constant [erg s]
      FILE *ffreq_table = fopen("./freq_table.txt", "r");
      freq_table.NewAthenaArray(nfreq-1);

      // std::cout << "### User-defined frequency groups in [Mesh::InitUserMeshData]"
      //           << std::endl;
      for (int i=0; i<nfreq-1; ++i) {
        fscanf(ffreq_table, "%le", &(freq_table(i)));
        // std::cout << "\tfreq_table(" << i << ") = " << std::scientific << freq_table(i)
        //           << " (raw input)" << std::endl;
        if (freq_table(i) < 0) {
          freq_table(i) *= -1;
        } else if (t_unit == 0.0) {
          msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]" << std::endl
              << "`T_unit` input parameter needed to convert `frequency_table` to code "
              << "units";
          ATHENA_ERROR(msg);

          return;
        } else {
          freq_table(i) /= k_b*t_unit/h;
        }
        // std::cout << "\tfreq_table(" << i << ") = " << std::scientific << freq_table(i)
        //           << " hnu/k_BT_0" << std::endl;
      }
      fclose(ffreq_table);
    }
  }

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);

    if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);

    if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);

  return;
}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//! used to initialize variables which are global to other functions in this file.
//! Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2+nfreq);

  // enroll user-defined opacity function
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    pnrrad->EnrollOpacityFunction(DiskOpacity);
    if (user_freq == 1)
      pnrrad->EnrollFrequencyFunction(SetFrequencies);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1, x2, x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        den = DenProfileCyl(rad,phi,z);
        vel = VelProfileCyl(rad,phi,z);
        if (porb->orbital_advection_defined)
          vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = den*vel;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0*den*vel;
        }

        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//! \brief Function called before generating output files
//========================================================================================
// void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
//   for(int k=ks; k<=ke; k++) {
//     for(int j=js; j<=je; j++) {
//       for(int i=is; i<=(ie+NGHOST); i++) {
//         for (int n=0; n<pnrrad->nang; ++n) {
//           user_out_var(n,k,j,i) = pnrrad->ir(k,j,i-NGHOST,n); // store intensities
//         }
//       }
//     }
//   }
// }

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }

  return;
}

//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates
Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  Real EULER = 2.7182818284590452;

  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow((rad + r0)/r0, dslope)\
                /(1 + std::exp(-std::exp(EULER)*(rad - r0)/r0));
  Real dentem = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;

  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates
Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);

  return poverr;
}

//----------------------------------------------------------------------------------------
//! computes rotational velocity in cylindrical coordinates
Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;

  return vel;
}

//----------------------------------------------------------------------------------------
//! Iteratively use binary search on a strictly increasing array of Reals to find the
//! index that right-brackets the target, i.e., array[index - 1] < target < array[index]
//! NOTE: CONSIDER REMOVING THIS CONVENIENCE FUNCTION
int BinarySearchIncreasing(AthenaArray<Real> &arr, int low, int high, const Real target) {
  int mid;
  std::stringstream msg;

  while (low <= high) {
    mid = low + (high - low)/2;
    if ((arr(mid-1) < target) && (target < arr(mid))) {
      return mid;
    } else if (arr(mid) < target) {
      low = mid;
    } else
      high = mid;
  }

  msg << "### FATAL ERROR in function [BinarySearchIncreasing]" << std::endl
      << "Array may not be strictly increasing";
  ATHENA_ERROR(msg);
  return -1;
}

//----------------------------------------------------------------------------------------
//! Linear interpolation between two given points
//! TODO: Check if src/utils/interp_table.cpp can be used instead of this
//!       (will need to #include "../utils/interp_table.hpp")
Real LinearInterpolation(Real x, Real x0, Real x1, Real y0, Real y1) {
  return y0 + (x - x0)*(y1 - y0)/(x1 - x0);
}

//----------------------------------------------------------------------------------------
//! Gets or interpolates the Rosseland and Planck mean absorption opacities for the given
//! frequency bin and the given temperature.
void GetOpacities(const Real temp, const int ifr, Real &kappa_af, Real &kappa_pf) {
  if (temp < temp_table(0)) {
    kappa_af = kappa_rf_table(0, ifr);
    kappa_pf = kappa_pf_table(0, ifr);
  } else if (temp > temp_table(ntemp-1)) {
    kappa_af = kappa_rf_table(ntemp-1, ifr);
    kappa_pf = kappa_pf_table(ntemp-1, ifr);
  } else {
    // For a non-log-uniform temp_table, use the following instead:
    // int i = BinarySearchIncreasing(temp_table, 0, ntemp-1, temp);
    int i = int((std::log10(temp) - std::log10(temp_table(0)))/dlog10T) + 1;
    kappa_af = LinearInterpolation(temp, temp_table(i-1), temp_table(i),
                                   kappa_rf_table(i-1, ifr), kappa_rf_table(i, ifr));
    kappa_pf = LinearInterpolation(temp, temp_table(i-1), temp_table(i),
                                   kappa_pf_table(i-1, ifr), kappa_pf_table(i, ifr));
  }

  return;
}

//----------------------------------------------------------------------------------------
//! Gets or interpolates the scattering and total Rosseland and Planck-absorption mean
//! opacities for the given frequency bin and the given temperature.
void GetScatteringOpacities(const Real temp, const int ifr, Real &kappa_pf,
                            Real &kappa_rf, Real &kappa_sf) {
  if (temp < temp_table(0)) {
    kappa_pf = kappa_pf_table(0, ifr);
    kappa_rf = kappa_rf_table(0, ifr);
    kappa_sf = kappa_sf_table(0, ifr);
  } else if (temp > temp_table(ntemp-1)) {
    kappa_pf = kappa_pf_table(ntemp-1, ifr);
    kappa_rf = kappa_rf_table(ntemp-1, ifr);
    kappa_sf = kappa_sf_table(ntemp-1, ifr);
  } else {
    int i = BinarySearchIncreasing(temp_table, 0, ntemp-1, temp);
    kappa_pf = LinearInterpolation(temp, temp_table(i-1), temp_table(i),
                                   kappa_pf_table(i-1, ifr), kappa_pf_table(i, ifr));
    kappa_rf = LinearInterpolation(temp, temp_table(i-1), temp_table(i),
                                   kappa_rf_table(i-1, ifr), kappa_rf_table(i, ifr));
    kappa_sf = LinearInterpolation(temp, temp_table(i-1), temp_table(i),
                                   kappa_sf_table(i-1, ifr), kappa_sf_table(i, ifr));
  }

  return;
}

//----------------------------------------------------------------------------------------
//! Sets the frequency grid to the provided frequency table
void SetFrequencies(NRRadiation *prad) {
  prad->nu_grid(0) = 0.0;

  for(int i=0; i<nfreq-1; ++i) {
    prad->nu_grid(i+1) = freq_table(i);
    prad->nu_cen(i) = (prad->nu_grid(i) + prad->nu_grid(i+1))/2;
    prad->delta_nu(i) = prad->nu_grid(i+1) - prad->nu_grid(i);
  }
}

//========================================================================================
//! \fn int GetMaxErfDnu(const AthenaArray<Real> erf_dnu)
//! \brief Returns the index of the finite frequency band with the maximum specific
//! radiation energy density.
//========================================================================================
int GetMaxErfDnu(const AthenaArray<Real> erf_dnu) {
  int f = 0;

  for(int ifr=1; ifr<nfreq-1; ++ifr) {
    if (erf_dnu(ifr) > erf_dnu(f))
      f = ifr;
  }

  return f;
}

//========================================================================================
//! \fn Real GetColorTemp(const Real nu_peak)
//! \brief Returns the color temperature for the peak of the spectrum, using the Wien
//! displacement law.
//========================================================================================
Real GetColorTemp(const Real nu_peak) {
  return nu_peak/2.82;
}

//========================================================================================
//! \fn Real GetNuPeak(const Real temp)
//! \brief Returns the frequency of the peak of the Planck law B_nu(temp), using the Wien
//! displacement law.
//========================================================================================
Real GetNuPeak(const Real temp) {
  return 2.82*temp;
}

//========================================================================================
//! \fn int GetFreqGroup(const NRRadiation *prad, const Real nu)
//! \brief Returns the frequency band that contains the given frequency `nu`.
//========================================================================================
int GetFreqGroup(const NRRadiation *prad, const Real nu) {
  int f;
  
  if (nu < prad->nu_grid(1)) {
    f = 0;
  } else if (nu > prad->nu_grid(nfreq-1)) {
    f = nfreq-1;
  } else {
    for (int i=1; i<nfreq-1; ++i) {  // USE BinarySearchIncreasing
      if (nu > prad->nu_grid(i) && nu < prad->nu_grid(i+1)) {
        f = i;
        break;
      }
    }
  }

  return f;
}

} // namespace

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real flux;
  int &nang = prad->nang;                                               // total angles

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        flux = std::pow(t_star, 4)*std::pow(r_star/pco->x1v(is-i), 2)/4;
        if (nfreq == 1) {                                               // gray approx
          for (int n=0; n<nang-2; ++n) {                                // non-radial rays
            if (prad->mu(0,k,j,is-i,n) < 0.0)                           // exiting rays
              ir(k,j,is-i,n) = ir(k,j,is,n)\
                               *std::pow(pco->x1v(is)/pco->x1v(is-i), 2);
            else                                                        // entering rays
              ir(k,j,is-i,n) = 0.0;                                     // disable
          }                                                             // radial rays
          ir(k,j,is-i,nang-2) = flux/prad->wmu(nang-2);                 // n-hat_r,out
          ir(k,j,is-i,nang-1) = ir(k,j,is,nang-1)\
                                *std::pow(pco->x1v(is)/pco->x1v(is-i), 2);
        } else {                                                        // multifrequency
          for (int ifr=0; ifr<nfreq; ++ifr) {                           // each band
            for (int n=0; n<nang-2; ++n) {                              // non-radial rays
              if (prad->mu(0,k,j,is-i,ifr*nang+n) < 0.0)                // exiting rays
                ir(k,j,is-i,ifr*nang+n) = ir(k,j,is,ifr*nang+n)\
                                          *std::pow(pco->x1v(is)/pco->x1v(is-i), 2);
              else                                                      // entering rays
                ir(k,j,is-i,ifr*nang+n) = 0.0;                          // disable
            }
            if (ifr < nfreq-1) {                                        // n-hat_r,out
              ir(k,j,is-i,ifr*nang+nang-2) = flux*prad->IntPlanckFunc(prad->nu_grid(ifr)\
                                             /t_star, prad->nu_grid(ifr+1)/t_star)\
                                             /prad->wmu(nang-2);
            } else {                                                    // n-hat_r,out
              ir(k,j,is-i,ifr*nang+nang-2) = flux*(1.0\
                                             - prad->FitIntPlanckFunc(prad->nu_grid(ifr)\
                                             /t_star))/prad->wmu(nang-2);
            }                                                           // n-hat_r,in
            ir(k,j,is-i,ifr*nang+nang-1) = ir(k,j,is,ifr*nang+nang-1)\
                                           *std::pow(pco->x1v(is)/pco->x1v(is-i), 2);
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int &nang = prad->nang;                                         // total angles

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        if (nfreq == 1) {                                         // gray approx
          for (int n=0; n<nang; ++n) {                            // non-radial rays
            if (prad->mu(0,k,j,ie+i,n) > 0.0)                     // exiting rays
              ir(k,j,ie+i,n) = ir(k,j,ie,n)\
                               *std::pow(pco->x1v(ie)/pco->x1v(ie+i), 2);
            else                                                  // entering rays
              ir(k,j,ie+i,n) = 0.0;                               // disable
          }                                                       // radial rays
        } else {                                                  // multifrequency
          for (int ifr=0; ifr<nfreq; ++ifr) {                     // each band
            for (int n=0; n<nang; ++n) {                          // non-radial rays
              if (prad->mu(0,k,j,ie+i,ifr*nang+n) > 0.0)          // exiting rays
                ir(k,j,ie+i,ifr*nang+n) = ir(k,j,ie,ifr*nang+n)\
                                          *std::pow(pco->x1v(ie)/pco->x1v(ie+i), 2);
              else                                                // entering rays
                ir(k,j,ie+i,ifr*nang+n) = 0.0;                    // disable
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadInnerX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang = prad->nang;       // total n-hat angles N

  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<nang; ++n) {
          ir(k,js-j,i,n) = 0.0;  // should only be for entering rays
          //! TODO: use 1/r^2 decay for exiting rays (disk emission)
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadOuterX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang = prad->nang;       // total n-hat angles N

  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<nang; ++n) {
          ir(k,je+j,i,n) = 0.0;  // should only be for entering rays
          //! TODO: use 1/r^2 decay for exiting rays (disk emission)
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang = prad->nang;       // total n-hat angles N

  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<nang; ++n) {
          ir(ks-k,j,i,n) = 0.0;  // should only be for entering rays
          //! TODO: use 1/r^2 decay for exiting rays (disk emission)
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang = prad->nang;       // total n-hat angles N

  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<nang; ++n) {
          ir(ke+k,j,i,n) = 0.0;  // should only be for entering rays
          //! TODO: use 1/r^2 decay for exiting rays (disk emission)
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = vel;
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = vel;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = vel;
          prim(IM3,k,jl-j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = 0.0;
          prim(IM3,k,jl-j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = vel;
          prim(IM3,k,ju+j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = 0.0;
          prim(IM3,k,ju+j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = vel;
          prim(IM3,kl-k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = 0.0;
          prim(IM3,kl-k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = vel;
          prim(IM3,ku+k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = 0.0;
          prim(IM3,ku+k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! Sets opacities throughout the domain
void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim) {
  NRRadiation *prad = pmb->pnrrad;
  AthenaArray<Real> erf_dnu;  // specific radiation energy density per frequency band
  int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;
  int &nang = prad->nang;     // total angles
  int f_irrad;                // frequency band with peak Er/\Delta\nu
  int f_therm;                // frequency band with peak thermal emission
  Real t_therm;               // thermal temperature
  Real t_c;                   // color temperature
  Real kappa_sf = kappa_s;    // Rosseland mean scattering opacity (gray by default)
  Real kappa_af = kappa_a;    // Rosseland mean absorption opacity (gray by default)
  Real kappa_pfe = kappa_a;   // J_0 (Plank mean) absorp. coef.    (gray by default)
  Real kappa_pf = kappa_a;    // epsilon_0 (Plank mean) abs. coef. (gray by default)
  Real kappa_rf;              // Rosseland mean total (abs+sca) opacity
  Real dummy;                 // dummy variable for GetOpacities

  // Include ghost zones in upper/lower directional limits
  // il -= NGHOST;
  // iu += NGHOST;
  // if(ju > jl){
  //   jl -= NGHOST;
  //   ju += NGHOST;
  // }
  // if(ku > kl){
  //   kl -= NGHOST;
  //   ku += NGHOST;
  // }

  if (nfreq > 1)
    erf_dnu.NewAthenaArray(nfreq-1);

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        if (nfreq > 1) {
          t_therm = prim(IEN,k,j,i)/prim(IDN,k,j,i);
          f_therm = GetFreqGroup(prad, GetNuPeak(t_therm));
          for (int ifr=0; ifr<nfreq-1; ++ifr) {
            erf_dnu(ifr) = 0.0;                                  // reset from last cell
            for (int n=0; n<nang; ++n) {
              erf_dnu(ifr) += prad->wmu(n)*prad->ir(k,j,i,ifr*nang+n);
            }
            erf_dnu(ifr) /= prad->delta_nu(ifr);
          }
          f_irrad = GetMaxErfDnu(erf_dnu);
          t_c = GetColorTemp(prad->nu_cen(f_irrad));
          pmb->user_out_var(0,k,j,i) = t_c;
          pmb->user_out_var(1,k,j,i) = 0.0;
          if (f_irrad != f_therm) { pmb->user_out_var(1,k,j,i) = 1.0; }
        }
        for (int ifr=0; ifr<nfreq; ++ifr) {
          if (nfreq > 1) {
            if (scattering) {
              GetScatteringOpacities(t_therm, ifr, kappa_pf, kappa_rf, kappa_sf);
              kappa_af = kappa_rf - kappa_sf;
            } else {
              GetOpacities(t_therm, ifr, kappa_af, kappa_pf);
            }
            if (f_irrad != f_therm) {                            // non-irrad.-therm. eq.
              GetOpacities(t_c, ifr, dummy, kappa_pfe);
            } else {                                             // irradiative-thermal eq.
              GetOpacities(t_therm, ifr, dummy, kappa_pfe);
            }
            pmb->user_out_var(2+ifr,k,j,i) = kappa_pfe;
          }
          prad->sigma_s(k,j,i,ifr) = prim(IDN,k,j,i)*kappa_sf;
          prad->sigma_a(k,j,i,ifr) = prim(IDN,k,j,i)*kappa_af;
          prad->sigma_pe(k,j,i,ifr) = prim(IDN,k,j,i)*kappa_pfe; // J_0 coefficient
          if (prad->set_source_flag == 0) {                      // Pure attenuation case
            prad->sigma_p(k,j,i,ifr) = 0;                        // \epsilon_0 coefficient
          } else {
            prad->sigma_p(k,j,i,ifr) = prim(IDN,k,j,i)*kappa_pf;
          }
        }
      }
    }
  }
}
