#!/usr/bin/env python3
#===============================================================================
# problem_setup.py
#
# Creates RADMC-3D input files for the protoplanetary disk model used in our
# analytic ray-tracing solutions. This script is adapted from
# `radmc3d-2.0/examples/run_ppdisk_simple_1/problem_setup.py`.
#
# Author: Stanley A. Baronett
# Created: 2024-12-16
# Updated: 2025-08-23
#===============================================================================
import athena_read
import numpy as np

# BEGIN `amr_grid.inp`==========================================================
# Read the Athena++ grid
root = '../athena/models/gray_opac/absorption'
run = 'multigroup'
problem_id = 'dsharp'
path = f'{root}/{run}'
athinput = athena_read.athinput(f'{path}/athinput.{problem_id}')
length_unit = athinput['radiation']['length_unit']
athdf = athena_read.athdf(f'{path}/athdf/{problem_id}.out1.00001.athdf')
nr = len(athdf['x1v'])
ntheta = len(athdf['x2v'])
nphi = len(athdf['x3v'])
ri = athdf['x1f']*length_unit
thetai = athdf['x2f']
phii = athdf['x3f']

# Write the grid file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-amr-grid-inp
with open('amr_grid.inp', 'w+') as f:
    f.write('1\n')                            # iformat
    f.write('0\n')                            # AMR grid style (0=regular grid)
    f.write('100\n')                          # Coordinate system: spherical
    f.write('0\n')                            # gridinfo
    f.write('1 1 0\n')                        # Include r,theta coordinates
    f.write(f'{nr:d} {ntheta:d} {nphi:d}\n')  # Size of grid
    for value in ri:
        f.write(f'{value:.16e}\n')            # X coordinates (cell walls)
    for value in thetai:
        f.write(f'{value:.16e}\n')            # Y coordinates (cell walls)
    for value in phii:
        f.write(f'{value:.16e}\n')            # Z coordinates (cell walls)
# END `amr_grid.inp`============================================================


# BEGIN `wavelength_micron.inp`=================================================
# Wavelength grid
lam1 = 0.1e0
lam2 = 7.0e0
lam3 = 25.e0
lam4 = 1.0e4
n12 = 20
n23 = 100
n34 = 30
lam12 = np.logspace(np.log10(lam1), np.log10(lam2), n12, endpoint=False)
lam23 = np.logspace(np.log10(lam2), np.log10(lam3), n23, endpoint=False)
lam34 = np.logspace(np.log10(lam3), np.log10(lam4), n34, endpoint=True)
lam = np.concatenate([lam12, lam23, lam34])
nlam = lam.size

# Write the wavelength file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-wavelength-micron-inp
with open('wavelength_micron.inp', 'w+') as f:
    f.write(f'{nlam:d}\n')
    for value in lam:
        f.write(f'{value:.16e}\n')
# END `wavelength_micron.inp`===================================================


# BEGIN `stars.inp`=============================================================
# Constants
ms = 1.98892e33                               # Solar mass [g]
ts = 5.78388e3                                # Solar temperature [K]
rs = 6.9368e10                                # Solar radius [cm]

# Star parameters
mstar = ms                                    # unimportant as of verion 2.0
rstar = rs                                    # for extended source treatment
tstar = ts                                    # blackbody temperature
pstar = np.array([0.,0.,0.])                  # position

# Write the stars.inp file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#sec-stars
with open('stars.inp', 'w+') as f:
    f.write('2\n')
    f.write(f'1 {nlam:d}\n\n')
    f.write(f'{rstar:.16e} {mstar:.16e} {pstar[0]:.16e} {pstar[1]:.16e} '\
            + f'{pstar[2]:.16e}\n\n')
    for value in lam:
        f.write(f'{value:.16e}\n')
    f.write(f'\n{-tstar:.16e}\n')
# END `stars.inp`===============================================================


# BEGIN `dust_density.inp`======================================================
# Make the cell-centered grid
density_unit = athinput['radiation']['density_unit']
rhod = athdf['rho']*density_unit

# Write the density file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-for-dust-transfer-dust-density-inp
with open('dust_density.inp', 'w+') as f:
    f.write('1\n')                            # Format number
    f.write(f'{nr*ntheta*nphi:d}\n')          # Nr of cells
    f.write('1\n')                            # Nr of dust species
    data = rhod.ravel()
    data.tofile(f, sep='\n', format='%.16e')
    f.write('\n')
# END `dust_density.inp`========================================================


# BEGIN `dustopac.inp`==========================================================
# Dust opacity control file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#the-dustopac-inp-file
with open('dustopac.inp', 'w+') as f:
    f.write('2               iformat (Format number of this file)\n')
    f.write('1               nspec (Nr of dust species)\n')
    f.write('============================================================================\n')
    f.write('1               inputstyle (Way in which this dust species is read)\n')
    f.write('0               iquantum (0=Thermal grain)\n')
    f.write('constant        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
# END `dustopac.inp`============================================================


# BEGIN `dustkappa_[name].inp`=================================================
name = 'constant'                             # dust species name (no spaces)
iformat = 1                                   # pure absorption (1)

kappa_a = athinput['problem']['kappa_a']      # absorption opacity [cm^2/g]

# Write the dust opacities file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#the-dustkappa-inp-files
with open(f'dustkappa_{name}.inp', 'w+') as f:
    f.write('# Opacity file for a constant dust opacity test\n')
    f.write('# Optical constants from...\n')
    f.write('# Please do not forget to cite in your publications the original paper of these optical constant measurements\n')
    f.write('# Made with the ... code by ...\n')
    f.write('# Grain size =  ? cm\n')
    f.write('# Material density =  ? g/cm^3\n')
    f.write(f'{iformat:d}\n')
    f.write(f'{nlam:d}\n\n')
    for value in lam:
        f.write(f'{value:.16e} {kappa_a:.16e}\n')
# END `dustkappa_[name].inp`===================================================


# BEGIN `radmc3d.inp`===========================================================
# Monte Carlo parameters
nphot = int(8e8)
nphot_mono = int(1e6)
countwrite = int(1e6)

# Write the radmc3d.inp control file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/clioptions.html#additional-arguments-general
with open('radmc3d.inp', 'w+') as f:
    f.write(f'countwrite = {countwrite}\n')   # nr photons between std. outputs
    f.write('iranfreqmode = 1\n')             # differ RNG seeds per OMP process
    f.write('istar_sphere = 0\n')             # point (0)/sphere (1) star source
    f.write(f'nphot = {nphot:d}\n')           # number of photons
    f.write(f'nphot_mono = {nphot_mono:d}\n') # number of monochromatic photons
    f.write('scattering_mode_max = 0\n')      # https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/dustradtrans.html?highlight=scattering_mode_max#five-modes-of-treating-scattering
    f.write('setthreads = 122\n')             # AMD Rome optimized
# END `radmc3d.inp`=============================================================


# BEGIN `mcmono_wavelength_micron.inp`=================================================
# Wavelength grid
lam1 = 0.1e0
lam2 = 7.0e0
lam3 = 25.e0
lam4 = 1.0e4
n12 = 20
n23 = 100
n34 = 30
lam12 = np.logspace(np.log10(lam1), np.log10(lam2), n12, endpoint=False)
lam23 = np.logspace(np.log10(lam2), np.log10(lam3), n23, endpoint=False)
lam34 = np.logspace(np.log10(lam3), np.log10(lam4), n34, endpoint=True)
lam = np.concatenate([lam12, lam23, lam34])
nlam = lam.size

# Write the wavelength file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/dustradtrans.html#special-purpose-feature-computing-the-local-radiation-field
with open('mcmono_wavelength_micron.inp', 'w+') as f:
    f.write(f'{nlam:d}\n')
    for value in lam:
        f.write(f'{value:.16e}\n')
# END `mcmono_wavelength_micron.inp`===================================================
