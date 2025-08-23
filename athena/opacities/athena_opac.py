#!/usr/bin/env python3
#==============================================================================
# athena_opac.py
#
# Creates multifrequency opacity tables for input into Athena++.
#
# This script takes as input (1) a RADMC-3D-formatted opacity table with the
# filename structure `dustkappa_*.inp` (see https://www.ita.uni-heidelberg.de/
# ~dullemond/software/radmc-3d/manual_radmc3d/
# inputoutputfiles.html#the-dustkappa-inp-files); (2) an Athena++-formatted
# input file with the following requred parameters:
#   <radiation>
#   nfreq            # no. of frequency groups
#   frequency_min    # [0, \nu_min) [k_BT_0/h] < 0 < [Hz]
#   frequency_max    # [\nu_max, \infty) [k_BT_0/h] < 0 < [Hz] (for nfreq > 2)
#
#   <problem>
#   n_temperature    # no. of temperature groups
#   temperature_min  # min mean opacity temperature [K]
#   temperature_max  # max mean opacity temperature [K]
#
# Author: Stanley A. Baronett
# Created: 2024-04-19
# Updated: 2025-03-01
#===============================================================================
import numpy as np
from pathlib import Path
from radmc3dPy import analyze
from radmc3dPy.natconst import *
from scipy import integrate
from scipy.constants import c, h, k
import sys
sys.path.insert(0,'/home/stanley/github/PrincetonUniversity/athena/vis/python')
import athena_read
import warnings

# Convert constants from SI to cgs
c *= 1e2
h *= 1e7
k *= 1e7

def GetBnu_table(Ts, nus):
    """Computes Planck's law for a table of temperatures and frequencies
    """
    table = np.zeros((len(Ts), len(nus)))

    for i, T in enumerate(Ts):
        for j, nu in enumerate(nus):
            prefactor = 2*(k*T)**3/(h*c)**2
            u = h*nu/k/T
            if u < 0.001:  # Rayleigh--Jeans Law
                table[i][j] = prefactor*u**2
            elif u > 15:   # Wien Law
                table[i][j] = prefactor*u**3*np.exp(-u)
            else:          # Planck Law
                table[i][j] = prefactor*u**3/(np.exp(u) - 1)

    return table

def GetdBnu_dT_table(Ts, nus, diag=False):
    """Computes the partial derivative of Planck's law with respect to
    temperature
    """
    table = np.zeros((len(Ts), len(nus)))

    for i, T in enumerate(Ts):
        for j, nu in enumerate(nus):
            prefactor = 2*k**3*(T/h/c)**2
            u = h*nu/k/T
            if u < 0.001:
                table[i][j] = prefactor*(u**2 - u**4)
            elif u > 15:
                table[i][j] = prefactor*u**4*np.exp(-u)
            else:
                table[i][j] = prefactor*u**4*np.exp(u)/(np.exp(u) - 1)**2

    return table

def BinarySearchIncreasing(arr, low, high, target):
    """Iterative binary search on a strictly increasing array.

    Iteratively use binary search on a strictly increasing array to find the
    index that right-brackets the target, arr[mid-1] < target < arr[mid]
    """    
    while (low <= high):
        mid = int(low + (high - low)//2)
        if ((arr[mid-1] < target) and (target < arr[mid])):
            return mid
        elif (arr[mid] < target):
            low = mid
        else:
            high = mid

    raise Exception("Array may not be strictly increasing")

def PlanckMeanOpacities(kappa_nu, Bnu, nu, temp_table):
    numer = integrate.simpson(kappa_nu*Bnu, x=nu)
    denom = integrate.simpson(Bnu, x=nu)
    kappa = numer/denom
    # TODO: supress warnings for invalid values
    # (e.g., division by zero, invalid value encountered in double_scalars)
    np.nan_to_num(kappa, copy=False, nan=kappa_nu[0])
    return kappa

def RosselandMeanOpacities(kappa_nu, dBnu_dT, nu):
    numer = integrate.simpson(dBnu_dT/kappa_nu, x=nu)
    denom = integrate.simpson(dBnu_dT, x=nu)
    kappa = denom/numer
    # TODO: supress warnings for invalid values
    # (e.g., division by zero, invalid value encountered in double_scalars)
    np.nan_to_num(kappa, copy=False, nan=kappa_nu[0], posinf=kappa_nu[0])
    return kappa

# Read absorption coefficient as a function of frequency
fname = list(Path('./').glob(f'dustkappa_*.inp'))[0].parts[0]
ext = fname[10:-4]
opac = analyze.readOpac(ext=['dsharp'])
opac_freq = np.flip(1e4*c/opac.wav[0])
opac_kabs = np.flip(opac.kabs[0])
opac_ksca = np.flip(opac.ksca[0])
scattering = False
if len(opac_ksca) == len(opac_kabs):
    scattering = True
    opac_comb = opac_kabs + opac_ksca

# Make tables to compute and save mean opacities
fname = list(Path('./').glob(f'athinput.*'))[0].parts[0]
athinput = athena_read.athinput(fname)
T_unit = athinput['radiation']['T_unit']                          # [K]
density_unit = athinput['radiation']['density_unit']              # [g/cm^3]
length_unit = athinput['radiation']['length_unit']                # [cm]
ntemp = athinput['problem']['n_temperature']
temperature_min = athinput['problem']['temperature_min']          # [K]
temperature_max = athinput['problem']['temperature_max']          # [K]
temp_table = np.logspace(np.log10(temperature_min), np.log10(temperature_max),
                         ntemp)
Bnu_table = GetBnu_table(temp_table, opac_freq)
dBnu_dT_table = GetdBnu_dT_table(temp_table, opac_freq)

# For a single or multiple frequency bands
nu_min, nu_max = opac_freq[0], opac_freq[-1]
nfreq = athinput['radiation']['n_frequency']
kappa_pf_table = np.zeros((ntemp, nfreq))
kappa_rf_table = np.zeros((ntemp, nfreq))
kappa_sf_table = np.zeros((ntemp, nfreq))
i_nu0 = 0                   # left frequency table index
i_nu1 = len(opac_freq) - 1  # right index
nu_grid = np.asarray([nu_min, nu_max])
if nfreq > 1:
    frequency_min = athinput['radiation']['frequency_min']        # [Hz]
    if frequency_min < 0:  # unit switch: code (<0) or cgs (>0)
        frequency_min *= -k*T_unit/h                              # [k_BT_0/h]
    nu_grid = np.asarray(frequency_min)   # frequency group f interfaces [Hz]
    
    if nfreq > 2:
        try:
            if (athinput['problem']['frequency_table'] == 1):
                fname = "freq_table.txt"
            nu_grid = np.loadtxt(fname)
            if len(nu_grid)+1 != nfreq:
                raise ValueError(f'{fname} inconsistent with `nfreq` '\
                                 +'input parameter')
            if np.all(nu_grid < 0):  # unit switch: code (<0) or cgs (>0)
                nu_grid *= -k*T_unit/h                            # [k_BT_0/h]
        except KeyError:
            frequency_max = athinput['radiation']['frequency_max'] # [Hz]
            if frequency_max < 0:  # unit switch: code (<0) or cgs (>0)
                frequency_max *= -k*T_unit/h                       # [k_BT_0/h]
            nu_grid = np.logspace(np.log10(frequency_min),
                                  np.log10(frequency_max), nfreq-1)
    nu_grid = np.insert(nu_grid, 0, nu_min)
    nu_grid = np.append(nu_grid, nu_max)
    if nu_grid[0] < nu_min:
        warnings.warn('Lowest frequency group is below the lowest frequency '\
                      +'given by the opacity table')
    if nu_grid[-1] > nu_max:
        warnings.warn('Highest frequency group is above the highest frequency '\
                      +'given by the opacity table')
    i_nu1 = BinarySearchIncreasing(opac_freq, 0, len(opac_freq)-1, nu_grid[1])

for i in range(nfreq):
    kappa_pf_table[:, i] = PlanckMeanOpacities(opac_kabs[i_nu0:i_nu1],
                                               Bnu_table[:, i_nu0:i_nu1],
                                               opac_freq[i_nu0:i_nu1],
                                               temp_table)
    if scattering:
        kappa_rf_table[:, i] = RosselandMeanOpacities(opac_comb[i_nu0:i_nu1],
                                                      dBnu_dT_table[:, i_nu0:i_nu1],
                                                      opac_freq[i_nu0:i_nu1])
        kappa_sf_table[:, i] = RosselandMeanOpacities(opac_ksca[i_nu0:i_nu1],
                                                      dBnu_dT_table[:, i_nu0:i_nu1],
                                                      opac_freq[i_nu0:i_nu1])
    else:
        kappa_rf_table[:, i] = RosselandMeanOpacities(opac_kabs[i_nu0:i_nu1],
                                                      dBnu_dT_table[:, i_nu0:i_nu1],
                                                      opac_freq[i_nu0:i_nu1])
    i_nu0 = i_nu1
    if i < (nfreq - 2):  # intermediate frequency group
        i_nu1 = BinarySearchIncreasing(opac_freq, 0, len(opac_freq)-1,
                                       nu_grid[i+2])
    else:                # (next-to-)last frequency group
        i_nu1 = len(opac_freq)-1

# Convert units from cgs to code
temp_table /= T_unit                                             # [T_0]
kappa_pf_table *= density_unit*length_unit                       # [\rho_0*L_0]
kappa_rf_table *= density_unit*length_unit                       # [\rho_0*L_0]
kappa_sf_table *= density_unit*length_unit                       # [\rho_0*L_0]

# Save tables to text files for Athena++ input
np.savetxt('temp_table.txt', temp_table)
np.savetxt('kappa_pf_table.txt', kappa_pf_table)
np.savetxt('kappa_rf_table.txt', kappa_rf_table)
if scattering:
    np.savetxt('kappa_sf_table.txt', kappa_sf_table)
