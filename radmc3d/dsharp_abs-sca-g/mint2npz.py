#!/usr/bin/env python3
#===============================================================================
# mint2npz.py
#
# Uses the radmc3dPy module to read `./mean_intensity.out`, integrate it over
# all (decreasing) frequencies, and save it to `./total_mean_intensity.npz`.
#
# Author: Stanley A. Baronett
# Created: 2024-07-25
# Updated: 2024-07-27
#===============================================================================
import numpy as np
from radmc3dPy import *
from scipy import integrate

data = analyze.readData(mint=True)
print('Integrating over all frequencies', flush=True)
totmeanint = np.abs(integrate.trapezoid(data.meanint, x=data.freq, axis=-1))
print('Saving to total_mean_intensity.npz', flush=True)
np.savez_compressed('total_mean_intensity', totmeanint=totmeanint)
