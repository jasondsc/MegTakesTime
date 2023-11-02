#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 09:59:37 2022

@author: admin
"""

# Import sim functions
import numpy as np
import pandas as pd

from scipy.linalg import norm
#from itertools import repeat

from neurodsp.sim import sim_powerlaw
#from neurodsp.utils import set_random_seed

# Import function to compute power spectra
#from neurodsp.spectral import compute_spectrum

# Import utilities for plotting data
from neurodsp.utils import create_times
#from neurodsp.plts.spectral import plot_power_spectra
#from neurodsp.plts.time_series import plot_time_series


def sim_peak_oscillation(sig_ap, fs, freq, bw, height):
    """Simulate a signal with an aperiodic component and a specific oscillation peak.
    Parameters
    ----------
    sig_ap : 1d array
        The timeseries of the aperiodic component.
    fs : float
        Sampling rate of ``sig_ap``.
    freq : float
        Central frequency for the gaussian peak in Hz.
    bw : float
        Bandwidth, or standard deviation, of gaussian peak in Hz.
    height : float
        Relative height of the gaussian peak at the central frequency ``freq``.
        Units of log10(power), over the aperiodic component.
    Returns
    -------
    sig : 1d array
        Time series with desired power spectrum.
    Notes
    -----
    - This function creates a time series whose power spectrum consists of an aperiodic component
    and a gaussian peak at ``freq`` with standard deviation ``bw`` and relative ``height``.
    - The periodic component of the signal will be sinusoidal.
    Examples
    --------
    Simulate a signal with aperiodic exponent of -2 & oscillation central frequency of 20 Hz:
    >>> from neurodsp.sim import sim_powerlaw
    >>> fs = 500
    >>> sig_ap = sim_powerlaw(n_seconds=10, fs=fs)
    >>> sig = sim_peak_oscillation(sig_ap, fs=fs, freq=20, bw=5, height=7)
    """

    sig_len = len(sig_ap)
    times = create_times(sig_len / fs, fs)

    # Construct the aperiodic component and compute its Fourier transform
    # Only use the first half of the frequencies from the FFT since the signal is real
    sig_ap_hat = np.fft.fft(sig_ap)[0:(sig_len // 2 + 1)]

    # Create the range of frequencies that appear in the power spectrum since these
    # will be the frequencies in the cosines we sum below
    freqs = np.linspace(0, fs / 2, num=sig_len // 2 + 1, endpoint=True)

    # Construct the array of relative heights above the aperiodic power spectrum
    rel_heights = np.array([height * np.exp(-(lin_freq - freq) ** 2 / (2 * bw ** 2)) \
        for lin_freq in freqs])

    # Build an array of the sum of squares of the cosines to use in the amplitude calculation
    cosine_norms = np.array([norm(np.cos(2 * np.pi * lin_freq * times), 2) ** 2 \
        for lin_freq in freqs])

    # Build an array of the amplitude coefficients
    cosine_coeffs = np.array([\
        (-np.real(sig_ap_hat[ell]) + np.sqrt(np.real(sig_ap_hat[ell]) ** 2 + \
        (10 ** rel_heights[ell] - 1) * np.abs(sig_ap_hat[ell]) ** 2)) / cosine_norms[ell] \
        for ell in range(cosine_norms.shape[0])])

    # Add cosines with the respective coefficients and with a random phase shift for each one
    sig_periodic = np.sum(np.array([cosine_coeffs[ell] * \
                                   np.cos(2 * np.pi * freqs[ell] * times + \
                                          2 * np.pi * np.random.rand()) \
                          for ell in range(cosine_norms.shape[0])]), axis=0)

    sig = sig_ap + sig_periodic

    return sig



parameters_used=pd.DataFrame({'SubID' : [], 'peak_freq' : [], 'peak_band' : [], 'peak_height' : [], 'aper_expon' : [], 'group': []})

#print('Begin simulations for group 1:')

for i in range(250):
    # Set some general settings, to be used across all simulations
    fs = 500 # sampling rate
    n_seconds = 6 # length of recording 
    times = create_times(n_seconds, fs)
    freq=np.random.normal(9.82, 1.37)
    while freq <= 0:
        freq=np.random.normal(9.82,1.37)
    bw=np.random.normal(2.68, 1.64) 
    while bw <= 0:
        bw=np.random.normal(2.68, 1.64)  
    height=np.random.normal(0.57, 0.23)
    while height <= 0:                                                          # avoid getting negative values of peak height
        height=np.random.normal(0.57, 0.23)  
    expon= np.random.normal(-0.91, 0.15)
    # Simulate aperiodic element
    print(i)
    for j in range(100):
        signal = sim_powerlaw(n_seconds, fs, expon)
        # Simulate periodic element 
        signal = sim_peak_oscillation(signal, fs, freq, bw, height)
        np.savetxt('/Users/admin/Downloads/SimulatedData/alpha_new/0.2_250/subject_' + str(i) + '_simulation_6sec_' + str(j) +'.csv', signal, delimiter=",")

    parameters_used= parameters_used.append({'SubID': i, 'peak_freq':freq, 'peak_band':bw, 'peak_height':height, 'aper_expon': expon, 'group': 1}, ignore_index=True)
    

# print(parameters_used)
# parameters_used.to_csv('/Users/admin/Downloads/SimulatedData/low_effect_size_parameters_group1.csv')

print('Begin simulations for group 2:')
for i in range(250,500):
    # Set some general settings, to be used across all simulations
    fs = 500 # sampling rate
    n_seconds = 6 # length of recording 
    times = create_times(n_seconds, fs)
    freq=np.random.normal(9.82, 1.37)
    while freq <= 0:
        freq=np.random.normal(9.82,1.37)
    bw=np.random.normal(2.68, 1.64) 
    while bw <= 0:
        bw=np.random.normal(2.68, 1.64) 
    height=np.random.normal(0.616, 0.23)
    while height <= 0:                                                          # avoid getting negative values of peak height
        height=np.random.normal(0.616, 0.23)
    expon= np.random.normal(-0.91, 0.15)
    # Simulate aperiodic element
    print(i)
    for j in range(100):
        signal = sim_powerlaw(n_seconds, fs, expon)
        # Simulate periodic element 
        signal = sim_peak_oscillation(signal, fs, freq, bw, height)
        np.savetxt('/Users/admin/Downloads/SimulatedData/alpha_new/0.2_250/subject_' + str(i) + '_simulation_6sec_' + str(j) +'.csv', signal, delimiter=",")

    parameters_used= parameters_used.append({'SubID': i, 'peak_freq':freq, 'peak_band':bw, 'peak_height':height, 'aper_expon': expon, 'group': 2}, ignore_index=True)

print(parameters_used)

parameters_used.to_csv('/Users/admin/Downloads/SimulatedData/alpha_new/0.2_250/effect_size_parameters.csv')

# Conduct a two-sided t-test
parameters_used = pd.read_csv('/Users/admin/Downloads/SimulatedData/alpha_new/0.2_250/effect_size_parameters.csv')
parameters_used = parameters_used.drop(parameters_used.columns[[0]], axis=1)

mask = parameters_used['group'] == 1
group_1 = parameters_used[mask]
group_2 = parameters_used[~mask] 
#print (group_1)
print (group_2)

from scipy.stats import ttest_ind
peak_freq_ttest = ttest_ind(group_1['peak_freq'],group_2['peak_freq'])
peak_band_ttest = ttest_ind(group_1['peak_band'],group_2['peak_band'])
peak_height_ttest = ttest_ind(group_1['peak_height'],group_2['peak_height'])
aper_expon_ttest = ttest_ind(group_1['aper_expon'],group_2['aper_expon'])
print(peak_band_ttest, peak_freq_ttest, peak_height_ttest, aper_expon_ttest)
