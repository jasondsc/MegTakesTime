# Import sim functions
import numpy as np
import pandas as pd

from scipy.linalg import norm
#from itertools import repeat

from neurodsp.sim import sim_powerlaw
from neurodsp.sim import sim_combined
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
    
    freq1=freq[0]
    freq2=freq[1]
    bw1=bw[0]
    bw2=bw[1]
    height1=height[0]
    height2=height[1]

    # Construct the aperiodic component and compute its Fourier transform
    # Only use the first half of the frequencies from the FFT since the signal is real
    sig_ap_hat = np.fft.fft(sig_ap)[0:(sig_len // 2 + 1)]

    # Create the range of frequencies that appear in the power spectrum since these
    # will be the frequencies in the cosines we sum below
    freqs = np.linspace(0, fs / 2, num=sig_len // 2 + 1, endpoint=True)

    # Construct the array of relative heights above the aperiodic power spectrum
    rel_heights1 = np.array([height1 * np.exp(-(lin_freq - freq1) ** 2 / (2 * bw1 ** 2)) \
        for lin_freq in freqs])
        
    # Construct the array of relative heights above the aperiodic power spectrum
    rel_heights2 = np.array([height2 * np.exp(-(lin_freq - freq2) ** 2 / (2 * bw2 ** 2)) \
            for lin_freq in freqs])
        
    rel_heights= rel_heights1+rel_heights2

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

#def freq_freq(freq):
    #return freq([freq_beta,freq_alpha])

parameters_used=pd.DataFrame({'SubID' : [], 'peak_freq_beta' : [], 'peak_band_beta' : [], 'peak_height_beta' : [], 'peak_freq_alpha' : [], 'peak_band_alpha' : [], 'peak_height_alpha' : [], 'aper_expon' : [], 'group': []})

#print('Begin simulations for group 1:')

for i in range(250):
    # Set some general settings, to be used across all simulations
    fs = 500 # sampling rate
    n_seconds = 6 # length of recording 
    times = create_times(n_seconds, fs)
    
    freq_beta=np.random.normal(19.33, 2.86)
    freq_alpha=np.random.normal(9.82, 1.37)
    
    bw_beta=np.random.normal(4.98, 1.07) 
    while bw_beta <= 0:
        bw_beta=np.random.normal(4.98, 1.07)
    bw_alpha=np.random.normal(2.09, 0.93) 
    while bw_alpha <= 0:
        bw_alpha=np.random.normal(2.09, 0.93) 
        
    height_beta=np.random.normal(0.48, 0.16)
    while height_beta <= 0:                                                          #avoid getting negative values of peak height
        height_beta=np.random.normal(0.48, 0.16)
    height_alpha=np.random.normal(0.57, 0.23)
    while height_alpha <= 0:                                                          #avoid getting negative values of peak height
        height_alpha=np.random.normal(0.57, 0.23)
    
    expon= np.random.normal(-0.91, 0.15)
    while expon >= 0:                                                          #avoid getting negative values of peak height
        expon=np.random.normal(-0.91, 0.15)
    # Simulate aperiodic element
    #freq = freq_beta, freq_beta
    #bw = bw_beta, bw_alpha
    #height = height_beta, height_alpha
    print(i)
    for j in range(100):
        signal = sim_powerlaw(n_seconds, fs, expon)
        #component = sim_components = {'sim_powerlaw': {'exponent':-2}, 'sim_oscillation': [{'freq': 10}]}
        # Simulate periodic element 
        signal = sim_peak_oscillation(signal, fs, [freq_beta,freq_alpha], [bw_beta,bw_alpha], [height_beta, height_alpha])
        np.savetxt('/Users/admin/Desktop/test/subject_' + str(i) + '_simulation_6sec_' + str(j) +'.csv', signal, delimiter=",")

    parameters_used= parameters_used.append({'SubID': i, 'peak_freq_beta':freq_beta, 'peak_band_beta':bw_beta, 'peak_height_beta':height_beta, 'peak_freq_alpha':freq_alpha, 'peak_band_alpha':bw_alpha, 'peak_height_alpha':height_alpha,'aper_expon': expon, 'group': 1}, ignore_index=True)
    
parameters_used.to_csv('/Users/admin/Desktop/test/low_effect_size_TWOPEAKS.csv')

# print(parameters_used)
# parameters_used.to_csv('/Users/admin/Downloads/SimulatedData/low_effect_size_parameters_group1.csv')

print('Begin simulations for group 2:')
for i in range(250,500):
    # Set some general settings, to be used across all simulations
    fs = 500 # sampling rate
    n_seconds = 6 # length of recording 
    times = create_times(n_seconds, fs)
    
    freq_beta=np.random.normal(19.33, 2.86)
    freq_alpha=np.random.normal(9.82, 1.37)
    
    bw_beta=np.random.normal(4.98, 1.07) 
    while bw_beta <= 0:
        bw_beta=np.random.normal(4.98, 1.07)
    bw_alpha=np.random.normal(2.09, 0.93) 
    while bw_alpha <= 0:
        bw_alpha=np.random.normal(2.09, 0.93) 
        
    height_beta=np.random.normal(0.48, 0.16)
    while height_beta <= 0:                                                          #avoid getting negative values of peak height
        height_beta=np.random.normal(0.48, 0.16)
    height_alpha=np.random.normal(0.57, 0.23)
    while height_alpha <= 0:                                                          #avoid getting negative values of peak height
        height_alpha=np.random.normal(0.57, 0.23)

    expon= np.random.normal(-0.88, 0.15)
    while expon >= 0:                                                          #avoid getting negative values of peak height
        expon=np.random.normal(-0.88, 0.15)
    # Simulate aperiodic element
    print(i)
    for j in range(100):
        signal = sim_powerlaw(n_seconds, fs, expon)
        # Simulate periodic element 
        signal = sim_peak_oscillation(signal, fs, [freq_beta,freq_alpha], [bw_beta,bw_alpha], [height_beta, height_alpha])
        np.savetxt('/Users/admin/Desktop/test/subject_' + str(i) + '_simulation_6sec_' + str(j) +'.csv', signal, delimiter=",")

    parameters_used= parameters_used.append({'SubID': i, 'peak_freq_beta':freq_beta, 'peak_band_beta':bw_beta, 'peak_height_beta':height_beta, 'peak_freq_alpha':freq_alpha, 'peak_band_alpha':bw_alpha, 'peak_height_alpha':height_alpha,'aper_expon': expon, 'group': 2}, ignore_index=True)
   

print(parameters_used)

parameters_used.to_csv('/Users/admin/Desktop/test/low_effect_size_TWOPEAKS.csv')

# Conduct a two-sided t-test
parameters_used = pd.read_csv('/Users/admin/Desktop/test/low_effect_size_TWOPEAKS.csv')
parameters_used = parameters_used.drop(parameters_used.columns[[0]], axis=1)

mask = parameters_used['group'] == 1
group_1 = parameters_used[mask]
group_2 = parameters_used[~mask] 
#print (group_1)
print (group_2)

from scipy.stats import ttest_ind
peak_freq_ttest_beta = ttest_ind(group_1['peak_freq_beta'],group_2['peak_freq_beta'])
peak_freq_ttest_alpha = ttest_ind(group_1['peak_freq_alpha'],group_2['peak_freq_alpha'])
peak_band_ttest_beta = ttest_ind(group_1['peak_band_beta'],group_2['peak_band_beta'])
peak_band_ttest_alpha = ttest_ind(group_1['peak_band_alpha'],group_2['peak_band_alpha'])
peak_height_ttest_beta = ttest_ind(group_1['peak_height_beta'],group_2['peak_height_beta'])
peak_height_ttest_alpha = ttest_ind(group_1['peak_height_alpha'],group_2['peak_height_alpha'])
aper_expon_ttest = ttest_ind(group_1['aper_expon'],group_2['aper_expon'])
print(peak_freq_ttest_beta, peak_freq_ttest_alpha, peak_band_ttest_beta, peak_band_ttest_alpha, peak_height_ttest_beta, peak_height_ttest_alpha, aper_expon_ttest)
