# MEA_analysis
Source code for analyzing multi-electrode array signal

Code is consists of 2 parts: Analysis of spike/burst and analysis of LFP event

Before the analysis, you may export each MEA channel response as .csv, read files from MATLAB, and save the data as .mat

# Analysis for spike/burst detection
"MEA_sample_data.mat" is for spike detection

Run v220523_MEA_SpikeAnalysis.m to run the analysis

To briefly explain the process:

1) Filter each electrode's signal with a gaussian filter 
2) Apply Notch filter
3) Set up a threshold with +N*(STD of signal) and -N*(STD of signal), and detect the time point when signal exceeds threshold. 
For threshold, typically ~5 values are used
4) Burst was defined with a series of spikes with ISI < 10ms (inter-spike-interval), the ISI threshold may vary across conditions
5) Firing rate, Burst rate, Number of burst spikes in a single burst, and ratio of burst spikes among all spikes were calculated.

# Analysis for LFP event detection
"MEA_sample_data_for_LFPevent.mat" is for LFP event detection

Run v220523_MEA_LFPAnalysis.m

To briefly explain the process:

1) Notch filter was applied
2) Low pass filter (150Hz, 4th order Chebyshev filter) was applied
3) Set up a negative threshold with -N*(STD of signal), and the signal exceeding the negative threshold was defined as LFP event. For threshold, 4Z was used.
4) In each electrode, all LFP events were averaged to measure shape of LFP event.
5) Amplitude of LFP event was defined as a difference of minimum and maximum signal.
6) The duration of LFP was defined as time between significantly fluctuating signals within -200ms ~ 200ms from LFP event;
significant fluctuation was defined when p<0.001 compared to a randomly sampled signal (bootstrap method)
