% Used in: Extended Data Figure 4 (a-c), Extended Data Figure 6 (d)
% Plots instantaneous firing rates of TRN cells and cortical LFP
% Spikes (20 000 Hz sample rate), and cortical LFP waveform (downsampled to 500 Hz) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger

%%

clear all

filename_EEG='EEG_RTLVE_9_sejt_2_baseline_500Hz.mat'; % copy here .mat file for cortical LFP waveform
load(filename_EEG);
v_EEG=EEG.values;
sv_EEG=((v_EEG-mean(v_EEG))./std(v_EEG));

[pxx,f]=pwelch(sv_EEG,2^12,[],[],500)
plot(f,pxx,'b'); 
xlim([0 4])

hold on

FFTmax=max(pxx)

FFT_0_to_4Hz=pxx(1:33)';

clearvars FFTmax FFT_0_to_4Hz