% Used in: Figure 6 (k), Extended Data Figure 4 (b)
% Plots instantaneous firing rates of TRN cells and cortical LFP
% Spikes (20 000 Hz sample rate), and cortical LFP waveform (downsampled to 10 000 Hz) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger

%%

clear all

filename_spike='all_spikes_RTLVE_9_sejt_2_baseline.mat'; % copy here .mat file for spike timestamps
load(filename_spike);
t_spikes=spikes.times; % in s
t_spikes_int=int32(t_spikes*10000)

filename_EEG='EEG_RTLVE_9_sejt_2_baseline_newvers.mat'; % copy here .mat file for cortical LFP waveform
load(filename_EEG);
v_EEG=EEG.values;
sv_EEG=((v_EEG-mean(v_EEG))./std(v_EEG)); %standardized EEG (Z-score)


diffiring=diff(t_spikes_int)
diffiring=double(diffiring)
diffiring=diffiring./80
plot(t_spikes_int(2:end),1./diffiring,'b');
ylim([-5 10])

hold on

hold on

plot(sv_EEG,'k')

hold on

z_treshold=0.8*ones(length(1:10000:length(sv_EEG)));%100 Hz

plot ((1:10000:length(sv_EEG)),z_treshold,'--r')

hold on

clearvars
