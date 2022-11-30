% Used in: Figure 5 (k)
% Calculates peak amplitudes of optogenetically evoked cortical LFP response averages 
% Plots optogenetically evoked cortical LFP response average waveforms
% Laser stimuli (20 000 Hz sample rate), and cortical LFP waveform (downsampled to 10 000 Hz) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger
%%

clear all

filename_stimulus='.mat'; % copy here .mat file for stimulus timestamps
load(filename_stimulus);
t_stimulus=stimulus.times; % in s

filename_EEG='.mat'; % copy here .mat file for cortical LFP waveform
load(filename_EEG);
v_EEG=EEG.values;
sv_EEG=((v_EEG-mean(v_EEG))./std(v_EEG)); %standardized EEG (Z-score)

% Response average waveform

EEG_window_halfwidth1=1000; % 100 ms at 10 000 Hz sample rate
EEG_window_halfwidth2=1000; % 100 ms at 10 000 Hz sample rate

trigger=t_stimulus*10000;

STA_matrix=zeros(length(trigger),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(trigger)
    
reference_point=trigger(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG)
    
STA_matrix (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

else

STA_matrix (i,:)=NaN;    
    
end

end

STA_vector=nanmean(STA_matrix);

figure(1)

x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));
plot(x_val,STA_vector)

% Response average max amplitude

STA_max_amplitude=max(STA_vector); % in arbitrary unit (au) 

clearvars -except STA_max_amplitude

