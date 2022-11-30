% Used in: Figure 6 (d, h, k), Figure 7 (d)  Extended Data Figure 3
% (h,i),Extended Data Figure 4 (e-h), Extended Data Figure 6 (a-c)
% Calculates firing rate, fraction of bursts, average intraburst frequencies (aIBF) -for all burst one by one and their average-, and numer of spikes / bursts for spontaneous firing. 
% Spikes (20 000 Hz sample rate) and cortical LFP waveform (downsampled to 10 000 Hz) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger

%%

clear all

filename_spike='.mat'; % copy here .mat file for spike timestamps
load(filename_spike);
t_spikes=spikes.times; % in s

filename_EEG='.mat'; % copy here .mat file for cortical LFP waveform
load(filename_EEG);
v_EEG=EEG.values;
sv_EEG=((v_EEG-mean(v_EEG))./std(v_EEG)); %standardized EEG (Z-score)

bk=100; %burst criteria (ISI in ms*10)
first_ISI_criteria=0.001; %(ms)

t_spikes_digital = zeros (1,max(int32(t_spikes*10000))); 
t_spikes_digital(int32(t_spikes*10000))=1; 
 

%generates a matrix for bursts and single spikes (spikes outside bursts) 
 
spikecount=int32(t_spikes(1)*10000);
bnr3=1;

while spikecount<max(int32(t_spikes*10000))
                
bnr=0;
bnr2=1;

while bnr<bk
    
    
if t_spikes_digital(spikecount+(bnr2-1))==1
        
baseline_bursts_matrix(bnr2,bnr3)=1;
bnr=0;

else
baseline_bursts_matrix(bnr2,bnr3)=0;
bnr=bnr+1;
end

if spikecount+(bnr2-1)== max(int32(t_spikes*10000))
bnr=bk+1;

else
bnr2=bnr2+1;

end
end

first_spike_time_correction_vector(bnr3)=spikecount;

bnr3=bnr3+1;

if int32(t_spikes(sum(sum(baseline_bursts_matrix)))*10000)==max(int32(t_spikes*10000))
spikecount=max(int32(t_spikes*10000))+1;
    
else    
spikecount=int32(t_spikes(sum(sum(baseline_bursts_matrix))+1)*10000);
end


end

% Spikes per bursts

baseline_spikes_per_event=sum(baseline_bursts_matrix); % takes into account but bursts anfd single spike events
baseline_spikes_per_burst=baseline_spikes_per_event(baseline_spikes_per_event~=1); % takes into account only bursts

t_spikes_int=int32(t_spikes*10000);

t_spikes_int=unique(t_spikes_int);


for fgh=1:length(first_spike_time_correction_vector)
      
first_ISI_second_spike_vector(fgh)=t_spikes_int(find(t_spikes_int==first_spike_time_correction_vector(fgh))+1);

end

first_ISI_vector_baseline=double(first_ISI_second_spike_vector-first_spike_time_correction_vector)/10000;
average_first_ISI_vector_baseline=mean(first_ISI_vector_baseline(find(first_ISI_vector_baseline<=first_ISI_criteria)));

spikes_per_bursts=mean(baseline_spikes_per_burst);

% aIBF

baseline_bursts_structure_matrix=NaN(length(first_spike_time_correction_vector),max(baseline_spikes_per_event));

for biiiii=(1:length(baseline_spikes_per_event))
 
baseline_bursts_structure_matrix(biiiii, (1:(sum(baseline_bursts_matrix(:,biiiii)))))=find(baseline_bursts_matrix(:,biiiii)>0)/10000;
    
end

for biiiiii=(1:length(baseline_spikes_per_event))

baseline_ISI_burst_matrix(biiiiii,:)=diff(baseline_bursts_structure_matrix(biiiiii,:));% sec

baseline_ISI_burst_average(biiiiii)=nanmean((baseline_ISI_burst_matrix(biiiiii,:)));

end

baseline_intraburst_frequency_matrix=1./baseline_ISI_burst_matrix; 

baseline_intraburst_frequency= nanmean(baseline_intraburst_frequency_matrix,2);

aIBF= nanmean (baseline_intraburst_frequency); % Hz

% Firing rate

firing_rate=length(t_spikes)/(length(sv_EEG)/10000); % Hz


% Fraction of bursts

fraction_of_bursts=length(baseline_spikes_per_burst)/length(baseline_spikes_per_event);

%

clearvars -except firing_rate aIBF fraction_of_bursts spikes_per_bursts baseline_intraburst_frequency