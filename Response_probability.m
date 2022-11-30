
% Used in: Figure 1 (q, r), Figure 5 (j), Extended Data Figure 2 (f),
% Extended Data Figure 3 (b)
% Calculates response probability, response delay and jitter for optogenetic stimulation data. 
% Plots PSTH
% Spikes (20 000 Hz sample rate) and laser stimuli (20 000 Hz sample rate) 
% were detected in Spike2, and exported to .mat format
%Code by Nora Hadinger

%%

clear all

filename_spike='.mat'; % copy here .mat file for spike timestamps
load(filename_spike);
t_spikes=spikes.times; % in s
filename_stimulus='.mat'; % copy here .mat file for stimulus timestamps
load(filename_stimulus);
t_stimulus=stimulus.times; % in s

delay_threshold=0.03; % time interval in s (from start of the stimulus) within which a spike is counted as response.

% PSTH

for i=(1:length(t_spikes))
    
 MM(i,:)=t_stimulus;
    
end

for ii=(1:length(t_stimulus))
    
 NN(:,ii)=t_spikes;
    
end

delays_matrix=NN-MM;
x_range=(-0.1 :0.002:0.1); % -100 to 100 ms

figure (1)

histogram(delays_matrix,x_range)
set(gca,'TickDir','out')

% Response probability

responses=find(delays_matrix<=delay_threshold&delays_matrix>=0);
stim_with_responses=unique(MM(responses),'stable'); % it takes burst responses as single events
number_of_responses= length (stim_with_responses); 
response_probability=number_of_responses/length(t_stimulus);

% Response delay (in case of bursts for the first spike)

all_responses=MM(responses);                                     

redundancy_in_delay_matrix=(find (diff(all_responses)==0))+1;

response_delays=delays_matrix(responses);

response_delays(redundancy_in_delay_matrix)=NaN; 

response_delays=(response_delays(~isnan(response_delays))); %

response_delay=nanmean(response_delays)*1000; %in ms

% Jitter

jitter=std(response_delays)*1000; %in ms

%

clearvars -except jitter response_delay response_probability