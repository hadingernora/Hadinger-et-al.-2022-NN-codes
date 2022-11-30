% Used in: Figure 6 (d, e), Extended Data Figure 3 (d-i),
% Calculates fraction of bursts, average intraburst frequencies (aIBF) and numer of spikes / bursts for evoked responses. 
% Spikes (20 000 Hz sample rate) and laser stimuli (20 000 Hz sample rate) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger

%%

clear all

filename_spike='.mat'; % copy here .mat file for spike timestamps
load(filename_spike);
t_spikes=spikes.times; % in s
filename_stimulus='.mat'; % copy here .mat file for stimulus timestamps
load(filename_stimulus);
t_stimulus=stimulus.times; % in s

delay_threshold=0.03; % time interval in s (from start of the stimulus) within which a spike is counted as response.
bk=100; %burst criteria (ISI in ms*10)
first_ISI_criteria=0.001; %(ms)

% Detecting evoked responses

for i=(1:length(t_spikes))
    
 MM(i,:)=t_stimulus;
    
end

for ii=(1:length(t_stimulus))
    
 NN(:,ii)=t_spikes;
    
end

delays_matrix=NN-MM;

responses=find(delays_matrix<=delay_threshold&delays_matrix>=0);
stim_with_responses=unique(MM(responses),'stable'); % it takes burst responses as single events

% Selecting first spikes per response


all_responses=MM(responses);                                     

redundancy_in_response_matrix=(find (diff(all_responses)==0))+1;

first_responses=NN(responses) ;

first_responses(redundancy_in_response_matrix)=NaN;

first_responses=(first_responses(~isnan(first_responses))); % selects first spikes of responses

% Delineates burst and single spike responses

t_spikes_digital = zeros (1,max(int32(t_spikes*10000)));% creates a vector where each time points of the recording have the value of zero 
t_spikes_digital(int32(t_spikes*10000))=1; % renders 1 to the timepoints of all spikes

first_responses_digital=int32(first_responses*10000);


for iii=(1:(length(first_responses)))

FR=first_responses_digital(iii);

nr=0;
nr2=1;

while nr<bk 

if t_spikes_digital(FR+(nr2-1))==1
response_matrix(nr2,iii)=1;
nr=0;

else
response_matrix(nr2,iii)=0;
nr=nr+1;
end

if FR==length(t_spikes_digital)
nr=bk+1;
    
else
nr2=nr2+1;

end
end
end



% Spikes per bursts

spikes_per_responses=sum(response_matrix);

spikes_per_bursts=mean(spikes_per_responses(find(spikes_per_responses>1)));

% Fraction of bursts

fraction_of_bursts=length(find(spikes_per_responses>1))/length(spikes_per_responses);

% aIBF

burst_structure_matrix=NaN(length(first_responses),max(spikes_per_responses));

for iiii=(1:length(spikes_per_responses))
 
burst_structure_matrix(iiii, (1:(sum(response_matrix(:,iiii)))))=find(response_matrix(:,iiii)>0)/10000;
    
end

for iiiii=(1:length(spikes_per_responses))

ISI_burst_matrix(iiiii,:)=diff(burst_structure_matrix(iiiii,:));% s

ISI_burst_average(iiiii)=nanmean(nonzeros (ISI_burst_matrix(iiiii,:)));

end

intraburst_frequency_matrix=1./ISI_burst_matrix; 
intraburst_frequency_average= nanmean(intraburst_frequency_matrix,2);


aIBF= nanmean (intraburst_frequency_average); %Hz

%

clearvars -except aIBF fraction_of_bursts spikes_per_bursts



