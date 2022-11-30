% Used in: Figure 6 (k-m), Figure 7 (e,f), Extended Data Figure 4 (a-b), Extended Data Figure 4 (j,k)
% Calculates correlation between aIBF of TRN bursts and the slope of the instantaneous LFP peak   
% Spikes (20 000 Hz sample rate), and cortical LFP waveform (downsampled to 10 000 Hz) 
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

% aIBF

baseline_spikes_per_event=sum(baseline_bursts_matrix); % takes into account but bursts anfd single spike events

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

% STA for all butsts

EEG_window_halfwidth1=1000; % in ms *10, here 100 ms
EEG_window_halfwidth2=1000; % in ms *10, here 100 ms

spike_trigger=first_spike_time_correction_vector(find (baseline_intraburst_frequency <=900 &  baseline_intraburst_frequency>=100))';

STA_matrix=zeros(length(spike_trigger),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(spike_trigger)
    
reference_point=spike_trigger(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG)
    
STA_matrix (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

else

STA_matrix (i,:)=NaN;    
    
end

end

STA_vector=nanmean(STA_matrix);

STA_max=find(STA_vector==max(STA_vector));

x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));

x_STA_max_ms_all_bursts=(x_val(STA_max)); % average lag between the first spikes and the average STA peak for all bursts

%% Calculates instantaneous LFP slope and correlates it to the aIBF of the corresponding burst

jjjj=250-(x_STA_max_ms_all_bursts*10);
jjj=50+(x_STA_max_ms_all_bursts*10);

spike_trigger=first_spike_time_correction_vector(find (baseline_intraburst_frequency<=900 &  baseline_intraburst_frequency>=100));

EEG_window_halfwidth1=jjjj; 
EEG_window_halfwidth2=jjj; 
x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));

R_sum=0;
P_sum=0;

R2_sum=0;
P2_sum=0;

for iter1=1 % for all but ArchT experiments
%for iter= 1:200 (in case of ArchT experiments, where a portion of bursts were randomly selected for each group)    
    
%spike_trigger_actual= randsample(spike_trigger,number); % in case of ArchT
%experiments, equal number (number) of bursts were selected for baseline and ArchT
%condition

spike_trigger_actual=spike_trigger;% taking all the bursts (all but ArchT experiments)
spike_trigger_actual=intersect(spike_trigger_actual,first_spike_time_correction_vector)

STA_matrix=zeros(length(spike_trigger_actual),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(spike_trigger_actual)
 
    
reference_point=spike_trigger_actual(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG);
    
STA_matrix (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

shift=double(reference_point-(EEG_window_halfwidth1));

iksz=(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));
       
end

end

intraburst_frequency_average_vector=[];
egy_per_intraburst_frequency_average_vector=[];
spikes_per_bursts_vector=[];
bursts_duration_vector=[];


    R_egyes_gorbek=[];
    slope=[];
    
    
    for i=1:size(STA_matrix,1)
        
    P = polyfit(x_val,STA_matrix(i,:),1);
    slope(i)=P(1);  
    
    intraburst_frequency_average_vector(i)=baseline_intraburst_frequency(find(first_spike_time_correction_vector==spike_trigger_actual(i)));
    egy_per_intraburst_frequency_average_vector(i)=1./intraburst_frequency_average_vector(i);
    spikes_per_bursts_vector(i)=baseline_spikes_per_event(find(first_spike_time_correction_vector==spike_trigger_actual(i)));
    bursts_duration_vector(i)=egy_per_intraburst_frequency_average_vector(i).*spikes_per_bursts_vector(i);
    
    
f = polyval(P, x_val);
Bbar = mean(STA_matrix(i,:));
SStot = sum((STA_matrix(i,:) - Bbar).^2);
SSreg = sum((f - Bbar).^2);
SSres = sum((STA_matrix(i,:) - f).^2);
R2 = 1 - SSres/SStot;
R_egyes_gorbek_pre = corrcoef(x_val,STA_matrix(i,:));
R_egyes_gorbek(i)=R_egyes_gorbek_pre(2);
R_egyes_gorbek_mean=mean(abs(R_egyes_gorbek));
R_egyes_gorbek_SD=std(abs(R_egyes_gorbek));
   
    end

    
[R,P] = corrcoef(slope,intraburst_frequency_average_vector);
    
R_sum=(R_sum+R);
P_sum=(P_sum+P);

R=R_sum/iter1
p=P_sum/iter1

hold on

end


figure(1)

scatter(intraburst_frequency_average_vector,slope,25,'k','.')

%

clearvars -except R p
