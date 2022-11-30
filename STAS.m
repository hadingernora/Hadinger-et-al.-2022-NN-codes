% Used in: Figure 6 (i,j), Figure 7 (g-j), Extended Data Figure 4 (d), Extended Data Figure 6 (i)
% Calculates peak amplitudes of optogenetically evoked cortical LFP response averages 
% Plots average waveforms (STAs), and colorplots for single instantaneous cortical LFP tarces (100-100 ms) around the single spikes or the first spikes of the bursts.  
% Calculates STA peak amplitudes and lags between individual first spikes per bursts and the corresponding individual instanteneous LFP peaks. 
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

%% STAs

EEG_window_halfwidth1=1000; % in ms *10, here 100 ms
EEG_window_halfwidth2=1000; % in ms *10, here 100 ms


% arranges bursts according to their aIBFs (descending order)

e1=first_spike_time_correction_vector(find (baseline_intraburst_frequency <=900 &  baseline_intraburst_frequency>=100))';
e2=baseline_intraburst_frequency(find (baseline_intraburst_frequency <=900 &  baseline_intraburst_frequency>=100));
mat_aIBF_firstspiketimecorrection=[e2, e1];
[~,idx] = sort(mat_aIBF_firstspiketimecorrection(:,1),'descend'); 
sortedmat = mat_aIBF_firstspiketimecorrection(idx,:); 
spike_trigger=sortedmat(:,2)


STA_matrix=zeros(length(spike_trigger),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(spike_trigger)
    
reference_point=spike_trigger(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG)
    
STA_matrix (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

else

STA_matrix (i,:)=NaN;    
    
end

end


fourth=find(e1>prctile(e1,75));
STA_vector_lfburst=nanmean(STA_matrix(fourth,:));
first=find(e1<prctile(e1,25));
third=find(e1<prctile(e1,75)& e1>prctile(e1,50));
STA_vector_third=nanmean(STA_matrix(third,:));
second=find(e1<prctile(e1,50)& e1>prctile(e1,25));
STA_vector_second=nanmean(STA_matrix(second,:));
STA_vector_hfburst=nanmean(STA_matrix(first,:));

% single spikes


ss_trigger=first_spike_time_correction_vector(isnan(baseline_intraburst_frequency));

STA_matrix_ss=zeros(length(ss_trigger),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(ss_trigger)
    
reference_point=ss_trigger(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG)
    
STA_matrix_ss (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

else

STA_matrix_ss (i,:)=NaN;    
    
end

end

STA_vector_ss=nanmean(STA_matrix_ss);

% STA plots for single spikes and the four burst categories

figure (1)

x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));

plot(x_val,STA_vector_ss)
hold on
plot(x_val,STA_vector_lfburst)
hold on
plot(x_val,STA_vector_third)
hold on
plot(x_val,STA_vector_second)
hold on
plot(x_val,STA_vector_hfburst)
hold on

legend('single spikes','lf bursts', '2nd 25%', '3rd 25%','hf bursts');
xlabel('ms');
ylabel('au'); % arbitrary unit


% colorplots

r1=size(STA_matrix(fourth,:),1);
rand1=randsample((1:(size(STA_matrix_ss,1))),r1); % randomly selects equal events for a better visual comparism 

figure (2)
set(gcf,'Position', [400 10 2000 5000])

x=[-100 100]
y=[1 r1]

subplot (1,6,1)
imagesc(x,y, STA_matrix_ss(rand1,:))
title ('single spikes')
xlabel('ms');
ylabel('#');
subplot(1,6,2)
imagesc(x,y,STA_matrix(fourth,:))
title ('lf bursts')
xlabel('ms');
ylabel('#');
subplot(1,6,3)
imagesc(x,y,STA_matrix(third,:))
title ('2nd 25%')
xlabel('ms');
ylabel('#');
subplot(1,6,4)
imagesc(x,y,STA_matrix(second,:))
title ('3nd 25%')
xlabel('ms');
ylabel('#');
subplot(1,6,5)
imagesc(x,y,STA_matrix(first,:))
title ('hf bursts')
xlabel('ms');
ylabel('#');
c=colorbar;
c.Position = [0.772 0.11 0.022 0.815];

caxis([-(3*(0.4491+1.0850)) 3*(0.4491+1.0850)]); %sets colormap axis

% STA peak amplitudes for lf and hf bursts

STA_peak_ampl_hf=max(STA_vector_hfburst);
STA_peak_ampl_lf=STA_vector_lfburst(find(STA_vector_hfburst==max(STA_vector_hfburst)));

% lags between individual first spikes per bursts and the corresponding individual instanteneous LFP peaks

for ii=1:length(fourth)
       
sta_peak_position_actual_1=find(max(STA_matrix(fourth(ii),:))==STA_matrix(fourth(ii),:));

if length (sta_peak_position_actual_1) >0

sta_peak_position_1(ii)=sta_peak_position_actual_1(1);
sta_peak_aplitude_actual_1=max(STA_matrix(fourth(ii),:));
sta_peak_aplitude_1(ii)=sta_peak_aplitude_actual_1(1);


else
    
sta_peak_position_1(ii)=NaN;
sta_peak_aplitude_1(ii)=NaN;   
    
end

end

average_lag_lf=-((nanmean(sta_peak_position_1)-1000)/10);

for iii=1:length(first)
       
sta_peak_position_actual_2=find(max(STA_matrix(first(iii),:))==STA_matrix(first(iii),:));

if length (sta_peak_position_actual_2) >0

sta_peak_position_2(iii)=sta_peak_position_actual_2(1);
sta_peak_aplitude_actual_2=max(STA_matrix(first(iii),:));
sta_peak_aplitude_2(iii)=sta_peak_aplitude_actual_2(1);


else
    
sta_peak_position_2(iii)=NaN;
sta_peak_aplitude_2(iii)=NaN;   
    
end

end

average_lag_hf=-((nanmean(sta_peak_position_2)-1000)/10);

%

clearvars -except STA_vector_lfburst STA_vector_hfburst STA_peak_ampl_lf STA_peak_ampl_hf average_lag_lf average_lag_hf 