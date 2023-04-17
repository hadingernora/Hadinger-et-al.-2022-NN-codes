% Used in: Extended Data Figure 5, Extended Data Figure 6 (f)
% Detects fast cortical LFP evets, generates PSTHs. Plots fast evet waveform.   
% Spikes (20 000 Hz sample rate), and cortical LFP waveform (downsampled to 10 000 Hz) 
% were detected in Spike2, and exported to .mat format
% Code by Nora Hadinger

%%

clear all

filename_spike='all_spikes_RTLVE_9_sejt_2_baseline.mat'; % copy here .mat file for spike timestamps
load(filename_spike);
t_spikes=spikes.times; % in s

filename_EEG='EEG_RTLVE_9_sejt_2_baseline_newvers.mat'; % copy here .mat file for cortical LFP waveform
load(filename_EEG);
v_EEG=EEG.values;
sv_EEG=((v_EEG-mean(v_EEG))./std(v_EEG)); %standardized EEG (Z-score)

% Calculates linear fit to randomly selected 30ms long LFP segments 

trigger=randsample((1:length(sv_EEG)), 200000);

EEG_window_halfwidth1=15; 
EEG_window_halfwidth2=15; 

x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));

STA_matrix=zeros(length(trigger),((EEG_window_halfwidth1+EEG_window_halfwidth2)+1));

for i=1:length(trigger)
    
reference_point=trigger(i);

if reference_point>=EEG_window_halfwidth1 && reference_point+EEG_window_halfwidth2<=length(sv_EEG)
    
STA_matrix (i,:)=sv_EEG(reference_point-(EEG_window_halfwidth1):reference_point+(EEG_window_halfwidth2));

else

STA_matrix (i,:)=NaN;    
    
end

end

    R_egyes_gorbek=[];
    slope=[];
    
    
    for i=1:size(STA_matrix,1)
        
    P = polyfit(x_val,STA_matrix(i,:),1);
    slope(i)=P(1);  
        
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

i

    end
  
% Selects the upper 2% of linear fits according to their steepness 
    
steepest_fits=prctile(slope,98)

TT=(slope>steepest_fits);

kk=trigger(TT);

sk=sort(kk);

finddiff=find(diff(sk)>500);

sharp_events=sk(finddiff);

sharp_events_pos=sharp_events(find(sv_EEG(sharp_events)>-2));

% Finds the peaks of the sharp events

halfwidth=200;

sharp_events_pos=sharp_events_pos(sharp_events_pos>halfwidth);

for sh=(1:length (sharp_events_pos))
 
    sharp_events_wave=(sv_EEG(sharp_events_pos(sh)-halfwidth:sharp_events_pos(sh)+halfwidth))';
    peak_sharp_events(sh)= max(sharp_events_wave);
    peak_position_rel(sh)=find(sharp_events_wave==peak_sharp_events(sh),1); 
    peak_position_b(sh)= sharp_events_pos(sh)+(peak_position_rel(sh)-(halfwidth+1));
    
end   


high_ampl_peaks=find(peak_sharp_events> mean(sv_EEG)+3*std(sv_EEG));
high_amp_peaks_pos=peak_position_b(high_ampl_peaks)
peak_position_b=high_amp_peaks_pos; 

% plots cortical LFP, fast event peaks (green), and TRN spikes (red)

figure (1)

plot(sv_EEG,'k')
hold on
plot(peak_position_b, sv_EEG(peak_position_b),'g.','MarkerSize',15)
plot(int32(t_spikes*10000),sv_EEG(int32(t_spikes*10000)),'r.','MarkerSize',15)

legend('cortical_LFP', 'fast_event_peak','TRN_spike')

% PSTH

clear('psth_vector', 'psth_vector_position', 'psth_vector_position_sum', 'psth_vector_delay',' psth_vector_delay_sum','t');
t=int32(t_spikes*10000);

psth_window=2000 % 200 ms
psth_vector_delay_sum=[];
psth_vector_position_sum=[];

peak_position=unique(peak_position_b); % valamikor 2x veszi ugyanazt a peaket az el?z? szekcióban, még nem jöttem rá, miért, de ezzel így kevésbé elegánsan megoldottam

for shh=1:length(peak_position)
    
   psth_vector_position= t(find (t>peak_position(shh)-psth_window & t<peak_position(shh)+psth_window));
   psth_vector_position_sum=[psth_vector_position_sum;psth_vector_position]
   psth_vector_delay=psth_vector_position-peak_position(shh);
   psth_vector_delay_sum=[psth_vector_delay_sum;psth_vector_delay];
    
end

figure (2)

H01.histogr=histogram (psth_vector_delay_sum','BinWidth',100,'Normalization','probability') % 10 ms bins
xlim([-2000 2000]) % -200 to 200 ms
legend('PSTH')

% plos spikes around fast event peaks (-200t o 200 ms)

figure (3)

plot (psth_vector_delay_sum',1:length(psth_vector_delay_sum),'.')
legend('spikes around peaks')

PSTH_values=H01.histogr.Values

% Fast event waveform (-200 ms to 200 ms)

EEG_window_halfwidth1=2000; %200 ms 
EEG_window_halfwidth2=2000; % 200 ms

trigger=peak_position_b

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
STA_SD=std(STA_matrix);

x_val=(-(EEG_window_halfwidth1/10):0.1:(EEG_window_halfwidth2/10));
plot(x_val,STA_vector,'b')
legend('fast event waveform')

clearvars -except PSTH_values
