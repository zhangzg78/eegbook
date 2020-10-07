%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to calculate ERS/ERD of multiple  
%               trials of laser-evoked potentials (LEP).
%
% Note:         This script is a supplementary file of Chapter6 in the book 
%               "EEG Signal Processing and Feature Extraction" (Springer)
%                                          
% Author:       ZHANG Zhiguo, zgzhang@szu.edu.cn
%               School of Biomedical Engineering, Shenzhen University, 
%               Shenzhen, China
%               Jan 2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and define parametersclear all
clear all; clc;

load data_lep
% data_lep.mat contains 3 variables
%   - x: the LEP trials (512 time points and 74 trials)
%   - Fs: the sampling rate (Fs = 256Hz)
%   - t: time index of the LEP trials (from -1 to 1 sec)
%   The description of the LEP experiment and data can be found in:
%   Zhang, Hu, et al. "Gamma-Band Oscillations in the Primary Somatosensory
%   Cortex¡ªA Direct and Obligatory Correlate of Subjective Pain Intensity",
%   Journal of Neuroscience, vol. 32, no. 22, pp. 7429-7438, May 2012.

[N_T, N_Trials]= size(x); % N_T: #time points = 512, N_Trials: #trials = 74 

%% STFT
nfft = 2^nextpow2(N_T); % the number of FFT points
winparam = round(Fs*0.4); % 400ms windows are used in STFT
[P, f, S] = subfunc_stft(x, winparam, nfft, Fs); % STFT
P_AVE = mean(P,3); % across-trial average of power values
F = S./sqrt(S.*conj(S)); % complex phase values
PLV = abs(mean(F,3));  % phase locking value (PLV)
% estimate the TFD of across-trial averaged LEP waveform for comparison 
P_LEP = subfunc_stft(mean(x,2), winparam, nfft, Fs);

%% Baseline correction
t_baseln_lim = [-0.5, -0.2]; % define the time interval of baseline correction (-0.5 to -0.2s in this example)
t_baseln_idx = find((t>=t_baseln_lim(1))&(t<=t_baseln_lim(2))); % the time indices of baseline samples

% baseline value is obtained from each trial, NOT all trials
P_Baseline_Mean_vec = mean(P(:,t_baseln_idx,:),2);
P_Baseline_Mean = repmat(P_Baseline_Mean_vec,[1,N_T,1]);
% P_lep_Baseline_Mean_Vec = mean(P_lep(:,t_baseln_idx),2);
% P_lep_Baseline_Mean = repmat(P_lep_Baseline_Mean_Vec,[1,N_T]);
% FOUR approaches for baseline correction 
% APPROACH 1: Substraction
    P_unit_1 = 'Power (\mu\rmV^2/Hz)';
    P_BC_1 = P - P_Baseline_Mean;
    P_BC_AVE_1 = mean(P_BC_1,3); % across-trial average of baseline-corrected TFD
% APPROACH 2: Relative Change (substract and divide)
    P_unit_2 = 'Percentage (100%)';
    P_BC_2 = (P - P_Baseline_Mean )./P_Baseline_Mean;
    P_BC_AVE_2 = mean(P_BC_2,3); % across-trial average of baseline-corrected TFD    
% APPROACH 3: Power Ratio (divide and log)
    P_unit_3 = 'Ratio (dB)';
    P_BC_3 = 20*log10(P./P_Baseline_Mean);
    P_BC_AVE_3 = mean(P_BC_3,3); % across-trial average of baseline-corrected TFD
% APPROACH 4: Z-score
    P_unit_4 = 'Z-score';
    P_Baseline_SD_vec = std(P(:,t_baseln_idx,:),[],2);
    P_Baseline_SD = repmat(P_Baseline_SD_vec,[1,N_T,1]);
    P_BC_4 = (P - P_Baseline_Mean )./P_Baseline_SD;
    P_BC_AVE_4 = mean(P_BC_4,3); % across-trial average of baseline-corrected TFD
    
%% Show results (LEP, TFD of LEP, averaged TFDs with/without basline correction, and PLV)
f_lim = [min(f(f>0)) 30]; % specify the frequency range to be shown (remove 0Hz)
f_idx = find((f<=f_lim(2))&(f>=f_lim(1)));
t_lim = [-0.5 1]; % specify the time range to be shown
t_idx = find((t<=t_lim(2))&(t>=t_lim(1)));

figure('units','normalized','position',[0   0.15    1    0.7])
% show LEP waveforms (all trials and the average)
subplot(2,3,1)
hold on; box on; 
plot(t(t_idx),x(t_idx,:),'linewidth',0.1,'color',[.5 .5 .5])
plot(t(t_idx),mean(x(t_idx,:),2),'linewidth',2,'color','k')
plot([0 0],[-50 50],'k:')
plot(t_lim,[0 0],'k:')
xlabel('Time (s)'); ylabel('Amplitude (\muV)')
set(gca,'ylim',[-50 50])
title('LEP Waveforms','fontsize',12)

% show TFD of averaged LEP
subplot(2,3,2)
imagesc(t(t_idx),f(f_idx),P_LEP(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,'\muV^2/Hz','rotation',90,'horizontalalignment','center','verticalalignment','top')
title('TFD of Averaged LEP','fontsize',12)
colorbar

% show averaged TFD of LEP (without baseline correction)
subplot(2,3,4)
imagesc(t(t_idx),f(f_idx),P_AVE(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,'\muV^2/Hz','rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Averaged TFD of LEP (w/o baseline correction)','fontsize',12)
colorbar

% show averaged TFD of LEP (with baseline correction/subtraction)
subplot(2,3,5)
imagesc(t(t_idx),f(f_idx),P_BC_AVE_1(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,'\muV^2/Hz','rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Averaged TFD of LEP (w/ baseline subtraction)','fontsize',12)
colorbar

% show PLV
subplot(2,3,6)
imagesc(t(t_idx),f(f_idx),PLV(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim,'clim',[0 1])
text(t_lim(2),f_lim(2)/2,'PLV','rotation',90,'horizontalalignment','center','verticalalignment','top')
title('PLV','fontsize',12)
colorbar

%% Show and compare four baseline-correction approaches for ERS/ERD estimation
f_lim = [min(f(f>0)), 30]; % specify the frequency range to be shown (remove 0Hz)
f_idx = find((f<=f_lim(2))&(f>=f_lim(1)));
t_lim = [-0.5, 1]; % specify the time range to be shown
t_idx = find((t<=t_lim(2))&(t>=t_lim(1)));

figure('units','normalized','position',[0.15   0.15    0.7    0.7])
% show baseline-correction APPROACH 1 (Subtraction)
subplot(2,2,1)
imagesc(t(t_idx),f(f_idx),P_BC_AVE_1(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,P_unit_1,'rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Baseline-corrected TFD (Subtraction)','fontsize',12)
colorbar

% show baseline-correction APPROACH 2 (Relative Change)
subplot(2,2,2)
imagesc(t(t_idx),f(f_idx),P_BC_AVE_2(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,P_unit_2,'rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Baseline-corrected TFD (Relative Change)','fontsize',12)
colorbar

% show baseline-correction APPROACH 3 (Power Ratio)
subplot(2,2,3)
imagesc(t(t_idx),f(f_idx),P_BC_AVE_3(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,P_unit_3,'rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Baseline-corrected TFD (Power Ratio)','fontsize',12)
colorbar

% show baseline-correction APPROACH 4 (Z-score)
subplot(2,2,4)
imagesc(t(t_idx),f(f_idx),P_BC_AVE_4(f_idx,t_idx))
axis xy; hold on;
plot([0 0],f_lim,'w--')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'xlim',t_lim,'ylim',f_lim)
text(t_lim(2),f_lim(2)/2,P_unit_4,'rotation',90,'horizontalalignment','center','verticalalignment','top')
title('Baseline-corrected TFD (Z-score)','fontsize',12)
colorbar
