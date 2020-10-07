%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to estimate the spectrum of an
%               EEG signal using the multitaper method.
%
% Note:         This script is a supplementary file of Chapter6 in the book 
%               "EEG Signal Processing and Feature Extraction" (Springer)
%                                          
% Author:       ZHANG Zhiguo, zgzhang@szu.edu.cn
%               School of Biomedical Engineering, Shenzhen University, 
%               Shenzhen, China
%               Jan 2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and define parameters
clear all; clc;

load data_eeg.mat
% data_eeg.mat contains 2 variables
%   - x: the EEG signal (with eyes-closed)
%   - Fs: the sampling rate (Fs = 160Hz)

N = length(x); % the number of samples (N=480)
x = detrend(x); % remove the low-frequency trend from EEG

%% spectral estimation (the multitaper method)
nfft = 2^nextpow2(N); % the number of FFT points
[P_per, f] = periodogram(x,[],nfft,Fs); % periodogram is also estimated for comparison
% check the help file to learn how to specify parameters in "pmtm.m"
% three parameter settings are used below
[P_mtm_1, f] = pmtm(x,4,nfft,Fs);
[P_mtm_2, f] = pmtm(x,2,nfft,Fs);
[P_mtm_3, f] = pmtm(x,6,nfft,Fs);

%% display spectral estimates
f_lim = f((f>0)&(f<=50)); % specify the frequency range to be shown

figure('units','normalized','position',[0.1    0.3    0.8    0.5])
subplot(1,2,1)
hold on; box on;
plot(f,10*log10(P_per),'k','linewidth',0.5)
plot(f,10*log10(P_mtm_1),'r','linewidth',2)
xlabel('Frequency (Hz)'); ylabel('Power (dB)')
hl = legend('Periodogram','Multitaper (L=4)');
set(hl,'box','off','location','southwest')
set(gca,'xlim',[min(f_lim),max(f_lim)])

subplot(1,2,2)
hold on; box on;
plot(f,10*log10(P_mtm_1),'r','linewidth',2)
plot(f,10*log10(P_mtm_2),'g','linewidth',1)
plot(f,10*log10(P_mtm_3),'b','linewidth',1)
xlabel('Frequency (Hz)'); ylabel('Power (dB)')
hl = legend('Multitaper (L=4)','Multitaper (L=2)','Multitaper (L=6)');
set(hl,'box','off','location','southwest')
set(gca,'xlim',[min(f_lim),max(f_lim)])
