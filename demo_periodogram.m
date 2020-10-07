%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to estimate the spectrum of an
%               EEG signal using the periodogram.
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

%% spectral estimation (periodogram)
nfft = 2^nextpow2(N); % the number of FFT points
% check the help file to learn how to specify parameters in "peridogram.m"
[P_per, f] = periodogram(x,[],nfft,Fs); 

%% display spectral estimates
f_lim = f((f>0)&(f<=50)); % specify the frequency range to be shown

figure('units','normalized','position',[0.1    0.3    0.8    0.5])
subplot(1,2,1) 
hold on; box on;
plot(f,P_per,'k','linewidth',1) % show the periodogram in a linear scale
xlabel('Frequency (Hz)'); ylabel('Power (\muV^2/Hz)')
title(['Periodogram (in a linear scale)'],'fontsize',12)
set(gca,'xlim',[min(f_lim),max(f_lim)])

subplot(1,2,2)
hold on; box on; 
plot(f,10*log10(P_per),'k','linewidth',1) % show the periodogram in a log scale
xlabel('Frequency (Hz)'); ylabel('Power (dB)')
title(['Periodogram (in a logarithmic scale)'],'fontsize',12)
set(gca,'xlim',[min(f_lim),max(f_lim)])
