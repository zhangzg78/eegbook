%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to estimate the spectrum of an
%               EEG signal based on the AR model and using the Yule-Walker method.
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
% three parameter settings (AR orders) are used below
p1 = 20;
p2 = 10;
p3 = 50;
[P_ar_1,f] = pyulear(x,p1,nfft,Fs);
[P_ar_2,f] = pyulear(x,p2,nfft,Fs);
[P_ar_3,f] = pyulear(x,p3,nfft,Fs);

%% display spectral estimates
f_lim = f((f>0)&(f<=50)); % specify the frequency range to be shown

figure('units','normalized','position',[0.1    0.3    0.8    0.5])
subplot(1,2,1)
hold on; box on;
plot(f,10*log10(P_per),'k','linewidth',0.5)
plot(f,10*log10(P_ar_1),'r','linewidth',2)
xlabel('Frequency (Hz)'); ylabel('Power (dB)')
hl = legend('Periodogram','AR (P=20)');
set(hl,'box','off','location','southwest')
set(gca,'xlim',[min(f_lim),max(f_lim)])

subplot(1,2,2)
hold on; box on;
plot(f,10*log10(P_ar_1),'r','linewidth',2)
plot(f,10*log10(P_ar_2),'g','linewidth',1)
plot(f,10*log10(P_ar_3),'b','linewidth',1)
xlabel('Frequency (Hz)'); ylabel('Power (dB)')
hl = legend('AR (P=20)','AR (P=10)','AR (P=50)');
set(hl,'box','off','location','southwest')
set(gca,'xlim',[min(f_lim),max(f_lim)])
