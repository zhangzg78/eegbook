%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to calculate the short-time Fourier 
%               transform (STFT) of a VEP signal.
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

load data_vep.mat
% data_vep.mat contains 3 variables
%   - x: the VEP signal (512 time points, averaged from multiple trials)
%   - Fs: the sampling rate (Fs = 250Hz)
%   - t: the time indices (256 pre-stimulus samples and 256 post-stimulus samples, with a time interval of 1/Fs)
%   The original data are from https://vis.caltech.edu/~rodri/data/cg_o1t.asc
%   Please see https://vis.caltech.edu/~rodri/data.htm for more details.

N = numel(x); % the number of time points
x = x - mean(x(t<0)); % baseline corrction for averaged VEP

%% perform time-frequency analysis using STFT
nfft = 2^nextpow2(N); % the number of FFT points
winsize = round([0.2 0.4 0.8]*Fs); % three window sizes (0.2s, 0.4s, 08s) are used in STFT for comparison
for n=1:numel(winsize)
    [P(:,:,n),f] = subfunc_stft(x, winsize(n), nfft, Fs);
end

%% display VEP and STFT results with different window sizes
f_lim = [min(f(f>0)), 30]; % specify the frequency range to be shown (remove 0Hz)
f_idx = find((f<=f_lim(2))&(f>=f_lim(1)));
t_lim = [-0.2, 1]; % specify the time range to be shown
t_idx = find((t<=t_lim(2))&(t>=t_lim(1)));

figure('units','normalized','position',[0.1    0.15    0.8    0.7])
subplot(2,2,1)
hold on; box on
plot(t(t_idx),x(t_idx),'k','linewidth',1);
plot([-1 1],[0 0],'k--')
plot([0 0],[-30 30],'k--')
set(gca,'xlim',[min(t_lim),max(t_lim)])
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title(['VEP'],'fontsize',12)
for n=1:numel(winsize)
    subplot(2,2,n+1)
    imagesc(t(t_idx),f(f_idx),P(f_idx,t_idx,n))
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(gca,'xlim',t_lim,'ylim',f_lim)
    axis xy; hold on;
    plot([0 0],[0 Fs/2],'w--')
    text(t_lim(2),f_lim(2)/2,'Power (dB)','rotation',90,'horizontalalignment','center','verticalalignment','top')
    title(['Spectrogram (winsize = ',num2str(winsize(n)/Fs,'%1.2g'),'s)'],'fontsize',12)
    colorbar
end
