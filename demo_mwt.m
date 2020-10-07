%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is the demo script to calculate the continuous wavelet
%               transform (CWT) of a VEP signal with the Morlet wavelet basis.
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

%% perform time-frequency analysis using MWT
nfft = 2^nextpow2(N); % the number of FFT points
[P_stft,f] = subfunc_stft(x,100, nfft, Fs); % STFT (with a window size of 100) is calculated for comparison
% note: the frequency index f is obtained from STFT in this demo, while it 
% can actually be defined by users.
N_F = length(f); % number of frequency bins

omega = [0.5 1];    % test two values for the parameter omega
sigma = [0.5 .2];   % test two values for the parameter sigma

ff = f/Fs;     % normalized frequency
L_hw = N; % filter length
for fi=1:N_F
    scaling_factor = omega(1)./ff(fi);
    u = (-[-L_hw:L_hw])./scaling_factor;
    w(:,fi) = sqrt(1/scaling_factor)*exp(-(u.^2)/(2*sigma(1).^2));
    hw(:,fi) = w(:,fi).'.* exp(1i*2*pi*omega(1)*u);
end

for n_omega=1:numel(omega)
    for n_sigma=1:numel(sigma)
        [P_mwt(:,:,n_omega,n_sigma)] = subfunc_mwt(x, f, Fs, omega(n_omega), sigma(n_sigma));
    end
end

%% display MWT results with different parameters
f_lim = [min(f(f>0)), 30]; % specify the frequency range to be shown (remove 0Hz)
f_idx = find((f<=f_lim(2))&(f>=f_lim(1)));
t_lim = [-0.2, 1]; % specify the time range to be shown
t_idx = find((t<=t_lim(2))&(t>=t_lim(1)));

figure('units','normalized','position',[0.1    0.15    0.8    0.7])
for n_omega=1:numel(omega)
    for n_sigma=1:numel(sigma)
        n = (n_omega-1)*2+n_sigma;
        subplot(2,2,n)
        imagesc(t(t_idx),f(f_idx),P_mwt(f_idx,t_idx,n_omega,n_sigma))
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        set(gca,'xlim',t_lim,'ylim',f_lim)
        axis xy; hold on;
        plot([0 0],[0 Fs/2],'w--')
        title(['Scaleogram (\omega = ',num2str(omega(n_omega),'%1.2g'),', \sigma = ',num2str(sigma(n_sigma),'%1.2g'),')'],'fontsize',12)
        text(t_lim(2),f_lim(2)/2,'Power (\muV^2/Hz)','rotation',90,'horizontalalignment','center','verticalalignment','top')
        colorbar 
    end
end
