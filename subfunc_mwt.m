function [P, S] = subfunc_mwt(x, f, Fs, omega, sigma)
% Morlet wavelet transform (continuous wavelet transform with Morlet basis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is a sub-function to calculate Morlet wavelet transform.
%
% Note:         This script is a supplementary file of Chapter6 in the book 
%               "EEG Signal Processing and Feature Extraction" (Springer)
%                                          
% Author:       ZHANG Zhiguo, zgzhang@szu.edu.cn
%               School of Biomedical Engineering, Shenzhen University, 
%               Shenzhen, China
%               Jan 2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% // Input // %
% x:            the original data samples (Time Points x Channels)
% Fs:           sampling rate
% omega:        a parameter to define the central frequency of Morlet wavelet
% sigma:        a parameter to define the spread of Morlet wavelet in time domain

% // Output // %
% P:            squared magnitude of MWT (scaleogram)
% S:            complex values of Morlet wavelet transform

fprintf('Calculating Morlet Wavelet Transform ... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pre-processing and Parameters
if size(x,2)==1; x = x.'; end
x = detrend(x,'linear'); % remove linear trends

N_F = length(f); % number of frequency bins
N_T = length(x); % number of time samples
f = f/Fs;     % normalized frequency

S = single(zeros(N_F,N_T)); % define the size of output 

%% Morlet wavelet transform
L_hw = N_T; % filter length
for fi=1:N_F
    scaling_factor = omega./f(fi); % the scaling factor
    u = (-[-L_hw:L_hw])./scaling_factor;
    hw = sqrt(1/scaling_factor)*exp(-(u.^2)/(2*sigma.^2)).* exp(1i*2*pi*omega*u);
    S_full = conv(x,conj(hw));
    S(fi,:) = S_full(L_hw+1:L_hw+N_T); % complex values
end
P = abs(S).^2;  % power values

fprintf('Done!\n')
