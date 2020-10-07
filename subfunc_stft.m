function [P, f, S] = subfunc_stft(x, winparam, nfft, fs)
% Short-time Fourier Transform (STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  This is a sub-function to calculate STFT.
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
% x:        the original data samples (Time Points x Channels/Trials)
% winparam: for an positive interger input, winparam is the window size (default window is hamming)
%           for a vector input, winparam is a window
% nfft:     number of fft points
% fs:       sampling rate

% // Output // %
% P:        squared magnitude of STFT (spectrogram)
% f:        evaluated frequency bins in STFT
% S:        complex time-frequency value of STFT

fprintf('\nShort-time Fourier Transform: ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify Parameters
if size(x,1)==1; x = x.'; end; % transpose data if the 1st dimension is 1
f = fs/2*linspace(0,1,nfft/2+1)';
N_Trials = size(x,2); % number of trials
N_T = size(x,1); % number of time samples
N_F = length(f); % number of frequency bins

fprintf('%d Time Points x %d Frequency Bins x %d Trial(s)\n',N_T,N_F,N_Trials);
S = single(zeros(N_F,N_T,N_Trials));
fprintf('Processing...     ')

%% Windowing and Padding
if length(winparam)==1 % a window size is specified
    if mod(winparam,2); h = winparam; %  window size (points)
    else h = winparam+1; %  window size (points); 
    end
    win = window('hamming',h); % window (one trial); default window type is hamming
else
    win = winparam;
    h = length(win);
end
W = repmat(win,1,N_Trials); % window (all trials)
U = win'*win;  % compensates for the power of the window

% Zero padding (default mode is "zero")
X = padarray(x,(h-1)/2);        % padding data
X = detrend(X,'linear');        % remove low-frequency trend 

%% STFT
for n=1:N_T
    fprintf('\b\b\b\b%3.0f%%',n/N_T*100)
    X_n = X(n+[0:h-1],:);   % windowed data segment
    X_n = detrend(X_n);     % remove low-frequency trend 
    S_n = fft(X_n.*W,nfft,1);   % FFT of windowed data
    S(:,n,:) = S_n(1:(nfft/2+1),:) / sqrt(U); % complex values
end
P = S.*conj(S)/fs;  % power values
fprintf('  Done!\n')

end