%% Script for processing power spectral densities of a matrix of data
% signals: n x m array, where m is the number of signals and n is the
% number of samples in the signals
%
% frequencies: 1D array of frequencies to compute the PSD at.
%
% fs: sampling frequency of the data
%
% returns:
%   f: 1D array, frequencies that PSD is computed at
%   Pxx: f x m array, Power spectral densities of the m signals.

function [f,Pxx] = matrix_PSD(signals,frequecies,fs)

window = hamming(size(signals,1)*0.5);
[Pxx, f] = pwelch(signals,window,[],frequecies,fs);
plot = 0;
if plot
    figure(1);
    semilogy(f,mean(Pxx,2))
end
end