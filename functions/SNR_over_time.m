%% SNR_over_time
% takes an input signal array of trials and computes SNR based on local vs
%   global variance. Breaks up signal based on intervals to track temporal changes
%   in SNR
%
% PARAMETERS
% ------------
% signals (array): m x n array of timeseries data. m is samples, n is
%   trials
% fs (int): sampling rate of the system
% intervals (int): intervals to split the original signals into for independent SNR calculation
% 
% OPTIONAL PARAMETERS
% binWidth (int): sets the bin width of the SNR calculation. Will supercede 
%       targetFreq. Default is 10 samples/bin.
% targetFreq (double): finds the minimum bin width to represent the
%       frequency of interest with 4 bins.
% RETURNS 
% -----------
% t (array): each time point specified in intervals
% SNR (array): SNR at each time point specified in intervals
function [t,SNR] = SNR_over_time(signals,fs,intervals,varargin)
t_all = linspace(0,size(signals,1)/fs,size(signals,1));

t = zeros(size(intervals));
SNR = zeros(length(intervals),1);
for i=1:length(intervals)
    onset = intervals(i,1);
    offset = intervals(i,2);
    sig_slice = signals(onset:offset,:);
    [avgSNR, ~] = timeseries_SNR(sig_slice,fs,varargin);
    t(i,:) = [t_all(onset),t_all(offset)];
    SNR(i) = avgSNR;
end

end