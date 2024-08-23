function res = res_crosscor(signals,fs,test_freqs,varargin)
% runs cross correlation on a matrix with pure sinusoids to test for
% correlation and phase shift
%
% signals n x m matrix where m is a single channel and n are the samples of
% m. This function operates column-wise on signals
%
% fs: sampling freq of data
%
% test_freqs: frequencies to assess
%
% varagin: pass numerical arguments. 
%       Varargin(1): maximum lag
%       Varargin(2): 


if ~isempty(varargin)
    maxlag = varargin{1};
else
    maxlag = fs;
end
res = struct;
t = linspace(0,size(signals,1)/fs,size(signals,1));
for i=1:length(test_freqs)
b = sin(2*pi*t*test_freqs(i));
lagout = zeros(maxlag*2+1,size(signals,2));
corrout = lagout;
maxcorr_out = zeros(size(signals,2),1);
maxlag_out = maxcorr_out;
for j=1:size(signals,2)
[corr,lags] = xcorr(signals(:,j),b,maxlag,'normalized');
corrout(:,j) = corr;
lagout(:,j) = lags;
[maxcorr_out(j), idx] = max(corr);
maxlag_out(j) = lags(idx);
end
res(i).freq = test_freqs(i);
res(i).corrs = corrout;
res(i).lags = lagout/fs;
res(i).max_corr = maxcorr_out;
res(i).phaseshift = maxlag_out/fs;
res(i).avg_corr = mean(maxcorr_out);
res(i).avg_lag = mean(maxlag_out/fs);
end
end