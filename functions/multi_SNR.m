% multi_SNR: compute the SNR of a signal across trials
% this method will trend toward zero when there is no phase coherence across trials 
% and high SNR when local variance across trials is much higher than global variance
% 
% PARAMETERS
% ---------------
% signals: n x m array of signals.
%     n: samples in the signal
%     m: each trial in the block
% targetFrequency: frequncy of interest that bins width will be designated from
% fs: sampling frequency of the data.
% bins_per_cycle: number of bins to assign for one oscillation cycle
% 
% RETURNS
% ----------------
% SNR: signal to noise ratio of the signals
function [SNR,resFreq] = multi_SNR(signals, targetFrequency, fs, bins_per_cylce)
STDEV = std(signals(:)); % global STD of all the signals
optimalBinlen = ceil(fs/targetFrequency/bins_per_cylce);
nbins_init = floor(size(signals,1)/optimalBinlen);
nbins_h = nbins_init;
while mod(size(signals,1),nbins_h) ~= 0
    nbins_h = nbins_h +1;
end
nbins_l = nbins_init;
while mod(size(signals,1),nbins_l) ~= 0
    nbins_l = nbins_l -1;
end
if nbins_init - nbins_l < nbins_h - nbins_init
    nbins = nbins_l;
else
    nbins = nbins_h;
end
% fprintf('%d bins, target bins: %d\n',nbins,nbins_init)
resFreq = 1 / ((size(signals,1) / nbins) / fs ) / bins_per_cylce;
binRows = nbins;
binCols = size(signals,2);
sig = reshape(signals.', binCols,size(signals,1)/binRows, nbins);
sig = permute(sig,[2,1,3]);
% reshape(sig(:,1,:),[],1) == signals(:,1) 
% tt = reshape(sig(:,1,:),[],1);
% figure; plot(mean(signals,2)),hold on, plot(tt)
binSTD = squeeze(std(sig,[],[1,2]));
SNR_bin = STDEV./binSTD;
SNR = mean(SNR_bin);
end