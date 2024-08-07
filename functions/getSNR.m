function [snr,p] = getSNR(signals,target_f,fs,bins_per_cycle)
% get single trial SNR distributions and averages at a target frequency
% target frequency determines the size of the bins, as bins >= 1/2 the
% oscillation time will have the same variance as the entire signal itself.
%
% signals: array of n_trials x n_datapoints
% target_f: frequency of interest, determines size of timeseries bins
% fs: sampling frequency of signals
snr = zeros(size(signals,1),1); 
siglen = size(signals,2); % in samples
cycle_len = fs/target_f;
binlen = cycle_len/bins_per_cycle; % in samples
n_bins = floor(siglen / binlen);
maxlen = siglen - mod(siglen,n_bins);
binlen = maxlen/n_bins;
sliding_adjust = uint16((siglen - maxlen)/2);
signals = signals(:,1+sliding_adjust:siglen-sliding_adjust);
t = linspace(0,size(signals,2)/fs,size(signals,2));
t = repmat(t,size(signals,1),1);
for i =1:size(signals,1)
    tot_var = var(signals(i,:));
    loc_var = zeros(n_bins,1);
    data = reshape(signals(i,:),[],n_bins);
    % t = reshape(t(i,:),[],n_bins);
    % plot_bins(t,data)
    for bin=1:n_bins
        loc_var(bin) = var(data(:,bin));
    end
    snr(i) = tot_var / mean(loc_var);
end

snr_avg = mean(snr);
p="NONE";

end


function plot_bins(t,y)
figure(1)
hold on
for i=1:size(t,2)
    plot(t(:,i),y(:,i));
    xline(t(end,i))
end
end