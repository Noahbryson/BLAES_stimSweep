function [snr, p] = oscillatory_SNR_tao(signals, fs, windowsize, n_permute)
% signals: a 2D array of trials x samples.
% fs: sampling frequency
% windowsize: duration of windows in msec to split the signals into.
% n_permute: iterations for permutation testing.
siglen = size(signals,2);
windowlen = (windowsize * fs / 1000);
n_bins = floor(siglen / windowlen);
perm_snr = zeros(n_permute);
p =  zeros(n_permute);
maxlen = siglen - mod(siglen,n_bins);
sliding_adjust = uint16((siglen - maxlen)/2);
signals = signals(:,1+sliding_adjust/2:siglen-sliding_adjust/2);
snr = single_chan_snr(signals,windowlen,n_bins);
% for pm=1:n_permute
%     perm_snr(pm) = single_chan_snr(signals,windowlen,n_bins);
% end
end


function snr = single_chan_snr(signal,bin_size,n_bins)
data = reshape(signal.',bin_size,n_bins,size(signal,1));
data = permute(data, [3, 1, 2]); % trials x data point x bin #
var_flat1 = (var(signal));
var_flat = reshape(var_flat1.',bin_size,n_bins);
var_bin = zeros(n_bins,1);
var_flat2bin = zeros(n_bins,1);


for bin=1:n_bins  
    t = mean(var(data(:,:,bin)));
    var_bin(bin) = t;
    var_flat2bin(bin) = mean(var_flat(:,bin));
end
var_total = mean(var(signal));
snr = var_total./var_bin;

end
