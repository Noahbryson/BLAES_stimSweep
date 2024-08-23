function [res,corr,lags] = res_autocorr(signals,fs,test_freqs,varargin)
if ~isempty(varargin)
    maxlag = varargin{1};
else
    maxlag = fs;
end

res = struct;
[corr,lags] = xcorr(signals,maxlag,'normalized');
corr = abs(corr);
center_idx = ceil(length(corr)/2);
sampleIDX = linspace(1,length(corr),length(corr));
neg_loc = sampleIDX(1:center_idx-1); 
pos_loc = flip(sampleIDX(center_idx+1:end)); 

for i=1:length(test_freqs)
% sample_dist = ceil(fs/test_freqs(i) * (maxlag/fs));
sample_dist = ceil(fs/test_freqs(i));

buffer = ceil(0.05*sample_dist* (maxlag/fs));
padded_corr = [nan(1, buffer) corr nan(1, buffer)];
padded_lag = [nan(1, buffer) lags nan(1, buffer)];
neg_loc1 = neg_loc(1:sample_dist:end)+buffer;
negIdx = bsxfun(@plus, neg_loc1', -buffer:buffer);
negSig = padded_corr(negIdx);
negLags = padded_lag(negIdx);

[negCorr, negMaxIdx] = max(negSig,[],2);
[row, col] = size(negLags);
rowIdx = (1:row)';
linIdx = sub2ind([row col], rowIdx, negMaxIdx);
negLag = negLags(linIdx);
pos_loc1 = flip(pos_loc(1:sample_dist:end));
posIdx = bsxfun(@plus, pos_loc1', -buffer:buffer);
posSig = padded_corr(posIdx);
posLags = padded_lag(posIdx);
[posCorr, posMaxIdx] = max(posSig, [],2);
[row, col] = size(posLags);
rowIdx = (1:row)';
linIdx = sub2ind([row col], rowIdx, posMaxIdx);
posLag = posLags(linIdx);
% temp = neg_loc1 < buffer+1; %##TODO: grab window of size buffer around samples and take max corr in window as value
% neg_loc1(temp) = neg_loc1(temp)+buffer;
% temp = pos_loc1 > length(corr) - buffer;
% neg_loc1(temp) = pos_loc1(temp)-buffer;

rescorr= [negCorr; posCorr];
resLag = [negLag' posLag'];
rmIdx = resLag == 0;
rescorr(rmIdx) = [];
resLag(rmIdx) = [];
res(i).freq = test_freqs(i);
res(i).lags = resLag;
res(i).corr = rescorr;
res(i).avg_corr = mean(rescorr);
res(i).max_corr = max(rescorr);
% close all
% figure(1)
% plot(lags/fs,corr)
% hold on
% scatter(resLag/fs,rescorr)
% title(sprintf('%d Hz', test_freqs(i)))
end
end

% neg_lags = flip(lags(1:center_idx-1),1); 
% neg_lags = neg_lags(1:sample_dist:end);
% pos_lags = flip(lags(center_idx+1:end)); 
% pos_lags = flip(pos_lags(1:sample_dist:end));
% neg_corr = flip(corr(1:center_idx-1),1); 
% neg_corr = neg_corr(1:sample_dist:end)';
% pos_corr = flip(corr(center_idx+1:end)); 
% pos_corr = flip(pos_corr(1:sample_dist:end))';
% reslag = [neg_lags pos_lags];
% rescorr = [neg_corr; pos_corr];