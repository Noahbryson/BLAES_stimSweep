close all
time = 1; fs = 1000;
t = linspace(0,time,fs*time);
targetFreq = 40;
y1 = sin(t*2*pi*targetFreq);
y1 = circshift(y1,ceil(fs/targetFreq/2));
y0=y1;
freqz = [15 33 44 50 21 60 80 66, 3, 18, 38, 90];
for i=1:length(freqz)
    temp = sin(t*2*pi*freqz(i)) * 0.9*rand(1);
    y1 = y1 + temp;
end
% y1 = y1 + sin(t*2*pi*60);
y1 = y1 / max(abs(y1));
f = linspace(2,150,149);
numLags = 500;

%% Cross Correlation
close all
% [maxVal, maxLoc] = max([theta_epochs.snr_av]);
maxLoc = 201;
targetChan = theta_epochs(maxLoc);
name = targetChan.channel;
yyy = targetChan.baseline';
xxx = targetChan.signals';
res = res_crosscor(yyy,fs,f,numLags);
agg_res = aggregate_corr_trials(res);
res2 = res_crosscor(xxx,fs,f,numLags);
agg_res2 = aggregate_corr_trials(res2);
nrows = 3; ncols = 2;
figure(1)
subplot(nrows,ncols,1)
plot(yyy)
title(sprintf('%s Baseline',name))
subplot(nrows,ncols,2)
plot(xxx)
title(sprintf('%s Stim',name))
subplot(nrows,ncols,3)
plot([res.freq],[res.avg_corr])
subplot(nrows,ncols,4)
plot([res2.freq],[res2.avg_corr])
subplot(nrows,ncols,5)
plot([res.freq],agg_res)
subplot(nrows,ncols,6)
plot([res2.freq],agg_res2)
% subplot(nrows,ncols,7)
% plot(yyy)
% subplot(nrows,ncols,8)
% plot(xxx)

%% Autocorrelation
[res,corr,lags] = res_autocorr(y0,fs,f,numLags);
[res2,corr2,lags2] = res_autocorr(y1,fs,f,numLags);
targetIdx = [res.freq] == targetFreq;
figure(2)
subplot(4,2,1)
plot(t,y0,'LineWidth',2)
subplot(4,2,2)
plot(t,y1,'LineWidth',2)
subplot(4,2,3)
hold on, xline(freqz)
hold on, xline(targetFreq,'LineWidth',4)
plot([res.freq],[res.avg_corr],'LineWidth',2)
subplot(4,2,4)
hold on, xline(freqz)
hold on, xline(targetFreq,'LineWidth',4)
plot([res2.freq],[res2.avg_corr],'LineWidth',2)
subplot(4,2,5)
% stem(lags,corr)
scatter(res(targetIdx).lags/fs,res(targetIdx).corr)
hold on
plot(lags/fs,corr)
subplot(4,2,6)
% stem(lags2,corr2)
scatter(res2(targetIdx).lags/fs,res2(targetIdx).corr)
hold on
plot(lags2/fs,corr2)

window = hamming(length(y0)*0.5);
subplot(4,2,7)
[Pxx, f] = pwelch(y0,window,[],[],fs);
semilogy(f,Pxx)
subplot(4,2,8)
[Pxx, f] = pwelch(y1,window,[],[],fs);
semilogy(f,Pxx)

%%
function X = aggregate_corr_trials(res_struct)
X = zeros(length(res_struct(1).max_corr),length(res_struct));
for i=1:length(res_struct)
    X(:,i) = res_struct(i).max_corr;
end


end