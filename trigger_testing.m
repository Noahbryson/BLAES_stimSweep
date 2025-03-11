%% Trigger Testing
close all
load("dat/triggertesting.mat")
channel = 1;
cname = dat_channelNames{1};
t = linspace(0,1,length(xxx));
xx_l = downsample(xxx,4);
t_l = linspace(0,1,length(xx_l));
trigger_l = downsample(trigger,4);
f=figure;
f.Position = [100 100 1540 1400];
subplot(3,1,1)
plot(t,xxx)
hold on
plot(t_l,xx_l)
hold off
subplot(3,1,2)
plot(t,abs(trigger))
hold on
plot(t_l,abs(trigger_l))
hold off
subplot(3,1,3)
trig_x = abs(trigger);
trig_x = trig_x/max(trig_x);
[peaks,peaklocs] = findpeaks(trig_x,"MinPeakHeight",4*std(trig_x),'MinPeakDistance',20);
% [peaks,peaklocs] = findpeaks(trig_x,"MinPeakHeight",0.7);
plot(t,trig_x)
hold on 
scatter(t(peaklocs),peaks)
