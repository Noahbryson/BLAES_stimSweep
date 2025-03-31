%% Entrainment_Methodology_Examples
fs = 2000;
time = 1.5;
stimTime=1;
freq = 8;
n = 50;
noise_intensity = 0.5; % signal amplitude is +/- 1, so noise intensity of 2 makes SNR=1
t = linspace(0,time,time*fs);
%% Representative Signal Example
figure(1)
testsig=repmat(sin(t*2*pi*freq),10,1) + 0.1*randn(10,3000)-0.5 ;
plot(t,testsig)

%% Entrainment Magnitude


%% Entrainment Phase
close all


postStimIDX = stimTime*fs+1;

dampening = flip(exp(-0.1*t));
dampening(1:postStimIDX) = dampening(1:postStimIDX).*t(1:postStimIDX);
postStimt=linspace(0,time-stimTime,fs*(time-stimTime));
postStimDamp = exp(-5*postStimt)';
dampening(postStimIDX:end) = dampening(postStimIDX-1)*postStimDamp;
signal1 = sin(t*2*pi*freq).*dampening;
signal2 = sin(t*2*pi*freq-25/(2*pi)).*dampening;
signal3 = sin(t*2*pi*freq-42/(2*pi)).*dampening;

[locs,stim] = stimtrain(signal1,time,stimTime,fs,freq,0.25,4);
stims = repmat(stim,1,n)';
y1 = repmat(signal1,n,1);
y2 = repmat(signal2,n,1);
y3 = repmat(signal3,n,1);
noise = noise_intensity*(rand(n,time*fs)-0.5).*flip(exp(-1.1*t));
yy1 = noise+y1+stims;
yy2 = noise+y2+stims;
yy3 = noise+y3+stims;
% plot(t,signal1,t,signal2,t,signal3,t,signal3+stim')

interval = floor(2 * fs/freq)-fs*0.1;
% interval = 490;
[t1, y1s] = entrainment_phase(yy1',locs(5:end),interval,fs);

[t2, y2s] = entrainment_phase(yy2',locs(5:end),interval,fs);

[t3, y3s] = entrainment_phase(yy3',locs(5:end),interval,fs);

figure(1)
% subplot(1,3,1)
a1=subplot(3,1,1);
plot(t1,y1s','Color',[0.4 0.4 0.4])
hold on
plot(t1,mean(y1s),'Color',[0.6 0 0.4])
ylabel('amplitude (au)')
title('Channel A')
% subplot(1,3,2)
a2=subplot(3,1,2);
plot(t2,y2s','Color',[0.4 0.4 0.4])
hold on
plot(t2,mean(y2s),'Color',[0.6 0 0.4])
ylabel('amplitude (au)')
title('Channel B')
% subplot(1,3,3)
a3=subplot(3,1,3);
plot(t3,y3s','Color',[0.4 0.4 0.4])
hold on
plot(t3,mean(y3s),'Color',[0.6 0 0.4])
xlabel('time (s)')
ylabel('amplitude (au)')
title('Channel C')


figure(4)
a21=subplot(3,1,1);

plot(t,(yy1+stims),'Color',[0.5,0.5,0.5]),hold on, plot(t,mean(yy1+stims))
ylabel('amplitude (au)')
title('Channel A')
xline(t(locs(5:end)-interval/2),'LineWidth',2.5,'Color',[0 1 0])
xline(t(locs(5:end)+interval/2),'LineWidth',2.5,'Color',[1 0 0])

a22=subplot(3,1,2);

plot(t,(yy2+stims),'Color',[0.5,0.5,0.5]),hold on, plot(t,mean(yy2+stims))
ylabel('amplitude (au)')
title('Channel B')
xline(t(locs(5:end)-interval/2),'LineWidth',2.5,'Color',[0 1 0])
xline(t(locs(5:end)+interval/2),'LineWidth',2.5,'Color',[1 0 0])

a23=subplot(3,1,3);

plot(t,(yy3+stims),'Color',[0.5,0.5,0.5]),hold on, plot(t,mean(yy3+stims))
xlabel('time (ms)')
ylabel('amplitude (au)')
xline(t(locs(5:end)-interval/2),'LineWidth',2.5,'Color',[0 1 0])
xline(t(locs(5:end)+interval/2),'LineWidth',2.5,'Color',[1 0 0])

title('Channel C')
linkaxes([a21,a22,a23],'x')
%% Post-Stim Feature


%% Entrainment Onset

%% Functions
function [idx,y] = stimtrain(X,sampTime,stimTime,fs,stimfreq,pulseWidth,pulsePerTrain)
pw = pulseWidth/1000*fs; %pulseWidth is in ms
samplen = fs*sampTime;% time in s, fs in hz
y = zeros(samplen,1);
samp_per_cycle = fs/stimfreq;
offset = floor(0.5*samp_per_cycle);
n_pulse = stimfreq*stimTime;
idx = zeros(n_pulse,1);
for i = 1:n_pulse
idx(i) = offset + (i-1)*floor(samp_per_cycle);
line = linspace(1.5*max(X),2*min(X),3);
lineOutput = repmat(line,1,pulsePerTrain);
y(idx(i):idx(i)+length(lineOutput)-1) = lineOutput;
end
end