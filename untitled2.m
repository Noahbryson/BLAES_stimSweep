on = 104001;
off = 105500;
% on = 62001;
% off = 64000;
% on = 36401;
% off = 37900;
testSignal = signals(on:off,25);
time = linspace(0,length(testSignal)/fs,length(testSignal));
noise = sin(time*2*pi*60)*0.1*abs(max(testSignal));
figure(1)
subplot(2,1,2)
plot(time,testSignal+noise','LineWidth',2);
subplot(2,1,1)
plot(time,testSignal);

%% CAR
testSignals = signals(on:off,21:30) + repmat(noise,10,1)';
figure(2)
subplot(2,1,1)
hold on
for i=1:10
    plot(time,testSignals(:,i) + 0.2*max(testSignal) * i',"Color",'black');
end
plot(time,mean(testSignals,2)-2*max(testSignal),'Color','red','LineWidth',2)
hold off
subplot(2,1,2)
plot(time,testSignal+noise','LineWidth',1);
hold on
sigres = testSignal+noise'-mean(testSignals,2);
plot(time,sigres);
Fn = fs/2;
[b,a] = butter(4,2/Fn,'high');
filt_data = filtfilt(b,a,sigres);
plot(time,filt_data)
hold off
legend({'raw' 'reref' 'filt'})
saveas(gcf,fullfile(rawfilesavepath,'CAR_ref.svg'))
%% Lin-interp 
[b,a] = butter(4,80/Fn,'high');
trigger = abs(filtfilt(b,a,testSignal));
figure(3)
tcl = tiledlayout(2,1);
ax = nexttile(tcl);
[pk,pkl] = findpeaks(trigger,'MinPeakHeight',.5*std(trigger),'MinPeakDistance',floor(1/80*fs));
hold on
plot(ax,time,trigger)
scatter(ax,time(pkl),pk)
spikesig = zeros(length(time),1);
spikesig(pkl) = 1;
interpres = interpolateSpikes(filt_data,0,spikesig,fs,10,10);
hold off
ax = nexttile(tcl);
plot(ax,time,testSignal)
hold on
plot(ax,time,interpres)
hold off
linkaxes(tcl.Children,'x')
saveas(gcf,fullfile(rawfilesavepath,'lin_interp.svg'))
%% Entrainment Synthetic Figures
t = linspace(0,1,1000);
y = sin(2*pi*25*t);
dy = linspace(.8,0.2,900);
dy = circshift(dy,floor(800*rand(1)));
Y = y; Y(101:end) = Y(101:end) .* dy; 
figure(4);
tcl = tiledlayout(3,1);
ax =nexttile(tcl);
plot(t,y)
ax =nexttile(tcl);
plot(t,Y)
ax=nexttile(tcl);
plot(t,square(t*2*pi*25));
ylim([-4,4])
linkaxes(tcl.Children,'y')
saveas(gcf,fullfile(rawfilesavepath,'phase_synth.svg'))

