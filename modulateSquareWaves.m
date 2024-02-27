function [t,y,y1,y2] = modulateSquareWaves(duration, freq, carrier, fs, amp,plotFlag)

t = linspace(0,duration,fs*duration); %tspan in seconds
y1 = square(freq*2*pi*t);
y1 = (y1 - min(y1));
y1 = y1/max(y1)*amp;
y2 = square(carrier*2*pi*t);
y2 = (y2 - min(y2));
y2 = y2/max(y2)*amp;
y=zeros(size(y1));
for i=1:numel(y1)
    if y2(i) == amp
        y(i) = y1(i);
    end
end
if plotFlag
tscale = 0.05;
figure(1)
subplot(3,1,1);
plot(t,y1)
xlim([-tscale duration+duration*tscale])
ylim([-0.1 1.1])
subplot(3,1,2);
plot(t,y2)
xlim([-tscale duration+duration*tscale])
ylim([-0.1 1.1])
subplot(3,1,3);
plot(t,y)
xlim([-tscale duration+duration*tscale])
ylim([-0.1 1.1])
end
end