function [time,y] =generate_theta_burst_waveform(frequency,amplitude)

f2 = frequency; %Hz
f1 = 8;
fs = 2000;
time = linspace(0,1,fs);
x = square(f1*2*pi*time);
x = (x+1) / 2;
y = amplitude*modulate_square(x,f1,f2,fs);
y(1)=0;

function x = modulate_square(waveform,base_f,mod_f,fs)
locs = find(waveform==max(waveform));
intervals = zeros(base_f,2);
intervals(1,1) = 1;
intervals(end,end) = locs(end-1);
idx = 1;
for i=1:numel(locs)-2
    if locs(i)+1 ~= locs(i+1)
        intervals(idx,2) = locs(i); %offset
        idx = idx+1;
        intervals(idx,1) = locs(i+1); %onset
        
    end
end
x = zeros(size(waveform));
for i=1:size(intervals,1)
    t = linspace(0,(intervals(i,2)-intervals(i,1)+1)/fs,intervals(i,2)-intervals(i,1)+1);
    temp = square(mod_f*2*pi*t);
    temp = (temp-min(temp));
    temp = temp/max(temp);
    x(intervals(i,1):intervals(i,2)) = temp;
end
end
end
