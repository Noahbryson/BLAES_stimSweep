% CCEP + DTI general projections from human connectome -> probabilistic
% model of projection from basolateral complex to anterior/posterior
% hippocampus. wont have this for all subjects but can probe this in some. 
% 
% why don't we have access to de-identified EHR/Neuropsych evaluations?
% 
close all
carrier = 8;%hz
inner = 32;%hz
num_inner = 4;
carrier_t = 1/carrier * 1000; %ms
inner_t = 1/inner * 1000; %ms
inner_t_total = inner_t*num_inner;
f2 = 80; %Hz
f1 = 8;
fs = 2000;
t = linspace(0,1,fs);
x = square(f1*2*pi*t);
x = (x+1) / 2;
y = modulate_square(x,f1,f2,fs);
figure(1);
% y = modulate(x,f2,fs);
% y = (y+1)/2;
% locs = find(x==0);
% y(locs) = 0;
% plot(x)
hold on
plot(t,y)
hold off
ylim([-0.5,1.5])
xlim([-0.01,1])

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





