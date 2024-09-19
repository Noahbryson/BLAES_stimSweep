% stim_predicition: locate ideal pulse locations during a stim train for
% interpolation
%  
% PARAMETERS
% ---------------
% signal: trigger signal during one stimulation trial.
% fs: sampling frequency of the data.
% fc_b: base frequency of indiviudal pulses
% fc_m: modulation frequency of the carrier wave, typically lower than fc_b
% n_peaks: number of stimulation pulses to assign.
% n_cycles: number of total repeats due to modulation frequency
% 
% RETURNS
% ----------------
% spike_idx: logical array of samples in the signal containing a
% theoretical spike based on stim parameters. 
function spike_idx = stim_prediction(signal, fs, fc_b, fc_m, n_peaks, duration, onsetSample)
onsetSample = floor(onsetSample);
n_cycles = duration * fc_m;
oscil_samps = floor(fs / fc_m);
half_oscil = floor(oscil_samps/2);
spike_dist = floor(fs / fc_b);
train_idx = zeros(length(signal),1);
spike_idx = zeros(length(signal),1);
loc = onsetSample;
for i=1:n_cycles
train_idx(loc:loc+half_oscil) = 1;
spike_idx(loc:spike_dist:loc+half_oscil) = 1;
loc = loc + oscil_samps;
% loc = loc + half_oscil;
% train_idx(loc:loc+half_oscil) = 0;
% loc = loc + half_oscil;
end

% spike_idx = (spike_idx -1)*1.5;
% plot(abs(signal))
% hold on
% xline(linspace(1,fc_m,fc_m)*fs/fc_m,'Color',[0 1 0])
% plot(train_idx * max(signal),'Color',[1 0 0]);
% plot(spike_idx * 500,'Color',[1 1 1],'LineWidth',0.5, 'LineStyle','-');
end