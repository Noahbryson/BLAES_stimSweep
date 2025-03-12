%% run_theta_entrainment
clear 
user = expanduser('~'); % Get local path for interoperability on different machines, function in my tools dir. 
if ispc
    boxpath = fullfile(user,"Box/Brunner Lab"); % Path to data
     BCI2KPath = "C:\BCI2000\BCI2000";
else
    boxpath =  fullfile(user,'Library/CloudStorage/Box-Box/Brunner Lab'); % Path to data
    BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
end
rootDataPath = fullfile(boxpath,"/DATA/BLAES/BLAES_param");
addpath(genpath(fullfile(user,'Documents/NCAN/code/BLAES_stimSweep')));
addpath(genpath(fullfile(user,'Documents/NCAN/code/MATLAB_tools')));
bci2ktools(BCI2KPath);
%%
all_subject_info = readtable(fullfile(rootDataPath,'Subject_Locations.xlsx'));
UtahSubs = {'UIC202407' 'UIC202412' 'UIC202414'};
BJHSubs = {'BJH050' 'BJH052' 'BJH056'};
allSubList = [BJHSubs UtahSubs];
% for subIdx = 1:length(allSubList)
Subject = allSubList{end};
[epochs,epoch_meta,epoch_corrs,pulseLocs]=load_data(fullfile(rootDataPath,Subject,'processed'));
fs = epochs.fs;
n_trials = 24;
%%
trainLen = size(epochs.signals,1);
train_intervals = get_train_intervals(pulseLocs,n_trials,fs,8,trainLen);

%% data viz
count = 1;
for i=111:112

loc = i;
reg = epoch_meta(loc).region;
stimInfo = epoch_meta(loc).label;
SNR = timeseries_SNR(epochs.signals(:,:,i),epochs.fs,'targetFreq',8,'testPlot',count);
figure
% axx1 = subplot(2,1,1);


[d,p] = compare_distributions(epoch_corrs.stim_corr(i,:),epoch_corrs.baseline_corr(i,:));
title(sprintf("%s %s\nSNR = %0.2f\n d=%0.2f, p=%0.3f", reg,stimInfo,SNR,d,p))

% axx2 = subplot(2,1,2);
hold on
full_trial = [epochs.baseline(:,:,i); epochs.signals(:,:,i); epochs.post_stim(:,:,i)];
plot(full_trial,'Color',[0.5 0.5 0.5])
plot(mean(full_trial,2),'Color',[1 0 0],'LineWidth',3)
hold off
% linkaxes([axx1,axx2],'x')
count =+1;
end


%%

signals = epochs.signals(:,:,401);
SNR = timeseries_SNR(signals,epoch_meta(401),'targetFreq',8,'testPlot');




%% SNR Method Testing
fs = 500;
endTime = 1;
t = linspace(0,1,fs*endTime);
testSig = cos(t*2*pi*8);
signals = repmat(testSig,24,1);
out = zeros(size(signals));
close all

for i=1:size(signals,1)
    x = circshift(signals(i,:),floor(pi*i/20));
    out(i,:) = x;
end
signals = rand(24,500);
SNR = timeseries_SNR(signals',fs,'targetFreq',8);


% figure(1)
% hold on
% plot(t,out')
% plot(t,testSig,'LineWidth',2)
% hold off

function [epochs,epoch_meta,epoch_corrs,pulseLocs]=load_data(dataPath)
epochs = load(fullfile(dataPath,'epochs_timeseries.mat'));
load(fullfile(dataPath,'epochs_metadata.mat'));
epoch_corrs = load(fullfile(dataPath,'epochs_pairwise_corrs.mat'));
pulseLocs = load(fullfile(dataPath,'epoch_pulseLocs.mat'));
pulseLocs = pulseLocs.pulseLocs;
end

function [intervals]=get_train_intervals(pulseLocs,n_trials,fs,trainFreq,stimLen)
n_trains = stimLen/fs * trainFreq;
N= length(pulseLocs);
stimLen= floor(fs / trainFreq*0.9); % number of samples in one cycle of the train
onsets  = zeros(n_trials,n_trains,N);
offsets = zeros(n_trials,n_trains,N);
for i=1:length(pulseLocs)
n = pulseLocs(i).n_pulse;
x = {pulseLocs(i).peaks.stim_array};
x = verify_symmetry(x);
if size(x,1) ~= 1
    x = x';
end
mat = cell2mat(x);
onsets(:,:,i) = mat(1:n:end,:)';
offsets(:,:,i) = mat(1:n:end,:)' + stimLen;
intervals = struct;
intervals.onsets = onsets;
intervals.offsets = offsets;
end

end



function out = verify_symmetry(dat)
lens = cellfun(@length,dat);
m = median(lens);
locs = find(lens ~= m);
out = dat;
for i=1:length(locs)
    subset = dat{locs(i)};
    diffs = diff(subset);
    [~, ml] = min(diffs);
    logIdx = ones(length(subset),1);
    logIdx(ml+1)=0;
    out{locs(i)} = subset(logical(logIdx));
    
end



end