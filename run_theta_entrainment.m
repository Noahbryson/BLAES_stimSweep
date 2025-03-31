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
Subject = allSubList{3};
[epochs,epoch_meta,epoch_corrs,pulseLocs]=load_data(fullfile(rootDataPath,Subject,'processed'));
fs = epochs.fs;
n_trials = 24;

stimLen = size(epochs.signals,1);
% train_intervals = get_train_intervals(pulseLocs,n_trials,fs,8,stimLen);
trainFreq = 8; % hz
intervals=make_train_intervals(fs,trainFreq,stimLen);
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


%% Temporal SNR
close all
idx = 119;
numFIg = 4;
% perm = randperm(length(epoch_meta)-numFIg,1);
perms = perm:1:perm+numFIg;
figure('Name',sprintf('%s_%d_temporal-SNR',Subject,perm))
% tcl = tiledlayout(1,6);
if numFIg <= 6
tcl = tiledlayout(1,numFIg,"Padding","compact");
else
tcl = tiledlayout('flow',"Padding","compact");    
end

for i=1:numFIg
    % ax = nexttile(tcl);
    sublayout = tiledlayout(tcl,2,1,'Padding','tight');
    sublayout.Layout.Tile=i;
    idx=perms(i);
    signals = epochs.signals(:,:,idx);
    lab = epoch_meta(idx).label;
    reg = epoch_meta(idx).region;
    chan = epoch_meta(idx).channel;
    tit = sprintf("During Stim Response\n%s %s\n%s",chan,reg,lab);
    % SNR = timeseries_SNR(signals,fs,'targetFreq',8,'testPlot');
    [t_t,SNR_t] = SNR_over_time(signals,fs,intervals,'targetFreq',8);
    SNR_t = SNR_t -1;
    [t_C,Coh_t] = temporal_channelCoherence(signals,fs,intervals);
    psignal = mean(signals,2);
    t = linspace(0,size(signals,1)/fs,size(signals,1));
    ba = nexttile(sublayout);
    plot(t,signals,'Color',[.5 .5 .5]);
    title(tit)
    ylabel('\muV')
    hold on
    plot(t,psignal,'LineWidth',2);
    hold off
    bq=nexttile(sublayout);
    hold on
    plot(mean(t_t,2),SNR_t,'r')
    scatter(mean(t_t,2),SNR_t,'r')
    plot(mean(t_C,2),Coh_t,'b')
    scatter(mean(t_C,2),Coh_t,'b')
    hold off
    title('SNR vs time')
    xlabel('time (s)')
    ylabel('SNR (a.u.)')

    minval=-.1;
    maxval=1.1;
    if min(SNR_t) < -0.1
        minval = min(SNR_t);
    end
    if max(SNR_t) > 1.1
        maxval = max(SNR_t);
    end
    
    ylim([minval,maxval])
    % if i==numFIg

    % end
            leg = legend({'adjusted SNR','','Pairwise-Corr',""},'Location','SouthOutside');
    % linkaxes([ba,bq],'x')
end



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

function intervals=get_train_intervals(pulseLocs,n_trials,fs,trainFreq,stimLen)
n_trains = stimLen/fs * trainFreq;
N= length(pulseLocs);
stimLen= floor(fs / trainFreq*0.9); % number of samples in one cycle of the train
onsets  = zeros(n_trials,n_trains,N);
offsets = zeros(n_trials,n_trains,N);
for i=1:length(pulseLocs)
    disp(i)
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
function intervals=make_train_intervals(fs,trainFreq,stimLen)
stimDuration = stimLen/fs;
trainLen = floor(fs/trainFreq);
trains_per_run = round(stimLen / trainLen);
intervals = zeros(trains_per_run,2);
remain = 0;
if mod(stimLen,trains_per_run) ~= 0
    remain = stimLen - trains_per_run * trainLen;
    % stimLen = stimLen-remain;
    stimLen = trains_per_run * trainLen;
end
mat = linspace(1,stimLen,stimLen);
subset = reshape(mat, [],trains_per_run);
intervals(:,1) = subset(1,:);
intervals(:,2) = subset(end,:);
intervals(end,2) = intervals(end,2)+remain;




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