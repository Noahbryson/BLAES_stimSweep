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
data = load(fullfile(rootDataPath,Subject,"theta_epochs.mat"));


%% data viz
count = 1;
for i=111:112

loc = i;
reg = data.theta_epochs(loc).region;
stimInfo = data.theta_epochs(loc).label;
SNR = timeseries_SNR(data.theta_epochs(loc).signals',data.fs,'targetFreq',8,'testPlot',count);
figure
axx1 = subplot(2,1,1);
title(sprintf("%s %s", reg,stimInfo))
hold on
plot(data.theta_epochs(loc).pre_stim_post','Color',[0.5 0.5 0.5])
plot(mean(data.theta_epochs(loc).pre_stim_post),'Color',[1 0 0],'LineWidth',3)
hold off 
axx2 = subplot(2,1,2);
title(sprintf('SNR = %0.2f',SNR))
hold on
plot(data.theta_epochs(loc).full_trial','Color',[0.5 0.5 0.5])
plot(mean(data.theta_epochs(loc).full_trial),'Color',[1 0 0],'LineWidth',3)
hold off
linkaxes([axx1,axx2],'x')
count =+1;
end


%%

signals = data.theta_epochs(110).signals;
SNR = timeseries_SNR(signals',data.fs,'targetFreq',8,'testPlot');




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
