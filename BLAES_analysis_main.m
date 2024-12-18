addpath(genpath('/Users/nkb/Documents/NCAN/code/BLAES_stimSweep'))
addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
Subject = 'BJH050';
if exist(fullfile(boxpath,'Subject_Locations.xlsx'),'file')
    subject_info = readtable(fullfile(boxpath,'Subject_Locations.xlsx'));
    subject_info = table2struct(subject_info(strcmp(subject_info.Subject, Subject),:));
    stim_info = struct();
    stim_info(1).pair = subject_info.Pair1;stim_info(1).val = subject_info.loc1;stim_info(2).pair = subject_info.Pair2;stim_info(2).val = subject_info.loc2;
    subject_info.Triggers = strsplit(subject_info.Triggers,',');
    channelSort = 1;
else
    channelSort = 0;
end
mainfigDir = fullfile(boxpath,Subject,"figures");
if ~exist(mainfigDir,'dir')
    [~,~]=mkdir(mainfigDir);
end
run = 1;
brain=load(fullfile(boxpath,Subject,sprintf("%s_MNI.mat",Subject))); % MNI Brain

load('colors.mat');
dataPath = dir(sprintf("%s/%s/*.dat",boxpath,Subject));
tic
[signals,states,params] = load_bcidat(sprintf('%s/%s',dataPath.folder,dataPath.name));
states = parseStates(states);
fs = params.SamplingRate.NumericValue;

stimMap = parse_stimsweep_stimcodes(params);
timeit = toc;
fprintf("loading time elapsed %f s",timeit)
preprocessFlag = 0;
%% Preprocessing
% Downsampling and Channel Idenficiation
if ~preprocessFlag
    [states,signals,fs] = downsample_seeg(signals,states,fs,1000);
    if exist(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject)),'file')
        load(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject))); %VERA_idx
        channelLocs = find(~isnan(VERA_idx));
        signals = signals(:,channelLocs);
        dat_channelNames = params.ChannelNames.Value(channelLocs);
    end
    preprocessFlag = 1;
end
[signals, trigger] = preprocess(signals,fs,dat_channelNames,subject_info.Triggers);
[epochLocs,intervals] =getAllIntervals(states.StimulusCode,stimMap);
fprintf('\ndone\n')
% Epoching
timeOffset = 0.5;
pulseLocs= triggerEpochs(trigger,stimMap,timeOffset,intervals,fs,8,1);
theta_epochs = theta_burst_epoch(signals,brain.electrodeNames,intervals,states.CereStimStimulation,timeOffset,fs,stimMap,pulseLocs);
configs  = {'current' 'pw' 'freq' 'loc'};
conditions = struct();
%% sort by stimulation channel
stimchannel_epochs = parseCondition(theta_epochs,'loc');
%% Base Frequency
local_dat = stimchannel_epochs(1).entry([stimchannel_epochs(1).entry.current]==1);
trajplot(local_dat,0,dir,fs,colors)



%% Charge Density
%% Single Channel
% close all
fig = figure(1);
chan = 1;
channel_name = brain.electrodeNames{chan};
local_data = theta_epochs(strcmp({theta_epochs.channel},channel_name));
tcl = tiledlayout(4,6);

for i=1:length(local_data)
    ax=nexttile(tcl);
    set(gca,'ButtonDownFcn',@fig_from_subplot)

    vals = local_data(i).avg;
    stdev = local_data(i).std;
    t = linspace(0,local_data(i).timeOffset+(length(vals)/fs),length(vals));
    hold on
    plot(ax,t,vals,'Color',[1 0 0]);
    plot(ax,t,local_data(i).signals','Color',[0,0,0 0.001])
    ylim(ax,[-1.01*max(abs(local_data(i).signals(:))), 1.01*max(abs(local_data(i).signals(:)))]);
    plotShading(ax,t,vals,stdev);
    title(local_data(i).label);
end

hold off
linkaxes(tcl.Children,'x');
sgtitle(channel_name,'Interpreter','none');
%% Each condition per trajectory
close all
saveFigs = 0;
dir = fullfile(mainfigDir,'trajs');
trajplot(stimchannel_epochs(2).entry,0,dir,fs,colors);

%%

function result = parseCondition(epoch_struct,condition_fieldname)
conditions = unique([epoch_struct.(condition_fieldname)]);
result = struct();
for i=1:length(conditions)
    locs = [epoch_struct.(condition_fieldname)]==conditions(i);
    result(i).type = condition_fieldname;
    result(i).val = conditions(i);
    result(i).entry = epoch_struct(locs);
end

end
function exportRawData(fpath,seeg,states,params,intervals)
seegdat=seeg.data;
channels = seeg.name;
fname = fullfile(fpath,'stimsweep.mat');
save(fname,"seegdat","states",'-mat','-v7.3');
save(fullfile(fpath,'channels.mat'),"channels", "-mat",'-v7');
save(fullfile(fpath,'params.mat'),"params", "-mat",'-v7');
save(fullfile(fpath,'intervals.mat'),"intervals", "-mat",'-v7');

end
function trig=triggerEpochs(trigger,stimMap,timeoffset,intervals,fs,base_frequency,stimDuration)
sampleOffset = timeoffset * fs;
trig = struct();
stimCodes = [stimMap.Code];
for i=1:length(intervals)
    
    onsets = intervals(i).start;
    offsets = intervals(i).stop;
    len = mean(offsets -onsets) + sampleOffset;
    temp = zeros(length(onsets),len+1);
    
    
    trig(i).code = intervals(i).code;
    label_loc = find(stimCodes == intervals(i).code);
    label_loc = label_loc(1);
    vals = [stimMap(label_loc).Config.Config];
    trig(i).current = str2double(vals{3})/1000;
    trig(i).pw = str2double(vals{5});
    trig(i).f = str2double(vals{7});
    trig(i).n_pulse = str2double(vals{2});
    peak_dist = floor(fs/trig(i).f) -1;
    num_peaks = base_frequency*stimDuration*trig(i).n_pulse;
    peaks = struct();
    for j=1:length(onsets)
        % TODO: add manual adjustment when number of peaks is off.
        % calculate number of peaks from train duration, frequency and num
        % pulses. 
        d = trigger(onsets(j):offsets(j)+sampleOffset);
        temp(j,:) = d;
        thresh = 1.5*std(abs(d));
        [pk,pkLoc] = findpeaks(abs(d(100:end)),'MinPeakHeight',thresh,'MinPeakDistance',peak_dist);
        pkLoc = pkLoc+100;
        annotate = 0;
        if length(pk) < num_peaks && annotate
            n = num_peaks - length(pk);
            fig=figure;
            hold on
            plot(abs(d),'LineWidth',2.2)
            yline(thresh)
            scatter(pkLoc,pk);
            title(sprintf('code %d iter %d\nadd %d points',i,j,n))
            hold off
            set(gcf,'Position',[800 600 2000 900])
            [pkLoc_n,pk_n] = ginput(n);
            close(fig)
            pkLoc = [pkLoc; pkLoc_n];
            pk = [pk; pk_n];
            peaks(j).n = pkLoc;

        elseif length(pk) > num_peaks
            n = length(pk)-num_peaks;
            fig=figure;
            hold on
            plot(abs(d),'LineWidth',2.2)
            yline(thresh)
            scatter(pkLoc,pk);
            title(sprintf('code %d iter %d\nremove %d points',i,j,n))
            hold off
            set(gcf,'Position',[800 600 2000 900])
            [pkLoc_n,pk_n] = ginput(n);
            pkLoc_n = pkLoc(knnsearch(pkLoc(:),pkLoc_n(:)));
            rm_loc = ismember(pkLoc, pkLoc_n);
            pkLoc = pkLoc(~rm_loc);
            pk = pk(~rm_loc);
            peaks(j).n = pkLoc;
            close(fig)
        else
            peaks(j).n = pkLoc;
        end
    end
    trig(i).signal = temp;
    trig(i).peaks = peaks;
end
end
function x = theta_burst_epoch(signals,chan_lab,intervals,stimulation_state,timeoffset,fs,stimMap,triggerEpochs)
sampleOffset = timeoffset * fs;
num_channels = size(signals,2);
x = struct();
stimCodes = [stimMap.Code];
idx = 1;
triggerCodes = [triggerEpochs.code];

for chan=1:num_channels
    % x = struct();
    for interval_idx=1:length(intervals)
        onsets = intervals(interval_idx).start;
        offsets = intervals(interval_idx).stop;
        len = mean(offsets -onsets) + sampleOffset;
        signal_holder = zeros(length(onsets),len+1);
        raw_sig_holder = zeros(length(onsets),len+1);
        % state = zeros(length(onsets),len+1);
        x(idx).code = intervals(interval_idx).code;
        stimloc = triggerEpochs(ismember(triggerCodes,intervals(interval_idx).code));
        for trial_num=1:length(onsets)
            peakset = stimloc.peaks(trial_num).n;
            d = signals(onsets(trial_num):offsets(trial_num)+sampleOffset,chan);
            d = getHighPassData(d,0.5,2,fs);
            e = interpolateSpikes(d,peakset,stimloc.signal(trial_num,:));
            % figure(1)
            % plot(d,'LineWidth',4)
            % hold on
            % plot(e,'LineWidth',2)
            % hold off
            raw_sig_holder(trial_num,:) = d;
            signal_holder(trial_num,:) = e;
            % state(trial_num,:) = stimulation_state(onsets(trial_num):offsets(trial_num)+sampleOffset);
        end
        label_loc = find(stimCodes == intervals(interval_idx).code);
        label_loc = label_loc(1);
        x(idx).label = makeLabelFromCereConfig(stimMap(label_loc).Config,stimMap(label_loc).Electrode);
        x(idx).timeOffset = timeoffset;
        x(idx).signals = signal_holder; % trials x samples
        x(idx).raw_sig = raw_sig_holder;
        % x(idx).stim_state = state;
        x(idx).avg = mean(signal_holder,1);
        x(idx).std = std(signal_holder,1);
        vals = [stimMap(label_loc).Config.Config];
        ma = str2double(vals{3})/1000;
        pw = str2double(vals{5});
        f = str2double(vals{7});
        n_pulse = str2double(vals{2});
        x(idx).current = ma;
        x(idx).pw = pw;
        x(idx).freq=f;
        x(idx).loc = stimMap(label_loc).Electrode;
        x(idx).channel = chan_lab{chan};
        electrodeSA = 5*0.01; % cm^2
        stimDur = 1; %s
        x(idx).charge = pw*ma*1000*n_pulse*f*stimDur/electrodeSA *1e-6;
        x(idx).chargePerPhase = pw*ma*1000/electrodeSA *1e-6;
        x(idx).chargeUnits = "uC/cm^2";
        % charge = us * uA * numPulses * Hz * s/cm^2 = pC/cm^2 -> pC/cm^2 * 1e-6 = uC/cm^2
        % charge/phase = us * uA/cm^2= pC/cm^2-> pC/cm^2 *1e-6 = uC/cm^2
        if ~isempty(strfind(chan_lab{chan},"'"))
        traj = sprintf('%s_L',chan_lab{chan}(1));
        else
            traj = sprintf('%s_R',chan_lab{chan}(1));
        end
        x(idx).shank = traj;
        [snr,~]= getSNR(signal_holder,8,fs,8);
        x(idx).snr = snr;
        x(idx).snr_av = mean(snr);
        x(idx).channel_idx = chan;
        idx = idx+1;
    end
    % x(k).epoch = x;
end
end

function [epochLocs,intervals] = getAllIntervals(stimcode,stimMap)
epochLocs = zeros(size(stimMap));
for i=1:length(epochLocs)
    epochLocs(i)=stimMap(i).Code;
end
epochLocs = unique(epochLocs);
intervals = struct;
rm_rows = [];
for i=1:length(epochLocs)
    intervals(i).code = epochLocs(i);
    temp = getInterval(stimcode,epochLocs(i));
    if ~isempty(temp)

        intervals(i).start = temp(:,1);
        intervals(i).stop = temp(:,2);
    else
        rm_rows = [rm_rows i];
    end
end
mask = true(size(epochLocs));
mask(rm_rows) = false;
epochLocs = epochLocs(mask);
intervals = intervals(mask);
end
function intervals = getInterval(stimCode,value)
A = find(stimCode==value);
differences = diff(A);

% Find indices where the difference is greater than 1
if ~isempty(differences)
    breaks = find(differences > 1);
    intervals = zeros(length(breaks),2);
    % Define start and end indices of each interval
    idx = A(1);
    for i=1:length(breaks)
        intervals(i,:) = [idx,A(breaks(i))];
        idx = A(breaks(i)+1);
    end
    % disp(intervals);
else
    intervals = [];
end
end
function [signals,trigger2] = preprocess(signals, fs,electrodeNames, triggerElectrodes)
[b,a] = butter(2,85/(fs/2),"high");
% filtDat = filtfilt(b,a,signals);
triggerLocs = ismember(electrodeNames,triggerElectrodes);
test_locs = ~cellfun(@isempty,strfind(electrodeNames,'AR'));

% close all
% figure
% subplot(2,1,1)
% plot(abs(trigger))
% hold on
% yline(0.5*std(abs(trigger)),'LineWidth',3,'Color','r')
% hold off

% subplot(2,1,2)
% plot(abs(trigger2))
% hold on
% yline(0.5*std(abs(trigger2)),'LineWidth',3,'Color','r')
% hold off
[b1,a1] = butter(6,65/(fs/2),"high");


common = getCommonAverage(signals(:,~triggerLocs));
common_agg = repmat(common,1,size(signals,2));
clear common
signals = signals - common_agg;
trigger = mean(signals(:,triggerLocs),2);
trigger = filtfilt(b,a,trigger);
trigger3 = filtfilt(b,a,signals(:,triggerLocs));

trigger2 = mean(signals(:,test_locs),2);
trigger2 = filtfilt(b,a,trigger2);
end
function x  = getCommonAverage(signals)
x = mean(signals,2);
end
function [states,signals,fs_out]=downsample_seeg(signals,states,fs,newfs)
ratio = round(fs/newfs,0);
fs_out = fs/ratio;
signals = signals(1:ratio:end,:);
fields = fieldnames(states);
for i=1:length(fields)
    states.(fields{i}) = states.(fields{i})(1:ratio:end);
end


end