addpath(genpath('/Users/nkb/Documents/NCAN/code/BLAES_stimSweep'))
addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
Subject = 'BJH050';
groupPath = fullfile(boxpath,'group');
sections = 0;
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
[signals,states,params] = load_bcidat(sprintf('%s/%s',dataPath.folder,dataPath.name));
states = parseStates(states);
fs = params.SamplingRate.NumericValue;

stimMap = parse_stimsweep_stimcodes(params);
fprintf('\ndone loading\n')
preprocessFlag = 0;
%% Preprocessing
% Downsampling and Channel Idenficiation
if ~preprocessFlag
    [states,signals,fs] = downsample_seeg(signals,states,fs,500);
    if exist(fullfile(boxpath,Subject,sprintf("%s_MNI_new.mat",Subject)),'file')
        disp('Found new VERA Struct')
        brain = load(fullfile(boxpath,Subject,sprintf("%s_MNI_new.mat",Subject)));
        channelKey = brain.electrodeNamesKey;
        channelLocs = ~cellfun(@isempty,{channelKey.VERANames});
        signals = signals(:,channelLocs);
        dat_channelNames = params.ChannelNames.Value(channelLocs);
        electrodeNames = brain.electrodes.Name;
        regions = brain.electrodes.Label;
        regions = cellfun(@(x) x{1},regions,'UniformOutput',false);

    elseif exist(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject)),'file')
        disp('Found VERA INDEX')
        load(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject))); %VERA_idx
        channelLocs = find(~isnan(VERA_idx));
        signals = signals(:,channelLocs);
        dat_channelNames = params.ChannelNames.Value(channelLocs);
        electrodeNames = brain.electrodeNames;
        regions = brain.SecondaryLabel;
        regions = cellfun(@(x) x{1},regions,'UniformOutput',false);
    else
        disp('No linking file between SEEG and Localization Found \n Removing ECG and Ref.')
        parseCell = {'REF1','REF2','EKG1','EKG2'};
        allchans = params.ChannelNames.Value;
        channelLocs = ~ismember(parseCell,allchans);
        signals = signals(:,channelLocs);
        dat_channelNames = params.ChannelNames.Value(channelLocs);
        elecrodeNames = dat_channelNames;
        regions = dat_channelNames;
    end
    preprocessFlag = 1;
end
timing_adjust = 80; % ms, due to software triggered stimulation state.
[signals, trigger,thresh] = preprocess(signals,fs,dat_channelNames,subject_info.Triggers);
[epochLocs,intervals] =getAllIntervals(states.StimulusCode,stimMap);
fprintf('\ndone preprocessing\n')
%% Epoching
pulseLocs= triggerEpochs(trigger,stimMap,timing_adjust,intervals,fs,8,1,thresh);
theta_epochs = theta_burst_epoch(signals,electrodeNames,intervals,timing_adjust,fs,stimMap,pulseLocs,regions);


theta_band = [4 10]; % Hz
smoothing_window = 0.1; % in seconds
filterOrder = computeStableFilter(theta_band,'bandpass',fs);
postStimStart = floor(2*fs);

parfor i=1:length(theta_epochs)
    [theta_epochs(i).theta_power,theta_epochs(i).theta_phase, theta_epochs(i).theta_filt] = timeseriesPower(theta_epochs(i).pre_stim_post',fs,theta_band,filterOrder, ...
        'smooth',smoothing_window,'baselineDuration',1);
    theta_epochs(i).theta_power = zscore(theta_epochs(i).theta_power);
    theta_epochs(i).theta_power = avg_baseline_correct(theta_epochs(i).theta_power,1,fs);
    theta_epochs(i).baseline_corr = channelCoherence(theta_epochs(i).baseline);
    theta_epochs(i).stim_corr = channelCoherence(theta_epochs(i).signals);
    theta_epochs(i).post_corr = channelCoherence(theta_epochs(i).pre_stim_post(:,postStimStart:end));
end
% configs  = {'current' 'pw' 'freq' 'loc'};
stimchannel_epochs = parseCondition(theta_epochs,'loc');
fprintf('\ndone epoching\n')
%% coherenence stats
export_labs = {'code' 'region' 'label' 'current' 'pw' 'freq' 'loc' 'channel' 'charge' 'chargePerPhase' 'chargeUnits' 'shank' 'snr_av' 'channel_idx' 'baseline_corr' 'stim_corr' 'post_corr'};

outstruct = struct();
N = length(theta_epochs);
parfor i=1:N
    for j=1:length(export_labs)
        field = export_labs{j};
        outstruct(i).(field) = theta_epochs(i).(field);
    end
    % stim vs baseline
    [d,p] = compare_distributions(theta_epochs(i).stim_corr,theta_epochs(i).baseline_corr);
    outstruct(i).stim_d = d;
    outstruct(i).stim_p = p;
    % post vs baseline
    [d,p] = compare_distributions(theta_epochs(i).post_corr,theta_epochs(i).baseline_corr);
    outstruct(i).post_d = d;
    outstruct(i).post_p = p;

end
savepath = fullfile(groupPath,sprintf('%s_cohens.mat',Subject));
save(savepath,"outstruct");
fprintf('\n coherence analyzed\n')
%% Fix Current Base Frequency
if sections
local_dat = stimchannel_epochs(1).entry([stimchannel_epochs(1).entry.current]==1);
analysis_name = "1 mA";
trajplot(local_dat,0,mainfigDir,fs,colors,analysis_name)
end

%% Charge Density
%% Single Channel
% close all
if sections
fig = figure(1);
chan = 100;
channel_name = brain.electrodeNames{chan};
local_data = theta_epochs(strcmp({theta_epochs.channel},channel_name));
tcl = tiledlayout(4,6);

for i=1:length(local_data)
    ax=nexttile(tcl);
    set(gca,'ButtonDownFcn',@fig_from_subplot)

    vals = local_data(i).avg;
    stdev = local_data(i).std;
    t = linspace(0,(length(vals)/fs),length(vals));
    hold on
    % plot(ax,t,vals,'Color',[1 0 0],'LineWidth',2);
    plot(ax,t,local_data(i).raw_sig','Color',[1,1,1 0.1])
    % val = max(abs(local_data(i).signals(:)));
    ylim(ax,[-1.01*val, 1.01*val]);
    plotShading(ax,t,vals,stdev);
    title(local_data(i).label);
end

hold off
linkaxes(tcl.Children,'x');
sgtitle(channel_name,'Interpreter','none');
end
%% Each condition per trajectory
if sections
close all
saveFigs = 0;
dir = fullfile(mainfigDir,'trajs');
trajplot(stimchannel_epochs(2).entry,0,dir,fs,colors);
end

%% Export Data for Method Testing
if sections
channels = [22 100 54 5 131];
exp_idx = ismember([theta_epochs.channel_idx],channels);
exportStruct = theta_epochs(exp_idx);
exportPath = '/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param/method_building';
save(fullfile(exportPath,'selected_data.mat'),"exportStruct");
save(fullfile(exportPath,'fs.mat'),'fs');
end
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
function trig=triggerEpochs(trigger,stimMap,timing_adjust,intervals,fs,base_frequency,stimDuration,thresh)
close all
timing_adjust_samps = floor((timing_adjust / 1000) * fs);
trig = struct();
stimCodes = [stimMap.Code];
% figure
% rows = 6;
% cols = ceil(length(intervals) / rows);
% tcl = tiledlayout(rows,cols);
for i=1:length(intervals)
    % ax = nexttile(tcl);
    onsets = intervals(i).start+timing_adjust_samps;
    offsets = intervals(i).stop+timing_adjust_samps;
    len = mean(offsets -onsets);
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
    % % wrote this to not repeat the loop, but preallocated loops are much
    % % faster than arrayfun
    segments = arrayfun(@(s, e) trigger(s:e), onsets, offsets, 'UniformOutput', false);
    segmentMatrix = vertcat(segments{:});
    segmentMatrix = reshape(segmentMatrix, [],length(segments));
    onsetPeaks = zeros(length(onsets),1);
    for j=1:length(onsets) % find the median first peak of the signal to index
        % TODO: add manual adjustment when number of peaks is off.
        % calculate number of peaks from train duration, frequency and num
        % pulses. 
        d = trigger(onsets(j):offsets(j));
        thresh = 0.5 * std(d);
        [pk,onsetPeaks(j)] = findpeaks(abs(d),'MinPeakDistance',peak_dist,'NPeaks',1,'MinPeakHeight',thresh);     
    end
    onsetPeak = mode(onsetPeaks);

    peaks = struct();
    for j=1:length(onsets)
        % TODO: add manual adjustment when number of peaks is off.
        % calculate number of peaks from train duration, frequency and num
        % pulses. 
        d = trigger(onsets(j):offsets(j));
        temp(j,:) = d;
        % thresh = 0.5*std(abs(d));
        % [pk,pkLoc] = findpeaks(abs(d),'MinPeakHeight',thresh,'MinPeakDistance',peak_dist,'NPeaks',num_peaks);
        % [pk,pkLoc] = findpeaks(abs(d),'MinPeakDistance',peak_dist);
        % [pk,pkLoc] = findpeaks(abs(d),'MinPeakDistance',peak_dist);
        % [pk,tloc] = maxk(pk,num_peaks);
        % pkLoc = pkLoc(tloc);
        
        pkLoc = stim_prediction(d,fs,trig(i).f,base_frequency,num_peaks, stimDuration,onsetPeak);

        % figure
        % plot(abs(d))
        % hold on
        % yline(thresh)
        % xline(linspace(1,8,8)*fs/8,'Color',[0 1 0])
        % scatter(pkLoc,pk);
        % hold off
        % pkLoc = pkLoc;
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

        elseif length(pk) > num_peaks && annotate
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
            peaks(j).n = find(pkLoc ==1);
            peaks(j).stim_array = pkLoc;
        end
    end
    trig(i).signal = temp;
    trig(i).peaks = peaks;
end
end
function x = theta_burst_epoch(signals,chan_lab,intervals,timing_adjust,fs,stimMap,triggerEpochs,regions)
timing_adjust_samps = floor((timing_adjust / 1000) * fs);
num_channels = size(signals,2);
total_combos = num_channels * length(intervals);
channelVector = repmat(linspace(1,num_channels,num_channels),length(intervals),1);
channelVector = channelVector(:);
indexVector = repmat(linspace(1,length(intervals),length(intervals)),num_channels,1);
x = struct();
stimCodes = [stimMap.Code];
idx = 1;
triggerCodes = [triggerEpochs.code];
frequencies = linspace(1,250,250);
for chan=1:num_channels
    region = regions{chan};
    % x = struct();
    for interval_idx=1:length(intervals)
        onsets = intervals(interval_idx).start+timing_adjust_samps;
        offsets = intervals(interval_idx).stop+timing_adjust_samps;
        len = mean(offsets -onsets);
        all_holder = zeros(length(onsets),3*(len+1));
        baseline_holder = zeros(length(onsets),len+1);
        signal_holder = zeros(length(onsets),len+1);
        raw_sig_holder = zeros(length(onsets),len+1);
        recombined_holder = zeros(length(onsets),3*(len+1));
        % state = zeros(length(onsets),len+1);
        x(idx).code = intervals(interval_idx).code;
        stimloc = triggerEpochs(ismember(triggerCodes,intervals(interval_idx).code));
        for trial_num=1:length(onsets)
            peakset = stimloc.peaks(trial_num).n;
            slice = signals(-len-1+onsets(trial_num):offsets(trial_num)+len+1,chan);
            slice_hp = getHighPassData(slice,2,4,fs);
            all_holder(trial_num,:) = slice_hp;
            d = signals(onsets(trial_num):offsets(trial_num),chan);
            baseline_window = [onsets(trial_num)-1-len,onsets(trial_num)-1];
            d = getHighPassData(d,2,4,fs);
            b = signals(baseline_window(1):baseline_window(2),chan);
            baseline = getHighPassData(b,2,4,fs);
            e = interpolateSpikes(d,peakset,stimloc.peaks(trial_num).stim_array,fs,10);
            br = baseline - baseline(end);
            sr = e - e(end);
            pr = slice_hp(floor(2*fs)+1:end) - slice_hp(floor(2*fs)+1);
            agg_denoised = [br; sr; pr];
            % figure(1)
            % plot(d,'LineWidth',4)
            % hold on
            % plot(e,'LineWidth',2)
            % hold off
            raw_sig_holder(trial_num,:) = d;
            signal_holder(trial_num,:) = e;
            baseline_holder(trial_num,:) = baseline;
            recombined_holder(trial_num,:) = agg_denoised;
            % state(trial_num,:) = stimulation_state(onsets(trial_num):offsets(trial_num));
        end
        label_loc = find(stimCodes == intervals(interval_idx).code);
        label_loc = label_loc(1);
        x(idx).label = makeLabelFromCereConfig(stimMap(label_loc).Config,stimMap(label_loc).Electrode);
        x(idx).region = region;
        x(idx).timing_adjust = timing_adjust;
        x(idx).signals = signal_holder; % trials x samples
        x(idx).raw_sig = raw_sig_holder;
        x(idx).baseline = baseline_holder;
        x(idx).pre_stim_post = all_holder;
        x(idx).full_trial = recombined_holder;
        % x(idx).stim_state = state;
        x(idx).avg = mean(signal_holder,1);
        x(idx).std = std(signal_holder,1);
        % [~,psd_stim] = matrix_PSD(signal_holder',frequencies,fs);
        % [~,psd_rest] = matrix_PSD(baseline_holder',frequencies,fs);
        % x(idx).psd_rest = psd_rest;
        % x(idx).psd_stim = psd_stim;
        % x(idx).freqz = frequencies;
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
        if isstrprop(chan_lab{chan}(2),'digit')
            traj = strcat(traj,chan_lab{chan}(2));
        end
        x(idx).shank = traj;
        [snr,~]= multi_SNR(signal_holder',8,fs,8);
        x(idx).snr_av = snr;
        % x(idx).snr_av = mean(snr);
        x(idx).channel_idx = chan;
        x(idx).stimloc = stimloc;
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
function [signals,triggerOut,thresh] = preprocess(signals, fs,electrodeNames, triggerElectrodes)
[b,a] = butter(6,75/(fs/2),"high");
% filtDat = filtfilt(b,a,signals);
triggerLocs = ismember(electrodeNames,triggerElectrodes);
test_locs = ~cellfun(@isempty,strfind(electrodeNames,'AR'));
common = getCommonAverage(signals(:,~triggerLocs));
common_agg = repmat(common,1,size(signals,2));
signals = signals - common_agg;
clear common_agg
trigger = abs(filtfilt(b,a,common));

trigger2 = mean(signals(:,triggerLocs),2); % this one seemingly works best, however will break down in bilateral stimulation. 
triggerOut = filtfilt(b,a,trigger2);
thresh = 0.5*std(triggerOut);
close all
% figure
% subplot(2,1,1)
% plot(abs(triggerOut))
% hold on
% yline(0.5*std(abs(trigger)),'LineWidth',3,'Color','r')
% hold off
% 
% subplot(2,1,2)
% plot(abs(trigger2))
% hold on
% yline(0.5*std(abs(trigger2)),'LineWidth',3,'Color','r')
% hold off
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