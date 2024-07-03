addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
Subject = 'BJH050';
subject_info = readtable(fullfile(boxpath,'Subject_Locations.xlsx'));
subject_info = table2struct(subject_info(strcmp(subject_info.Subject, Subject),:));
stim_info = struct();
stim_info(1).pair = subject_info.Pair1;stim_info(1).val = subject_info.loc1;stim_info(2).pair = subject_info.Pair2;stim_info(2).val = subject_info.loc2;
mainfigDir = fullfile(boxpath,Subject,"figures"); 
if ~exist(mainfigDir,'dir')
    [~,~]=mkdir(mainfigDir);
end
run = 1;
brain=load(fullfile(boxpath,Subject,sprintf("%s_MNI.mat",Subject))); % MNI Brain
load(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject))); %VERA_idx
load('colors.mat');
dataPath = dir(sprintf("%s/%s/*.dat",boxpath,Subject));
tic
[signals,states,params] = load_bcidat(sprintf('%s/%s',dataPath.folder,dataPath.name));
states = parseStates(states);
fs = params.SamplingRate.NumericValue;

stimMap = parse_stimsweep_stimcodes(params);
timeit = toc;
fprintf("loading time elapsed %f s",timeit)

%% Downsampling and Channel Idenficiation
[states,signals,fs] = downsample_seeg(signals,states,fs,1000);
channelLocs = find(~isnan(VERA_idx));
signals = signals(:,channelLocs);
signals = preprocess(signals);
[epochLocs,intervals] =getAllIntervals(states.StimulusCode,stimMap);
disp('done')
%% Epoching
timeOffset = 0.5;
tic;
epochs = epochData(signals,brain.electrodeNames,intervals,states.CereStimStimulation,timeOffset,fs,stimMap);
timeit = toc;
fprintf("epoch time elapsed %f s \n",timeit)
% exportRawData(fullfile(boxpath,Subject),seeg,states,params,intervals);

%% Single Channel
close all
fig = figure(1);
chan = 69;
for i=1:length(epochs(chan).epoch)
    ax=subplot(4,6,i);
    set(gca,'ButtonDownFcn',@fig_from_subplot)

    vals = epochs(chan).epoch(i).avg;
    std = epochs(chan).epoch(i).std;
    t = linspace(-epochs(chan).epoch(i).timeOffset,epochs(chan).epoch(i).timeOffset+(length(vals)/fs/2),length(vals));
    plot(ax,t,vals);
    hold on
    plotShading(ax,t,vals,std);
    title(epochs(chan).epoch(i).label);
end
hold off
sgtitle(epochs(chan).channel,'Interpreter','none');
%% Each condition per trajectory
close all
dir = fullfile(mainfigDir,'trajs');
if ~exist(dir,'dir')
    mkdir(dir);
end
trajectories = unique({epochs.shank});
traj_idx = {epochs.shank};
for tr=1:length(trajectories)
% for tr=1:1

    locs = find(strcmp(traj_idx,trajectories{tr}));
    fig=figure('Visible','off');
    set(fig,"PaperSize",[8 11]);
    fig.PaperPosition = [0 0 8 11];
    set(gcf,'Position',[100 100 3000 2000])
    tcl = tiledlayout(4,4);
    for l=1:length(locs)
        chan = epochs(locs(l)).channel;
        data = epochs(locs(l)).epoch;
        % ax = subplot(4,4,l);
        ax=nexttile(tcl);
        set(gca,'ButtonDownFcn',@fig_from_subplot)
        hold on
        for j=1:length(data)
            y = data(j).avg + j;
            std = data(j).std;
            t = linspace(-data(j).timeOffset,data(j).timeOffset + length(y)/fs/2,length(y));
            plot(ax,t,y,'LineWidth',2,'DisplayName',data(j).label,'Color',colors(j,:));
            % ax = plotShading(ax,t,y,std);
        end
        hold off
        tit = title(chan,'Interpreter','none');
        fontsize(tit,24,'points')
    end
    hL = legend({data.label});
    fontsize(hL,20,'points')
    hL.Layout.Tile='East';
    exportgraphics(fig,fullfile(dir,sprintf("%s traj.pdf",trajectories{tr})),"Append",false);
    exportgraphics(fig,fullfile(dir,'trajs.pdf'),"Append",true);

    % hL.Location='eastoutside';
end



%%
function exportRawData(fpath,seeg,states,params,intervals)
seegdat=seeg.data;
channels = seeg.name;
fname = fullfile(fpath,'stimsweep.mat');
save(fname,"seegdat","states",'-mat','-v7.3');
save(fullfile(fpath,'channels.mat'),"channels", "-mat",'-v7');
save(fullfile(fpath,'params.mat'),"params", "-mat",'-v7');
save(fullfile(fpath,'intervals.mat'),"intervals", "-mat",'-v7');

end

function X = epochData(signals,chan_lab,intervals,stimulation_state,timeoffset,fs,stimMap)
sampleOffset = timeoffset * fs;
num_channels = size(signals,2);
X = struct();
stimCodes = [stimMap.Code];
for k=1:num_channels
    x = struct();
    for i=1:length(intervals)
        onsets = intervals(i).start;
        offsets = intervals(i).stop;
        len = mean(offsets -onsets) + sampleOffset*2;
        temp = zeros(length(onsets),len+1);
        state = zeros(length(onsets),len+1);
        for j=1:length(onsets)
            temp(j,:) = signals(onsets(j)-sampleOffset:offsets(j)+sampleOffset,k);
            state(j,:) = stimulation_state(onsets(j)-sampleOffset:offsets(j)+sampleOffset);
        end
        x(i).code = intervals(i).code;
        label_loc = find(stimCodes == intervals(i).code);
        label_loc = label_loc(1);
        x(i).label = makeLabelFromCereConfig(stimMap(label_loc).Config,stimMap(label_loc).Electrode);
        x(i).timeOffset = timeoffset;
        x(i).signals = temp; % trials x samples
        x(i).stim_state = state;
        x(i).avg = mean(temp,1);
        x(i).std = std(temp,1);
        vals = [stimMap(label_loc).Config.Config];
        ma = str2num(vals{3})/1000;
        pw = str2num(vals{5});
        f = str2num(vals{7});
        x(i).current = ma;
        x(i).pw = pw;
        x(i).freq=f;
        x(i).loc = stimMap(label_loc).Electrode;
    end
    X(k).epoch = x;
    X(k).channel = chan_lab{k};
    if ~isempty(strfind(chan_lab{k},"'"))
        traj = sprintf('%s_L',chan_lab{k}(1));
    else
        traj = sprintf('%s_R',chan_lab{k}(1));
    end

    X(k).shank = traj;

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



function signals = preprocess(signals,fs)
% [b,a] = butter(2,0.5/(fs/2),"high");
% filtDat = filtfilt(b,a,signals);
common = getCommonAverage(signals);
common_agg = repmat(common,1,size(signals,2));
clear common
signals = signals - common_agg;
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