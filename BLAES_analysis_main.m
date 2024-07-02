addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
Subject = 'BJH050';
run = 1;
brain=load(fullfile(boxpath,Subject,sprintf("%s_MNI.mat",Subject))); % MNI Brain
load(fullfile(boxpath,Subject,sprintf("%s_VERA_idx.mat",Subject))); %VERA_idx
dataPath = dir(sprintf("%s/%s/*.dat",boxpath,Subject));
tic
[signals,states,params] = load_bcidat(sprintf('%s/%s',dataPath.folder,dataPath.name));
states = parseStates(states);
fs = params.SamplingRate.NumericValue;

stimMap = parse_stimsweep_stimcodes(params);
timeit = toc;
fprintf("loading time elapsed %f s",timeit)

%% Downsampling and Channel Idenficiation
[states,signals,fs] = downsample_seeg(signals,states,fs,500);
channelLocs = find(~isnan(VERA_idx));
signals = signals(:,channelLocs);
signals = preprocess(signals);
[epochLocs,intervals] =getAllIntervals(states.StimulusCode,stimMap);
disp('done')
%%
timeOffset = 0.1;
tic;
epochs = epochData(signals,brain.electrodeNames,intervals,states.CereStimStimulation,timeOffset,fs,stimMap);
timeit = toc;
fprintf("epoch time elapsed %f s \n",timeit)
% exportRawData(fullfile(boxpath,Subject),seeg,states,params,intervals);

%%
fig = figure(1);
chan = 100;
for i=1:length(epochs(chan).epoch)
    ax=subplot(4,6,i);
    vals = epochs(chan).epoch(i).avg;
    std = epochs(chan).epoch(i).std;
    t = linspace(-1,1,length(vals));
    plot(ax,t,vals);
    hold on
    plotShading(ax,t,vals,std);
    title(epochs(chan).epoch(i).code);
end
hold off
sgtitle(epochs(chan).channel);




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
        x(i).timeOffset = timeoffset;
        x(i).signals = temp; % trials x samples
        x(i).stim_state = state;
        x(i).avg = mean(temp,1);
        x(i).std = std(temp,1);
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