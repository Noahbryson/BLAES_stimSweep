%% BCI2000 BLAES Stim Parameter Sweep .prm Generator
%
% A sequence is created to alternate the fixation cross stimuli with the
% image stimuli.
%
% The stimuli and meaningful parameters are written into a param
% variable and stored as a *.prm file using the convert_bciprm function.
% 
% The User Parameters Section is the only thing that should be edited in this
% script
% The current version requires two channels to be used as input, future
% versions will allow for 1-48 channels to be selected (although 48
% channels would be a bit excessive, and is the max supported by the cerestim)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Noah Bryson <noahbryson@neurotechcenter.org>
%
% $BEGIN_BCI2000_LICENSE$
%
% This file is part of BCI2000, a platform for real-time bio-signal research.
% [ Copyright (C) 2000-2021: BCI2000 team and many external contributors ]
%
% BCI2000 is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% BCI2000 is distributed in the hope that it will be useful, but
%                         WITHOUT ANY WARRANTY
% - without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
%
% $END_BCI2000_LICENSE$
% http://www.bci2000.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TODO
% test parameter files with cerestim to ensure stimulation is happening
%   when and how I predict it will. 
%   - it may be useful to do a recording on the gUSB amp to measure this lol.
%
% make CCEP configurations dynamic based on the number of channels.
%   Currently CCEPs are statically set for two channels, this should not be
%   hard to implement I am just tired of working on this.
%
% make including CCEPs optional, and add an option for only including CCEPs
% on every other block. 
%
% make an option to include user feedback so I can use the same script to
%   generate Tao's stimulation titration experiment
%
% make train frequency and duration configurable
%
% make number of pulses able to be a user entry component, dynamically
%   assign if not

%% Notes
% For best performance, use a sample block of 0.1 x sampling rate (100 ms)

%% User Parameters

% timing
%   1 Block takes (6s x number of configurations)
%   This is preserved with and without CCEPs. CCEPs add an additional 2
%   seconds, so the baseline period between 
%   NOTE: The start of a new BLOCK will have a 10s delay.

% exclusions
%   excluding a 1 mA condition will cause this to break unfortunately,
%   since the pseduorandomization uses recursive methods to ensure rules
%   about amplitude and channels are not violated. 
                                    
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools'; % Set the path of the BCI2000 main directory here, need this to write params
% noah windows path 'C:\BCI2000\BCI2000'
% noah mac path: '/Users/nkb/Documents/BCI2000tools'
% utah path
% bigcart path
% smallcart path
% enter stim channels, different at WU vs UU as WU plugs in hardware
% manually due to number of recording channels.
cathodeChannels = [1, 3]; % cathode leading stim channels passed to stimulator 
anodeChannels   = [2, 4]; % anode leading stim channels passed to stimulator
%       be sure cathode and anode channels are aligned with each other for
%       bipolar pairs. 
stimAmps = [1000,2000]; %uA
pulseWidth = [250, 500]; %us, single phase, double for entire biphasic stim
frequencies = [33, 50, 80]; %Hz
num_pulses =  [3,  4,  6]; % optional parameter, fill in for number of pulses associated with each primary frequency. 
%       If the numel between the two is not equal, num_pulses will be dynamically assigned
electrodeSurfaceArea = 5; %mm^2, should be the same for all sites, but verify to be sure. 
carrier_freq = 8; 
num_blocks = 60; % number of trials for each condition, each block contains one trial of each condition. 
video_num = 1; % video index, name and number below
%  1:   1 Year of Growing Food - A whole season of vegetable gardening.mp4'
%  2:   2 Hours of Digital Art.mp4'
%  3:   500 Days Survival And Build In A Rain Forest - NO FOOD, NO WATER, NO SHELTER.mp4'
%  4:   Bob Ross - One Hour Special - The Grandeur of Summer.mp4'
%  5:   y2mate.com - Great Inventions  60 Minutes Full Episodes_720p.mp4'                }

conditions2remove = [];% add the numeric value as they appear on the testing image
generateTest = 1; % set to 1 to regenerate testing sequence and image, note it will not generate if any conditons are excluded. 
rmHighestChargeCondition = 1; % set to 1 if removing highest amplitude-pulsewidth pair
allowCCEPs = 1; % set to 1 for CCEPs leading and lagging the stimulation pulses, 0 for no CCEPs
CCEP_amplitude = 1000; % uA
CCEP_ISI = 1; % sec, duration following a CCEP 
% Pathing
% root = 'C:\Paradigms'; % windows path
% stimuliDir = fullfile(root,'\tasks\BLAES\BLAES_param_sweep\stimuli'); % path to folder containing the videos used in this experiment
% checkDir(stimuliDir);
% parmDir  = fullfile(root,'\parms\BLAES\_BLAES_param_sweep'); % path to write parameter file to
% checkDir(parmDir);

% mac pathing adjust
root = '/Users/nkb/Documents/Paradigms'; % mac path
stimuliDir = fullfile(root,'tasks','BLAES_param_sweep','stimuli'); % path to folder containing the videos used in this experiment
checkDir(stimuliDir);
parmDir  = fullfile(root,'parms','BLAES','_BLAES_param_sweep'); % path to write parameter file to
checkDir(parmDir);
videos = dir(fullfile(stimuliDir,'*.mp4'));
% 
parmName_dec = videos(video_num).name(1:20);
videoPath = fullfile(videos(video_num).folder,videos(video_num).name);

% parmName_dec = 'test';
% videoPath = 'test';
bci2ktools(BCI2KPath) % sets up access to BCI2000 functions and .mex files, just change local dir above




%% Keyboard Map
keyboard = struct;
temp = cellfun(@upper,{'q' 'w' 'e' 'r' 't' 'y' 'u' 'i' 'o' 'p' 'a' 's' 'd' 'f' 'g' 'h' 'j' 'k' 'l' 'z' 'x' 'c' 'v' 'b' 'n' 'm'},'UniformOutput', false);
keyboard.keys= temp;
keyboard.X = [0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7]*10;
keyboard.Y = [2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]*20;
load("keys2numbers.mat") % load a keyboard mapping file which is apart of this repository. 
test = cellfun(@num2str, keyboard.keys, 'UniformOutput', false);
test2 = cellfun(@num2str, fullKeyMap(1,:), 'UniformOutput', false);
logicalIdx = ismember(test2,test);
keyMap = fullKeyMap(:,logicalIdx);
%% Video Information
video = VideoReader(videoPath);
vid_duration = video.Duration; % duration of the video in seconds
clear video
%% Generate Potential Conditions
% 
% 
% 
if numel(frequencies) ~= numel(num_pulses)
    num_pulses = zeros(size(frequencies));
    for k=1:length(frequencies)
        num_pulses(k) = floor(frequencies(k)/(2*carrier_freq))+1;
    end
end
stimMat_cathode = [];
stimLabel = [];
stimMat_cols = {'cathodeFirst' 'numPulses' 'pulseWidth' 'amp' 'channel' 'freq','stimCode','config_ID'};
% Channel Condition
channelCol = 5;
% stim Amp Condition
ampCol = 4;
stimCode = 9; % stimCode 1&2 reserved for AV start sequence
configID = 1;

counter =0;
for loc=1:length(pulseWidth)
    for j=1:length(stimAmps)
        for k=1:length(frequencies)
            for q=1:length(cathodeChannels)

                if pulseWidth(loc) > 250 && stimAmps(j) >= max(stimAmps) && rmHighestChargeCondition
                    % do not pair 2 mA with 500 us PW
                else
                    temp = [1,num_pulses(k),pulseWidth(loc),stimAmps(j),cathodeChannels(q),frequencies(k),stimCode,configID];
                    counter = counter+1;
                    stimCode = stimCode+1;
                    stimMat_cathode = [stimMat_cathode; temp];
                    stimLabel = [stimLabel; counter];
                    if q==2
                        configID = configID+ 1;
                    end
                end
            end
        end
    end
end
stimMat_cathode = rm_combinations(stimMat_cathode,conditions2remove,1);
stimLabel = rm_combinations(stimLabel,conditions2remove,0);
stimMat_anode = stimMat_cathode;
stimMat_anode(:,1) = stimMat_anode(:,1) - 1; %switch to anodic leading
stimMat_anode(:,8) = stimMat_anode(:,8) + max(stimMat_cathode(:,8)); %adjust configIDs


% adjust channel to make bipolar pair
for loc=1:size(stimMat_anode,1)
    idx = cathodeChannels == stimMat_anode(loc,channelCol);
    stimMat_anode(loc,channelCol) = anodeChannels(idx);
end
numConfigs = size(stimMat_cathode,1);
numParams = size(stimMat_cathode,2);
seq = pseduoRandReplace(stimMat_cathode,0,0,0,ampCol,channelCol)';
checkRepeats(seq,channelCol); % channel column, ensures no channels are repeated.
checkRepeatVal(seq,ampCol,max(stimAmps)); % amplitude column, ensures 2 mA is not repeated. 
%% Generate Trial Sequence
start = 1;
stop = numConfigs;
trial_seq = zeros(numConfigs*num_blocks,numParams);
for loc=1:num_blocks
    
    if loc == 1
        multi = 0;
        amp_init = 0;
        chan_init = 0;
    else
        multi = 1;
        amp_init = trial_seq(stop-numConfigs,ampCol);
        chan_init = trial_seq(stop-numConfigs,channelCol);
    end

    % trial_seq(start:stop,:) = pseduoRandReplace(stimMat_cathode,multi,amp_init,chan_init,ampCol,channelCol);
    trial_seq(start:stop,:) = pseduoRandReplace(stimMat_cathode,1,0,chan_init,ampCol,channelCol);
    % if i < numTrials
    %     trial_seq(stop+1,:) = [0 0 0 0 0 0 stimCode 0];
    % end
    start = stop+1;
    stop = start+numConfigs-1;
    
    
end
checkRepeats(trial_seq,channelCol); % channel column
checkRepeatVal(trial_seq,ampCol,max(stimAmps)); % amplitude column

%% Trial Sweep Parm Generation
param = struct;
% set up Stim Configs
rowLabs = {'Cathode First';'Number of pulses';'Phase 1 amp (uA)';'Phase 2 amp (uA)';'Phase 1 duration (us)';'Phase 2 duration (us)';'Frequency (Hz)';'Interphase duration (us)';'Train Duration (s)';'Train Frequency (Hz)'};
stim_configs = getUniqueConfigs(stimMat_cathode,stimMat_anode);
param.StimulationConfigurations.Section = 'CereStim';
param.StimulationConfigurations.Type = 'matrix';
param.StimulationConfigurations.DefaultValue = '';
param.StimulationConfigurations.LowRange = '';
param.StimulationConfigurations.HighRange = '';
param.StimulationConfigurations.Comment = 'Configurations for BLAES PARAM SWEEP';
param.StimulationConfigurations.Value = cell(length(rowLabs),size(stim_configs,1)+2);

for loc=1:size(stim_configs,1)
    param.StimulationConfigurations.Value{1, loc}  = sprintf('%d',stim_configs(loc,1));%cathode bool
    param.StimulationConfigurations.Value{2, loc}  = sprintf('%d',stim_configs(loc,2));% n pulses
    param.StimulationConfigurations.Value{3, loc}  = sprintf('%d',stim_configs(loc,4));% amp 1
    param.StimulationConfigurations.Value{4, loc}  = sprintf('%d',stim_configs(loc,4));% amp 2
    param.StimulationConfigurations.Value{5, loc}  = sprintf('%d',stim_configs(loc,3));%pw 1
    param.StimulationConfigurations.Value{6, loc}  = sprintf('%d',stim_configs(loc,3));%pw 2
    param.StimulationConfigurations.Value{7, loc}  = sprintf('%d',stim_configs(loc,6));%stim freq 
    param.StimulationConfigurations.Value{8, loc}  = sprintf('%d',53);%inter pulse duration (53 us min)
    param.StimulationConfigurations.Value{9, loc}  = sprintf('%d',1);% train duration 
    param.StimulationConfigurations.Value{10,loc}  = sprintf('%d',carrier_freq);% train freq 
    param.StimulationConfigurations.ColumnLabels{loc,1} = sprintf('Configuration %d',loc);
end
param.StimulationConfigurations.Value{1, end-1}  = sprintf('%d', 1);                %cathode bool
param.StimulationConfigurations.Value{2, end-1}  = sprintf('%d', 1);                % n pulses
param.StimulationConfigurations.Value{3, end-1}  = sprintf('%d', CCEP_amplitude); % amp 1
param.StimulationConfigurations.Value{4, end-1}  = sprintf('%d', CCEP_amplitude); % amp 2
param.StimulationConfigurations.Value{5, end-1}  = sprintf('%d', 250);              %pw 1
param.StimulationConfigurations.Value{6, end-1}  = sprintf('%d', 250);              %pw 2
param.StimulationConfigurations.Value{7, end-1}  = sprintf('%d', 50);               %stim freq (does not matter for SPES)
param.StimulationConfigurations.Value{8, end-1}  = sprintf('%d', 53);               %inter pulse duration (53 us min)
param.StimulationConfigurations.Value{9, end-1}  = sprintf('%d', 0);                % train duration (zero for SPES)
param.StimulationConfigurations.Value{10,end-1}  = sprintf('%d', 0);                % train freq (zero for SPES)
param.StimulationConfigurations.ColumnLabels{end+1,1} = sprintf('CCEP Cathode');
param.StimulationConfigurations.Value{1, end}  = sprintf('%d', 0);                  %cathode bool
param.StimulationConfigurations.Value{2, end}  = sprintf('%d', 1);                  % n pulses 
param.StimulationConfigurations.Value{3, end}  = sprintf('%d', CCEP_amplitude);   % amp 1;
param.StimulationConfigurations.Value{4, end}  = sprintf('%d', CCEP_amplitude);   % amp 2;
param.StimulationConfigurations.Value{5, end}  = sprintf('%d', 250);                %pw 1;
param.StimulationConfigurations.Value{6, end}  = sprintf('%d', 250);                %pw 2;
param.StimulationConfigurations.Value{7, end}  = sprintf('%d', 50);                 %stim freq (does not matter for SPES) 
param.StimulationConfigurations.Value{8, end}  = sprintf('%d', 53);                 %inter pulse duration (53 us min) 
param.StimulationConfigurations.Value{9, end}  = sprintf('%d', 0);                  % train duration (zero for SPES) 
param.StimulationConfigurations.Value{10,end}  = sprintf('%d', 0);                  % train freq (zero for SPES) 
param.StimulationConfigurations.ColumnLabels{end+1,1} = sprintf('CCEP Anode');
param.StimulationConfigurations.RowLabels = rowLabs;
n_configs =size(stim_configs,1);



%% Set up Stimuli
n_jitters = n_configs * num_blocks;
stimRowLabs = {'caption';'icon';'av';'EarlyOffsetExpression';'StimulusDuration'};
param.Stimuli.Section      = 'Application';
param.Stimuli.Type         = 'matrix';
param.Stimuli.DefaultValue = '';
param.Stimuli.LowRange     = '';
param.Stimuli.HighRange    = '';
param.Stimuli.Comment      = 'captions and icons to be displayed, sounds to be played for different stimuli';
param.Stimuli.Value        = cell(5,n_configs+n_jitters+8); % 
% param.Stimuli.Value      = cell(4,2*n_configs+2); %if using empty blocks for isi
onset = 0;
% vidIDX = 1;
% v=VideoReader(video);
% [vidSlice,vidIDX,onset]=sliceVideo(stimuliDir,v,onset,duration,vidIDX);
param.Stimuli.Value{1,1}   = 'Press Space To Begin'; %stim code 1 = start screen
param.Stimuli.Value{2,1}   = '';
param.Stimuli.Value{3,1}   = '';
param.Stimuli.Value{4,1}   = 'Keydown==32';
param.Stimuli.Value{5,1}   = '40s';

videoStart_duration = 0.1; % seconds
param.Stimuli.Value{1,2}   = ''; % stim code 2 = begin video
param.Stimuli.Value{2,2}   = '';
param.Stimuli.Value{3,2}   = videoPath;
param.Stimuli.Value{4,2}   = '';
param.Stimuli.Value{5,2}   = sprintf('%ds',videoStart_duration);
videoStimCode = 2;
param.Stimuli.ColumnLabels = cell(n_configs+n_jitters+8,1);
param.Stimuli.ColumnLabels{1} = 'Start Screen';
param.Stimuli.ColumnLabels{2} = 'Begin Video';

% pause between blocks stimuli
loc = 3;
blockPause = 10; % seconds
blockPause_stimCode = loc;
param.Stimuli.Value{1,loc}   = ''; %Pause Between Blocks
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = '';
param.Stimuli.Value{4,loc}   = '';
param.Stimuli.Value{5,loc}   = '10s';
param.Stimuli.ColumnLabels{loc,1} = sprintf('Block Pause');

% CCEP stimuli 1
CCEP_duration = .2; % seconds
loc = loc+1;
CCEP_stimcode_1 = loc;
param.Stimuli.Value{1,loc}   = ''; %Trigger CCEP
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = '';
param.Stimuli.Value{4,loc}   = '';
param.Stimuli.Value{5,loc}   = '200ms';
param.Stimuli.ColumnLabels{loc,1} = sprintf('CCEP');

% CCEP stimuli 2
loc = loc+1;
CCEP_stimcode_2 = loc;
param.Stimuli.Value{1,loc}   = ''; %Trigger CCEP
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = '';
param.Stimuli.Value{4,loc}   = '';
param.Stimuli.Value{5,loc}   = '200ms';
param.Stimuli.ColumnLabels{loc,1} = sprintf('CCEP');

% CCEP Pause
loc = loc+1;
CCEP_pause_stimCode = loc;
param.Stimuli.Value{1,loc}   = ''; %Pause after leading CCEP
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = '';
param.Stimuli.Value{4,loc}   = '';
param.Stimuli.Value{5,loc}   = sprintf('%ds',CCEP_ISI);
param.Stimuli.ColumnLabels{loc,1} = sprintf('CCEP Pause');

%Reset Video
loc = loc+1;
vid_reset_duration = 0.5; % seconds
vidReset_stimCode = loc;
param.Stimuli.Value{1,loc}   = ''; 
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = videoPath; %Reset Video
param.Stimuli.Value{4,loc}   = '';
param.Stimuli.Value{5,loc}   = sprintf('500ms');
param.Stimuli.ColumnLabels{loc,1} = sprintf('Video Reset');


% end of run stimuli
loc = loc+1;
endRun_stimcode = loc;
param.Stimuli.Value{1,loc}   = 'End of Run';
param.Stimuli.Value{2,loc}   = '';
param.Stimuli.Value{3,loc}   = '';
param.Stimuli.Value{4,loc}   = 'Keydown==32';
param.Stimuli.Value{5,loc}   = '40s';
param.Stimuli.ColumnLabels{loc} = sprintf('End Run');

stimuli_duration = 1;
for idx=1:n_configs % stim codes 3:n_configs+2 are stim trials, 2 seconds long with the first second 
    loc=idx+8;
    param.Stimuli.Value{1,loc}   = ''; %caption
    param.Stimuli.Value{2,loc}   = ''; %icon
    param.Stimuli.Value{3,loc}   = ''; %av
    param.Stimuli.Value{4,loc}   = ''; %interrupt
    param.Stimuli.Value{5,loc}   = sprintf('%ds',stimuli_duration); %duration   
    param.Stimuli.ColumnLabels{loc,1} = sprintf('stim config %d',stimLabel(idx));
end

param.Stimuli.RowLabels    = stimRowLabs;

jitter_idx = zeros(n_jitters,1);
jitter_vals = zeros(n_jitters,1);
if allowCCEPs
    baseTime = 2;
else
    baseTime = 2 + 2*CCEP_ISI + 0.4; % 2s base plus both CCEP pauses + CCEPs durations (400 ms)
end
for j= loc+1:size(param.Stimuli.Value,2) % setup random jitters
    jitter = round((2*(rand)-1),1);
    param.Stimuli.Value{1,j}   = ''; %caption
    param.Stimuli.Value{2,j}   = ''; %icon
    param.Stimuli.Value{3,j}   = ''; %av
    param.Stimuli.Value{4,j}   = ''; %interrupt
    param.Stimuli.Value{5,j}   = sprintf('%ds',baseTime+jitter); %duration
    param.Stimuli.ColumnLabels{j,1} = sprintf('Jitter %d',j);
    jitter_idx(j-loc) = j;
    jitter_vals(j-loc) = baseTime+jitter;
end


res = checkStructFields(param.Stimuli);

%% sequence type
param.SequenceType.Section       = 'Application';
param.SequenceType.Type          = 'int';
param.SequenceType.DefaultValue  = '';
param.SequenceType.LowRange      = '';
param.SequenceType.HighRange     = '';
param.SequenceType.Comment       = 'Sequence of stimuli is 0 deterministic, 1 random, 2 P3Speller compatible (enumeration)';
param.SequenceType.Value         = {'0'};


%% set up Stim Triggers
TriggerExp = 'StimulusCode';
param.StimulationTriggers.Section = 'CereStim';
param.StimulationTriggers.Type = 'matrix';
param.StimulationTriggers.DefaultValue = '';
param.StimulationTriggers.LowRange = '';
param.StimulationTriggers.HighRange = '';
param.StimulationTriggers.Comment = 'Trigger Library for BLAES PARAM SWEEP';
for loc=1:numConfigs
    param.StimulationTriggers.Value{1,2*loc-1} = sprintf('%s == %d',TriggerExp,stimMat_cathode(loc,7)); %Expression
    param.StimulationTriggers.Value{2,2*loc-1} = sprintf('%d',stimMat_cathode(loc,8));%Config
    param.StimulationTriggers.Value{3,2*loc-1} = sprintf('%d',stimMat_cathode(loc,5));%Electrode
    param.StimulationTriggers.Value{1,2*loc}   = sprintf('%s == %d',TriggerExp,stimMat_anode(loc,7)); %Expression
    param.StimulationTriggers.Value{2,2*loc}   = sprintf('%d',stimMat_anode(loc,8));%Config
    param.StimulationTriggers.Value{3,2*loc}   = sprintf('%d',stimMat_anode(loc,5));%Electrode

    param.StimulationTriggers.ColumnLabels{2*loc-1,1} = sprintf('Trigger %d',2*loc-1);
    param.StimulationTriggers.ColumnLabels{2*loc,1}   = sprintf('Trigger %d',2*loc);
end
%CCEPs, need a trigger for each channel pair.
% channel pair 1
config = size(param.StimulationConfigurations.Value,2)-1;
param.StimulationTriggers.Value{1,end+1} = sprintf('%s == %d',TriggerExp,CCEP_stimcode_1); %Expression
param.StimulationTriggers.Value{2,end} = sprintf('%d',config);%Config
param.StimulationTriggers.Value{3,end} = sprintf('%d',cathodeChannels(1));%Electrode

config = size(param.StimulationConfigurations.Value,2);
param.StimulationTriggers.Value{1,end+1}   = sprintf('%s == %d',TriggerExp,CCEP_stimcode_1); %Expression
param.StimulationTriggers.Value{2,end}   = sprintf('%d',config);%Config
param.StimulationTriggers.Value{3,end}   = sprintf('%d',anodeChannels(1));%Electrode
param.StimulationTriggers.ColumnLabels{end+1,1} = 'CCEP 1 Cathode';
param.StimulationTriggers.ColumnLabels{end+1,1}   = 'CCEP 1 Anode';
% channel pair 2
config = size(param.StimulationConfigurations.Value,2)-1;
param.StimulationTriggers.Value{1,end+1} = sprintf('%s == %d',TriggerExp,CCEP_stimcode_2); %Expression
param.StimulationTriggers.Value{2,end} = sprintf('%d',config);%Config
param.StimulationTriggers.Value{3,end} = sprintf('%d',cathodeChannels(2));%Electrode

config = size(param.StimulationConfigurations.Value,2);
param.StimulationTriggers.Value{1,end+1}   = sprintf('%s == %d',TriggerExp,CCEP_stimcode_2); %Expression
param.StimulationTriggers.Value{2,end}   = sprintf('%d',config);%Config
param.StimulationTriggers.Value{3,end}   = sprintf('%d',anodeChannels(2));%Electrode
param.StimulationTriggers.ColumnLabels{end+1,1} = 'CCEP 2 Cathode';
param.StimulationTriggers.ColumnLabels{end+1,1}   = 'CCEP 2 Anode';

param.StimulationTriggers.RowLabels = {'Expression';'Config ID'; 'Electrode(s)'};





%% Setup Dynamic Mode
param.DynamicConfiguration.Section = 'CereStim';
param.DynamicConfiguration.Type = 'int';
param.DynamicConfiguration.DefaultValue = '0';
param.DynamicConfiguration.LowRange = '0';
param.DynamicConfiguration.HighRange = '1';
param.DynamicConfiguration.Comment= '';
param.DynamicConfiguration.Value = {'1'};

%% set up sequence and write to file
for i=1:length(videos)

parmName_dec = videos(i).name(1:20);
videoPath = fullfile(videos(i).folder,videos(i).name);
video = VideoReader(videoPath);
vid_duration = video.Duration; % duration of the video in seconds
clear video
param.Stimuli.Value{3,2}   = videoPath;
param.Stimuli.Value{3,vidReset_stimCode}   = videoPath; %Reset Video
param.Sequence.Section       = 'Application';
param.Sequence.Type          = 'intlist';
param.Sequence.DefaultValue  = '';
param.Sequence.LowRange      = '';
param.Sequence.HighRange     = '';
param.Sequence.Comment       = '';
param.Sequence.Value         = '';
param.Sequence.Value{1,1}    = '1'; % startup screen
param.Sequence.Value{2,1}    = '2'; % begin movie
idx = 3; % incrementer for trial sequence
jitterLocs = [];
experimentTime = videoStart_duration;
vidResetTracker = experimentTime;
numResets = 0;
for loc=1:size(trial_seq,1)
    % fprintf('%d exp time %d vid duration\n',experimentTime,vid_duration)

    if vidResetTracker > vid_duration
        numResets = numResets +1
        param.Sequence.Value{idx,1}  = sprintf('%d',vidReset_stimCode) % restart video
        vidResetTracker = 0;
        idx = idx +1;
    end


    channelPair = find(cathodeChannels== trial_seq(loc,5));
    if channelPair == 1
        CCEP_stimcode = CCEP_stimcode_1;
    elseif channelPair == 2
        CCEP_stimcode = CCEP_stimcode_2;
    end

    if allowCCEPs
        param.Sequence.Value{idx,1}  = sprintf('%d',CCEP_stimcode); % leading CCEP
        vidResetTracker = vidResetTracker + CCEP_duration;
        experimentTime = experimentTime + CCEP_duration;
        idx = idx +1;
        param.Sequence.Value{idx,1}  = sprintf('%d',CCEP_pause_stimCode); % CCEP Pause
        vidResetTracker = vidResetTracker + CCEP_ISI;
        experimentTime = experimentTime + CCEP_ISI;
        idx = idx +1;
    end

    param.Sequence.Value{idx,1}  = sprintf('%d',trial_seq(loc,7)); % Stimulation
    vidResetTracker = vidResetTracker + stimuli_duration;
    experimentTime = experimentTime + stimuli_duration;
    idx = idx +1;

    if allowCCEPs
        param.Sequence.Value{idx,1}  = sprintf('%d',CCEP_pause_stimCode); % CCEP Pause
        idx = idx +1;
        vidResetTracker = vidResetTracker + CCEP_ISI;
        experimentTime = experimentTime + CCEP_ISI;

        param.Sequence.Value{idx,1}  = sprintf('%d',CCEP_stimcode); % lagging CCEP
        idx = idx +1;
        vidResetTracker = vidResetTracker + CCEP_duration;
        experimentTime = experimentTime + CCEP_duration;
    end

    param.Sequence.Value{idx,1}  = sprintf('%d',jitter_idx(loc)); % jittered ISI
    jitterLocs = [jitterLocs; idx];
    vidResetTracker = vidResetTracker + jitter_vals(loc);
    experimentTime = experimentTime + jitter_vals(loc);
    idx = idx+1;

    if vidResetTracker > vid_duration
        experimentTime = experimentTime + videoStart_duration
        numResets = numResets +1
        param.Sequence.Value{idx,1}  = sprintf('%d',videoStimCode) % restart video
        vidResetTracker = 0;
        idx = idx +1;
    end


    if mod(loc,n_configs)==0 && loc~=size(trial_seq,1) % block pauses
        param.Sequence.Value{idx,1}  = sprintf('%d',blockPause_stimCode);
        vidResetTracker = vidResetTracker + blockPause;
        experimentTime = experimentTime + blockPause;
        idx = idx +1;
    end

    if vidResetTracker > vid_duration
        numResets = numResets +1
        param.Sequence.Value{idx,1}  = sprintf('%d',vidReset_stimCode) % restart video
        vidResetTracker = 0;
        idx = idx +1;
    end
    
end
param.Sequence.Value{end+1,1} = sprintf('%d',endRun_stimcode); % end
fprintf('%d exp time %d vid duration',experimentTime,vid_duration)

% Write Trial Sweep Parm
fname = sprintf('BLAES_stimsweep_run_%s.prm',parmName_dec);
filename = fullfile(parmDir,fname);
parameter_lines = convert_bciprm( param );
fid = fopen(filename, 'w');

for loc=1:length(parameter_lines)
    fprintf( fid, '%s', parameter_lines{loc} );
    fprintf( fid, '\r\n' );
end
fclose(fid);
fprintf('\nwrote %s to %s\n', fname,parmDir);
end
%%
%
%
%
%
%
%
%
%
%
%
%
%
%% Stimulation Testing Parm Generation
if numel(conditions2remove) == 0 && generateTest 
    % set up Stim Configs
    stim_configs = getUniqueConfigs(stimMat_cathode,stimMat_anode);
    testing_param.StimulationConfigurations.Section = 'CereStim';
    testing_param.StimulationConfigurations.Type = 'matrix';
    testing_param.StimulationConfigurations.DefaultValue = '';
    testing_param.StimulationConfigurations.LowRange = '';
    testing_param.StimulationConfigurations.HighRange = '';
    testing_param.StimulationConfigurations.Comment = 'Configurations for BLAES PARAM SWEEP';
    testing_param.StimulationConfigurations.Value = cell(length(rowLabs),size(stim_configs,1));
    
    for loc=1:size(stim_configs,1)
        testing_param.StimulationConfigurations.Value{1, loc}  = sprintf('%d',stim_configs(loc,1));
        testing_param.StimulationConfigurations.Value{2, loc}  = sprintf('%d',stim_configs(loc,2));
        testing_param.StimulationConfigurations.Value{3, loc}  = sprintf('%d',stim_configs(loc,4));
        testing_param.StimulationConfigurations.Value{4, loc}  = sprintf('%d',stim_configs(loc,4));
        testing_param.StimulationConfigurations.Value{5, loc}  = sprintf('%d',stim_configs(loc,3));
        testing_param.StimulationConfigurations.Value{6, loc}  = sprintf('%d',stim_configs(loc,3));
        testing_param.StimulationConfigurations.Value{7, loc}  = sprintf('%d',stim_configs(loc,6));
        testing_param.StimulationConfigurations.Value{8, loc}  = sprintf('%d',53);
        testing_param.StimulationConfigurations.Value{9, loc}  = sprintf('%d',1);
        testing_param.StimulationConfigurations.Value{10,loc}  = sprintf('%d',8);
        testing_param.StimulationConfigurations.ColumnLabels{loc,1} = sprintf('Configuration %d',loc);
    end
    testing_param.StimulationConfigurations.RowLabels = rowLabs;
    
    
    
    n_configs =size(stim_configs,1);
    %% Set up Stimuli
    cerestim_map_path = strcat(stimuliDir,'\cerestim_map.png');
    stimRowLabs = {'caption';'icon';'av';'EarlyOffsetExpression';'StimulusDuration'};
    testing_param.Stimuli.Section      = 'Application';
    testing_param.Stimuli.Type         = 'matrix';
    testing_param.Stimuli.DefaultValue = '';
    testing_param.Stimuli.LowRange     = '';
    testing_param.Stimuli.HighRange    = '';
    testing_param.Stimuli.Comment      = 'captions and icons to be displayed, sounds to be played for different stimuli';
    
    testing_param.Stimuli.Value{1,1}   = '';
    testing_param.Stimuli.Value{2,1}   = cerestim_map_path;
    testing_param.Stimuli.Value{3,1}   = '';
    testing_param.Stimuli.Value{4,1}   = 'Keydown==32';
    testing_param.Stimuli.Value{5,1}   = '600s';
    
    
    testing_param.Stimuli.ColumnLabels{1} = '1';
    
    testing_param.Stimuli.RowLabels    = stimRowLabs;
    
    res = checkStructFields(testing_param.Stimuli);
    %% set up sequence;
    testing_param.Sequence.Section       = 'Application';
    testing_param.Sequence.Type          = 'intlist';
    testing_param.Sequence.DefaultValue  = '';
    testing_param.Sequence.LowRange      = '';
    testing_param.Sequence.HighRange     = '';
    testing_param.Sequence.Comment       = '';
    A = ones(100,1);
    B = arrayfun(@num2str, A, 'UniformOutput', 0); % make a cell array of 100 ones to ensure experiment does not timeout. 
    testing_param.Sequence.Value         = B; % startup screen
    
    
    %% sequence type
    testing_param.SequenceType.Section       = 'Application';
    testing_param.SequenceType.Type          = 'int';
    testing_param.SequenceType.DefaultValue  = '';
    testing_param.SequenceType.LowRange      = '';
    testing_param.SequenceType.HighRange     = '';
    testing_param.SequenceType.Comment       = 'Sequence of stimuli is 0 deterministic, 1 random, 2 P3Speller compatible (enumeration)';
    testing_param.SequenceType.Value         = {'0'};
    
    
    %% set up Stim Triggers
    TriggerExp = 'Keydown';
    testing_param.StimulationTriggers.Section = 'CereStim';
    testing_param.StimulationTriggers.Type = 'matrix';
    testing_param.StimulationTriggers.DefaultValue = '';
    testing_param.StimulationTriggers.LowRange = '';
    testing_param.StimulationTriggers.HighRange = '';
    testing_param.StimulationTriggers.Comment = 'Trigger Library for BLAES PARAM SWEEP';
    stimDescription = num2cell(stimMat_cathode);
    stimDescription_labels = stimMat_cols;
    stimDescription_labels{end+1} = 'key';
    stimDescription_labels{end+1} = 'keypress_value';
    fig = figure(1);
    set(gcf, 'Position', [1 1 2560 1080])
    sgtitle(sprintf('Load Configuration with Indicated Key\nPress ENTER to Deliver Stimulation'),'FontSize',18)
    chargeDensity = stimMat_cathode(:,3).*stimMat_cathode(:,4).*stimMat_cathode(:,2)*carrier_freq*stimuli_duration * 1e-6 / (electrodeSurfaceArea*0.01); % us * uA * numPulses * Hz * s/cm^2 = pC/cm^2 -> pC/cm^2 * 1e-6 = uC/cm^2
    chargeDensityPerPhase = stimMat_cathode(:,3).*stimMat_cathode(:,4)* 1e-6 / (electrodeSurfaceArea*0.01); % us * uA/cm^2= pC/cm^2-> pC/cm^2 *1e-6 = uC/cm^2
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5386002/
    aggMat_cath = [stimMat_cathode stimLabel chargeDensity chargeDensityPerPhase];
    [sort_cathode, sortIdx] = sortrows(aggMat_cath, [size(aggMat_cath,2)-1,4,size(aggMat_cath,2)-2],'ascend');
    sort_anode = stimMat_anode(sortIdx,:);
    sort_stimLabel = stimLabel(sortIdx);
    sort_keyMap = keyMap';
    sort_keyMap = sort_keyMap(sortIdx,:);
    for loc=1:numConfigs
        subplot(3,numConfigs/3,loc)
        f_g = sort_cathode(loc,6);
        cath_g = sort_cathode(loc,5);
        anode_g = sort_anode(loc,5);
        amp_g = sort_cathode(loc,4);
        pw_g = sort_cathode(loc,3);
        np_g = sort_cathode(loc,2);
        charge_g = sort_cathode(loc,end-1);
        chargePerPhase_g = sort_cathode(loc,end);
        [t,y] = generate_theta_burst_waveform(f_g,amp_g);
        plot(1000*t,y,'Linewidth',2,'Color','r')
        xlim([0,360])
        ylim([-0.5,max(stimAmps)])
        ylabel('amp (uA)')
        xlabel('time (ms)')
        titleblock = sprintf('Combo #%d, Key %s pair:[%d, %d]\n amp: %duA, f:%dHz pw: %dus\n%d uC/cm^2 per 1s, \n%d uC/cm^2 per phase',sort_stimLabel(loc),sort_keyMap{loc,2},cath_g,anode_g,amp_g,f_g,pw_g,charge_g,chargePerPhase_g);
        title(titleblock)
        
    end
    fig2 = figure(2);
    % % set(gcf, 'Position', get(0, 'Screensize')); %fullscreen fig generation
    set(gcf, 'Position', [1 1 1980 1080])
    for loc=1:numConfigs
        testing_param.StimulationTriggers.Value{1,2*loc-1} = sprintf('%s == %d',TriggerExp,keyMap{2,loc}); %Expression
        testing_param.StimulationTriggers.Value{2,2*loc-1} = sprintf('%d',stimMat_cathode(loc,8));%Config
        testing_param.StimulationTriggers.Value{3,2*loc-1} = sprintf('%d',stimMat_cathode(loc,5));%Electrode
        subplot(3,numConfigs/3,loc)
        f_g = stimMat_cathode(loc,6);
        cath_g = stimMat_cathode(loc,5);
        anode_g = stimMat_anode(loc,5);
        amp_g = stimMat_cathode(loc,4);
        pw_g = stimMat_cathode(loc,3);
        np_g = stimMat_cathode(loc,2);
        [t,y] = generate_theta_burst_waveform(f_g,amp_g);
        plot(1000*t,y,'Linewidth',2,'Color','r')
        xlim([0,360])
        ylim([-0.5,max(stimAmps)])
        ylabel('amp (uA)')
        xlabel('time (ms)')
        titleblock = sprintf('Combination #%d, Key %s \n [Cathode, Anode]: [%d, %d]\namp: %d mA, f: %d Hz\npw: %d us, np: %d',stimLabel(loc),keyMap{2,loc},cath_g,anode_g,amp_g,f_g,pw_g,np_g);
        title(titleblock)
        stimDescription{loc,9} = keyMap{1,loc};
        stimDescription{loc,10} = keyMap{2,loc};
        testing_param.StimulationTriggers.Value{1,2*loc}   = sprintf('%s == %d',TriggerExp,keyMap{2,loc}); %Expression
        testing_param.StimulationTriggers.Value{2,2*loc}   = sprintf('%d',stimMat_anode(loc,8));%Config
        testing_param.StimulationTriggers.Value{3,2*loc}   = sprintf('%d',stimMat_anode(loc,5));%Electrode
    
        testing_param.StimulationTriggers.ColumnLabels{2*loc-1,1} = sprintf('Trigger %d',2*loc-1);
        testing_param.StimulationTriggers.ColumnLabels{2*loc,1}   = sprintf('Trigger %d',2*loc);
    end
    testing_param.StimulationTriggers.RowLabels = {'Expression';'Config ID'; 'Electrode(s)'};
    locDescriptions = cell2struct(stimDescription,stimDescription_labels,2);
    % saveas(1,cerestim_map_path);
    exportgraphics(fig,cerestim_map_path,'Resolution',300);
    %% Setup Dynamic Mode
    testing_param.DynamicConfiguration.Section = 'CereStim';
    testing_param.DynamicConfiguration.Type = 'int';
    testing_param.DynamicConfiguration.DefaultValue = '0';
    testing_param.DynamicConfiguration.LowRange = '0';
    testing_param.DynamicConfiguration.HighRange = '1';
    testing_param.DynamicConfiguration.Comment= '';
    testing_param.DynamicConfiguration.Value = {'1'};
    
    testing_param.StartExpression.Section = 'CereStim';
    testing_param.StartExpression.Type = 'string';
    testing_param.StartExpression.DefaultValue = '';
    testing_param.StartExpression.LowRange = '';
    testing_param.StartExpression.HighRange = '';
    testing_param.StartExpression.Comment= 'BCI2000 expression to start stim';
    testing_param.StartExpression.Value = {'Keydown==13'};
    
    
    %% Write Trial Test Parm
    fname = 'BLAES_stimsweep_test.prm';
    filename = fullfile(parmDir,fname);
    parameter_lines = convert_bciprm( testing_param );
    fid = fopen(filename, 'w');
    
    for loc=1:length(parameter_lines)
        fprintf( fid, '%s', parameter_lines{loc} );
        fprintf( fid, '\r\n' );
    end
    fclose(fid);
else
    fprintf('Testing param not generated\nGenerate Test: %d\nNum Excluded Conditions: %d',generateTest,numel(conditions2remove))
end

%% Functions
function X = pseduoRandReplace(A,multi,amp_init,channel_init,ampCol,channelCol)
% takes all stimulation conditions and shuffles them pseudorandomly based
% the following constraints
% constraint 1: 2's in column 2 (stimAmp) are not allowed to be repeated
% constraint 2: no repeats in column 3 (stim site)
[X,rows] = initPermutation(A,multi,amp_init,channel_init,ampCol,channelCol);

temp = A(rows(2:end),:);
rowsAdded = [rows(1)];
unused = [rows(2:end)]';
fail = 0;
for j=2:length(rows)
    % stim Amp Condition
    cond1 = [];
    checkVal = X(end,ampCol);
    if checkVal == min(X(:,ampCol))
        cond1 = unused;
    else
        for i=1:length(unused)
            if temp(i,ampCol) ~= checkVal
                cond1 = [cond1 unused(i)];
            end
        end
    end
    % Channel Condition
    cond2 = [];
    singularityFlag = numel(unique(X(:,channelCol)));
    checkVal = X(end,channelCol);
    for i=1:length(unused)
        if temp(i,channelCol) ~= checkVal || singularityFlag ==1 
            cond2 = [cond2 unused(i)];
        end
    end
    targets = intersect(cond1, cond2);
    potentialRows = targets(randperm(length(targets)));
    if isempty(potentialRows)
        fail = 1;
        break
    end
    rowInsert = potentialRows(1);
    loc = find(unused == rowInsert);
    unused(loc) = [];
    rowsAdded = [rowsAdded; loc];
    X = [X; temp(loc,:)];
    temp(loc,:) = [];
end
if fail
    X = pseduoRandReplace(A,multi,amp_init,channel_init,ampCol,channelCol);
end
end

function [X,rows] = initPermutation(A,multi,amp_init,channel_init,ampCol,channelCol)
rows = randperm(size(A,1));
X = A(rows(1),:);
if multi
    v2 = X(ampCol);
    v3 = X(channelCol);
    if v2 == amp_init || v3 == channel_init
        [X,rows] = initPermutation(A,multi,amp_init,channel_init,ampCol,channelCol);
    end
end
end

function flag = checkRepeats(A,col,~)
% returns true if repeats are detected, false otherwise
col_target = A(:, col);
flag = any(col_target(1:end-1) == col_target(2:end));
end

function flag = checkRepeatVal(A,col,val)
% returns true if the specified value is repeated, false otherwise
col_target = A(:, col);
flag = any(diff(col_target) == 0 & col_target(1:end-1) == val);
end

function configs = getUniqueConfigs(cathode,anode)
joint = [cathode; anode];
configs = [];
ID = 0;
for i=1:size(joint,1)
    val = joint(i,8);
    if val ~= ID
        configs = [configs; [joint(i,1:6),val]];
        ID = val;
    end
end
end
function isAllStringsOrCells = checkStructFields(structure)
    % Initialize the flag as true
    isAllStringsOrCells = true;
    
    % Get field names of the structure
    fields = fieldnames(structure);
    
    % Iterate over each field
    for i = 1:length(fields)
        fieldValue = structure.(fields{i});
        
        % Check if the field is a structure itself and recurse
        if isstruct(fieldValue)
            isAllStringsOrCells = checkStructFields(fieldValue);
            if ~isAllStringsOrCells
                return; % Early return if any subfield is not a string or cell
            end
        % Check if the field is a string or cell
        elseif ~(isstring(fieldValue) || iscell(fieldValue)|| ischar(fieldValue))
            isAllStringsOrCells = false;
            return; % Early return if the field is not a string or cell
        end
    end
end

% function bci2ktools(BCI2KPath)
% olddir = pwd;
% cd(BCI2KPath);
% cd tools, cd matlab;
% bci2000path -AddToMatlabPath tools/matlab;
% bci2000path -AddToMatlabPath tools/mex;
% bci2000path -AddToSystemPath tools/cmdline;   % required so that BCI2000CHAIN can call the command-line tools
% cd(olddir); % change directory back to where we were before
% clear olddir;
% end

function x = rm_combinations(input, idx2remove,redoIndexing)

x = input;
if redoIndexing
stimcodeStart = min(x(:,7));
stimConfigStart = min(x(:,8));
end
if numel(idx2remove) > 0
    x(idx2remove,:) = [];
end
if redoIndexing
    newStimCodes = stimcodeStart:1:stimcodeStart+size(x,1)-1;
    x(:,7) = newStimCodes;
    newConfigs = zeros(size(x,1),1);
    num = stimConfigStart;
    for i=1:size(x,1)/2
        newConfigs(2*i-1) = num;
        newConfigs(2*i) = num;
        num = num+1;
    end
    x(:,8) = newConfigs;
end
end
