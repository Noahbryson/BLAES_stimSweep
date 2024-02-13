%% BCI2000 BLAES Stim Parameter Sweep .prm Generator
%
% Change the n_rows and n_stimuli variables to store more information with
% the stimuli or add additional stimuli. Best practice is to separate
% stimuli into banks (e.g. 1-25, 101-125, etc) for easy evaluation later.
%
% Note that every stimulus needs to have an index for every row desired,
% even if that row label is not meaningful for the stimulus.
%
% A sequence is created to alternate the fixation cross stimuli with the
% image stimuli.
%
% The stimuli and meaningful parameters are written into a param
% variable and stored as a *.prm file using the convert_bciprm function.

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


%% User Parameters
% Set the path of the BCI2000 main directory here
loadBCI2kTools;
BCI2KPath = 'C:\BCI2000\BCI2000';
load('C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\stimConfigRowLabs.mat')
videoPath = '';
% enter stim channels, different at WU vs UU as WU plugs in hardware
% manually due to number of recording channels.

stimChannels = [1,3]; %stim channels passed to stimulator (cerestim?)
stimAmps = [1,2]; %mA
pulseWidth = [250, 500]; %us, single phase, double for entire biphasic stim
frequencies = [28, 50, 80]; %Hz
numTrials = 10;
%% Generate Potential Conditions
% note for evetual configurations, only half as many configs need to be
% loaded since configs are repeated at two locations. triggers must be set
% up to be unique, however
stimMat_cathode = [];
cols = {'cathodeFirst' 'numPulses' 'pulseWidth' 'amp' 'channel' 'freq','stimCode','config ID'};
% Channel Condition
channelCol = 5;
% stim Amp Condition
ampCol = 4;
stimCode = 3; % stimCode 1 reserved for AV start sequence
configID = 1;
for i=1:length(pulseWidth)
    for j=1:length(stimAmps)
        for k=1:length(frequencies)
            for q=1:length(stimChannels)

                if pulseWidth(i) > 250 && stimAmps(j) > 1

                else
                    temp = [1,4,pulseWidth(i),stimAmps(j),stimChannels(q),frequencies(k),stimCode,configID];
                    stimCode = stimCode+1;
                    stimMat_cathode = [stimMat_cathode; temp];
                    if q==2
                        configID = configID+ 1;
                    end
                end
            end
        end
    end
end
stimMat_anode = stimMat_cathode;
stimMat_anode(:,1) = stimMat_anode(:,1) - 1; %switch to anodic leading
stimMat_anode(:,end) = stimMat_anode(:,end) + configID-1; %adjust configIDs
stimMat_anode(:,channelCol) = stimMat_anode(:,channelCol) + 1; % adjust channel to make bipolar pair
numConfigs = size(stimMat_cathode,1);
numParams = size(stimMat_cathode,2);
seq = pseduoRandReplace(stimMat_cathode,0,0,0,ampCol,channelCol)';
checkRepeats(seq,channelCol) % channel column
checkRepeatVal(seq,ampCol,2) % amplitude column
%% Generate Trial Sequence
trial_seq = zeros(numConfigs*numTrials,numParams);
for i=1:numTrials
    stop = numConfigs*i;
    start = stop-(numConfigs-1);
    if i == 1
        multi = 0;
        amp_init = 0;
        chan_init = 0;
    else
        multi = 1;
        amp_init = trial_seq(stop-numConfigs,ampCol);
        chan_init = trial_seq(stop-numConfigs,channelCol);
    end

    trial_seq(start:stop,:) = pseduoRandReplace(stimMat_cathode,multi,amp_init,chan_init,ampCol,channelCol);
end
checkRepeats(trial_seq,channelCol) % channel column
checkRepeatVal(trial_seq,ampCol,2) % amplitude column

%% Parm Generation
% set up Stim Configs
stim_configs = getUniqueConfigs(stimMat_cathode,stimMat_anode);
param.StimulationConfigurations.Section = 'CereStim';
param.StimulationConfigurations.Type = 'matrix';
param.StimulationConfigurations.DefaultValue = '';
param.StimulationConfigurations.LowRange = '';
param.StimulationConfigurations.HighRange = '';
param.StimulationConfigurations.Comment = 'Configurations for BLAES PARAM SWEEP';
param.StimulationConfigurations.Value = cell(length(rowLabs),size(stim_configs,1));

for i=1:size(stim_configs,1)
    param.StimulationConfigurations.Value{1, i}  = sprintf('%d',stim_configs(i,1));
    param.StimulationConfigurations.Value{2, i}  = sprintf('%d',stim_configs(i,2));
    param.StimulationConfigurations.Value{3, i}  = sprintf('%d',stim_configs(i,4));
    param.StimulationConfigurations.Value{4, i}  = sprintf('%d',stim_configs(i,4));
    param.StimulationConfigurations.Value{5, i}  = sprintf('%d',stim_configs(i,3));
    param.StimulationConfigurations.Value{6, i}  = sprintf('%d',stim_configs(i,3));
    param.StimulationConfigurations.Value{7, i}  = sprintf('%d',stim_configs(i,6));
    param.StimulationConfigurations.Value{8, i}  = sprintf('%d',53);
    param.StimulationConfigurations.Value{9, i}  = sprintf('%d',8);
    param.StimulationConfigurations.Value{10,i}  = sprintf('%d',1);
    param.StimulationConfigurations.ColumnLabels{i,1} = sprintf('Configuration %d',i);
end
param.StimulationConfigurations.RowLabels = rowLabs;



n_configs =size(stim_configs,1);
% Set up Stimuli
stimRowLabs = {'caption';'icon';'av';'EarlyOffsetExpression';'StimulusDuration'};
param.Stimuli.Section      = 'Application';
param.Stimuli.Type         = 'matrix';
param.Stimuli.DefaultValue = '';
param.Stimuli.LowRange     = '';
param.Stimuli.HighRange    = '';
param.Stimuli.Comment      = 'captions and icons to be displayed, sounds to be played for different stimuli';
param.Stimuli.Value        = cell(4,n_configs+2); % if isi built into stim blocks
% param.Stimuli.Value      = cell(4,2*n_configs+2); %if using empty blocks for isi
param.Stimuli.Value(:,1)   = {'Press Space To Begin';'';'';'KeyDown==32';'40s'};
param.Stimuli.Value(:,2)   = {'';'';videoPath;'';'5s'};
for i=3:n_configs+2
    param.Stimuli.Value(1,i)   = {''}; %caption
    param.Stimuli.Value(2,i)   = {''}; %icon
    param.Stimuli.Value(3,i)   = {''}; %av
    param.Stimuli.Value(4,i)   = {''}; %interrupt
    param.Stimuli.Value(5,i)   = {'1.5'}; %duration
end
param.Stimuli.Value(:,end+1)   = {'End of Run';'';'';'KeyDown==32';'40s'};
% for i=n_configs+3:2*n_configs+2  %if using empty blocks for isi
% param.Stimuli.Value(1,i)        = {''};
% param.Stimuli.Value(2,i)        = {''};
% param.Stimuli.Value(3,i)        = {''};
% param.Stimuli.Value(4,i)        = {'2'}; %isi
% end

param.Stimuli.RowLabels    = stimRowLabs;
param.Stimuli.ColumnLabels    = cell(n_configs+2,1);

% set up sequence;
param.Sequence.Section       = 'Application';
param.Sequence.Type          = 'intlist';
param.Sequence.DefaultValue  = '';
param.Sequence.LowRange      = '';
param.Sequence.HighRange     = '';
param.Sequence.Comment       = '';
param.Sequence.Value         = '';
param.Sequence.Value{1,1} = '1'; % startup screen
param.Sequence.Value{2,1} = '2'; % begin movie
for i=1:size(trial_seq,1)
    param.Sequence.Value{i+2,1}  = sprintf('%d',trial_seq(i,7));
    param.Sequence.NumericValue(i+2,1)  = trial_seq(i,7);
end
param.Sequence.Value{end+1,1} = sprintf('%d',max(trial_seq(:,7))+3); % end

% sequence type
param.SequenceType.Section       = 'Application';
param.SequenceType.Type          = 'int';
param.SequenceType.DefaultValue  = '';
param.SequenceType.LowRange      = '';
param.SequenceType.HighRange     = '';
param.SequenceType.Comment       = 'Sequence of stimuli is 0 deterministic, 1 random, 2 P3Speller compatible (enumeration)';
param.SequenceType.Value         = '0';
param.SequenceType.NumericValue  =  0;


% set up Stim Triggers
TriggerExp = 'StimulusCode';
param.StimulationTriggers.Section = 'CereStim';
param.StimulationTriggers.Type = 'matrix';
param.StimulationTriggers.DefaultValue = '';
param.StimulationTriggers.LowRange = '';
param.StimulationTriggers.HighRange = '';
param.StimulationTriggers.Comment = 'Trigger Library for BLAES PARAM SWEEP';
for i=1:numConfigs
    param.StimulationTriggers.Value{1,2*i-1} = sprintf('%s == %d',TriggerExp,stimMat_cathode(i,7)); %Expression
    param.StimulationTriggers.Value{2,2*i-1} = sprintf('%d',stimMat_cathode(i,8));%Config
    param.StimulationTriggers.Value{3,2*i-1} = sprintf('%d',stimMat_cathode(i,5));%Electrode
    param.StimulationTriggers.Value{1,2*i}   = sprintf('%s == %d',TriggerExp,stimMat_anode(i,7)); %Expression
    param.StimulationTriggers.Value{2,2*i}   = sprintf('%d',stimMat_anode(i,8));%Config
    param.StimulationTriggers.Value{3,2*i}   = sprintf('%d',stimMat_anode(i,5));%Electrode

    param.StimulationTriggers.ColumnLabels{2*i-1} = sprintf('Trigger %d',2*i-1);
    param.StimulationTriggers.ColumnLabels{2*i}   = sprintf('Trigger %d',2*i);
end
param.StimulationTriggers.RowLabels = {'Expression';'Config ID'; 'Electrode(s)'};













%% Write Parm
% filename = "C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\parms\BLAES_stimsweep_run.prm";
% parameter_lines = convert_bciprm( param );
% fid = fopen(filename, 'w');
%
% for i=1:length(parameter_lines)
%     fprintf( fid, '%s', parameter_lines{i} );
%     fprintf( fid, '\r\n' );
% end
% fclose(fid);










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
    if checkVal == 1
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
    checkVal = X(end,channelCol);
    for i=1:length(unused)
        if temp(i,channelCol) ~= checkVal
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
    val = joint(i,end);
    if val ~= ID
        configs = [configs; [joint(i,1:end-2),val]];
        ID = val;
    end
end
end
