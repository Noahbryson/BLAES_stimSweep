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
BCI2KPath = 'C:\BCI2000\BCI2000'; % Set the path of the BCI2000 main directory here
bci2ktools(BCI2KPath) % sets up access to BCI2000 functions and .mex files, just change local dir above



% enter stim channels, different at WU vs UU as WU plugs in hardware
% manually due to number of recording channels.
cathodeChannels = [1,3]; % cathode leading stim channels passed to stimulator 
anodeChannels   = [2,4]; % anode leading stim channels passed to stimulator
    % be sure cathode and anode channels are aligned with each other for
    % bipolar pairs. 
stimAmps = [1,2]; %mA, note 2 mA is not paired with 500 us pulse width.
pulseWidth = [250, 500]; %us, single phase, double for entire biphasic stim
frequencies = [33, 50, 80]; %Hz
numTrials = 10;

conditions2remove = [];% add the numeric value as they appear on the testing image

%% Pathing
disp('select stimuli main folder')
% stimuliRoot = uigetdir();
disp('select saving root folder (i.e. C:\Paradigms')
% saveRoot = uigetdir();
root = 'C:\Paradigms';
stimuliDir = fullfile(root,'\tasks\BLAES\BLAES_param_sweep\stimuli');
checkDir(stimuliDir);
parmDir  = fullfile(root,'\parms\BLAES\_BLAES_param_sweep');
checkDir(parmDir);


videos = dir(strcat(stimuliDir,'\','*.mp4'));
video_num = 3;
parmName_dec = videos(video_num).name(1:20);
video = fullfile(videos(video_num).folder,videos(video_num).name);




%% Keyboard Map
keyboard = struct;
temp = cellfun(@upper,{'q' 'w' 'e' 'r' 't' 'y' 'u' 'i' 'o' 'p' 'a' 's' 'd' 'f' 'g' 'h' 'j' 'k' 'l' 'z' 'x' 'c' 'v' 'b' 'n' 'm'},'UniformOutput', false);
keyboard.keys= temp;
keyboard.X = [0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7]*10;
keyboard.Y = [2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]*20;
load("C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\keys2numbers.mat")
test = cellfun(@num2str, keyboard.keys, 'UniformOutput', false);
test2 = cellfun(@num2str, fullKeyMap(1,:), 'UniformOutput', false);
logicalIdx = ismember(test2,test);
keyMap = fullKeyMap(:,logicalIdx);
%% Generate Potential Conditions
% note for evetual configurations, only half as many configs need to be
% loaded since configs are repeated at two locations. triggers must be set
% up to be unique, however
stimMat_cathode = [];
stimLabel = [];
stimMat_cols = {'cathodeFirst' 'numPulses' 'pulseWidth' 'amp' 'channel' 'freq','stimCode','config_ID'};
% Channel Condition
channelCol = 5;
% stim Amp Condition
ampCol = 4;
stimCode = 3; % stimCode 1 reserved for AV start sequence
configID = 1;
am_freq = 8;
counter =0;
for i=1:length(pulseWidth)
    for j=1:length(stimAmps)
        for k=1:length(frequencies)
            for q=1:length(cathodeChannels)

                if pulseWidth(i) > 250 && stimAmps(j) >= max(stimAmps)+1
                    % do not pair 2 mA with 500 us PW
                else
                    temp = [1,floor(frequencies(k)/(2*am_freq))+1,pulseWidth(i),stimAmps(j),cathodeChannels(q),frequencies(k),stimCode,configID];
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
stimMat_cathode = rm_combinations(stimMat_cathode,conditions2remove);
stimLabel = rm_combinations(stimLabel,conditions2remove);
stimMat_anode = stimMat_cathode;
stimMat_anode(:,1) = stimMat_anode(:,1) - 1; %switch to anodic leading
stimMat_anode(:,8) = stimMat_anode(:,8) + configID-1; %adjust configIDs


% adjust channel to make bipolar pair
for i=1:size(stimMat_anode,1)
    idx = find(cathodeChannels == stimMat_anode(i,channelCol));
    stimMat_anode(i,channelCol) = anodeChannels(idx);
end
numConfigs = size(stimMat_cathode,1);
numParams = size(stimMat_cathode,2);
seq = pseduoRandReplace(stimMat_cathode,0,0,0,ampCol,channelCol)';
checkRepeats(seq,channelCol); % channel column, ensures no channels are repeated.
checkRepeatVal(seq,ampCol,max(stimAmps)); % amplitude column, ensures 2 mA is not repeated. 
%% Generate Trial Sequence
start = 1;
stop = numConfigs;
trial_seq = zeros(numConfigs*numTrials + (numTrials-2),numParams);
for i=1:numTrials
    
    if i == 1
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
    if i < numTrials
        trial_seq(stop+1,:) = [0 0 0 0 0 0 stimCode 0];
    end
    start = stop+2;
    stop = start+numConfigs-1;
    
    
end
checkRepeats(trial_seq,channelCol); % channel column
checkRepeatVal(trial_seq,ampCol,max(stimAmps)); % amplitude column

%% Trial Sweep Parm Generation

% set up Stim Configs
rowLabs = {'Cathode First';'Number of pulses';'Phase 1 amp (uA)';'Phase 2 amp (uA)';'Phase 1 duration (us)';'Phase 2 duration (us)';'Frequency (Hz)';'Interphase duration (us)';'Train Duration (s)';'Train Frequency (Hz)'};
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
    param.StimulationConfigurations.Value{9, i}  = sprintf('%d',1);
    param.StimulationConfigurations.Value{10,i}  = sprintf('%d',8);
    param.StimulationConfigurations.ColumnLabels{i,1} = sprintf('Configuration %d',i);
end
param.StimulationConfigurations.RowLabels = rowLabs;



n_configs =size(stim_configs,1);



%% Set up Stimuli
stimRowLabs = {'caption';'icon';'av';'EarlyOffsetExpression';'StimulusDuration'};
param.Stimuli.Section      = 'Application';
param.Stimuli.Type         = 'matrix';
param.Stimuli.DefaultValue = '';
param.Stimuli.LowRange     = '';
param.Stimuli.HighRange    = '';
param.Stimuli.Comment      = 'captions and icons to be displayed, sounds to be played for different stimuli';
param.Stimuli.Value        = cell(5,n_configs+4); % if isi built into stim blocks
% param.Stimuli.Value      = cell(4,2*n_configs+2); %if using empty blocks for isi

param.Stimuli.Value{1,1}   = 'Press Space To Begin';
param.Stimuli.Value{2,1}   = '';
param.Stimuli.Value{3,1}   = '';
param.Stimuli.Value{4,1}   = 'Keydown==32';
param.Stimuli.Value{5,1}   = '40s';
param.Stimuli.Value{1,2}   = '';
param.Stimuli.Value{2,2}   = '';
param.Stimuli.Value{3,2}   = video;
param.Stimuli.Value{4,2}   = '';
param.Stimuli.Value{5,2}   = '5s';

param.Stimuli.ColumnLabels = cell(n_configs+4,1);
param.Stimuli.ColumnLabels{1} = '1';
param.Stimuli.ColumnLabels{2} = '2';
for i=3:n_configs+2
    param.Stimuli.Value{1,i}   = ''; %caption
    param.Stimuli.Value{2,i}   = ''; %icon
    param.Stimuli.Value{3,i}   = ''; %av
    param.Stimuli.Value{4,i}   = ''; %interrupt
    jitter = round((2*(rand)-1),2);
    duration = 4 + jitter; % stimuli duration centered about 4s +/- 1s jitter. Note the first second is stimulation, and the rest is rest. 
    param.Stimuli.Value{5,i}   = sprintf('%ds',duration); %duration
    param.Stimuli.ColumnLabels{i,1} = sprintf('%d',i);
end

% pause between blocks stimuli
param.Stimuli.Value{1,end-1}   = 'Pause Between Blocks';
param.Stimuli.Value{2,end-1}   = '';
param.Stimuli.Value{3,end-1}   = '';
param.Stimuli.Value{4,end-1}   = '';
param.Stimuli.Value{5,end-1}   = '10s';
param.Stimuli.ColumnLabels{end-1,1} = sprintf('%d',i+1);

% end of run stimuli 
param.Stimuli.Value{1,end}   = 'End of Run';
param.Stimuli.Value{2,end}   = '';
param.Stimuli.Value{3,end}   = '';
param.Stimuli.Value{4,end}   = 'Keydown==32';
param.Stimuli.Value{5,end}   = '40s';
param.Stimuli.ColumnLabels{end,1} = sprintf('%d',i+2);
% for i=n_configs+3:2*n_configs+2  %if using empty blocks for isi
% param.Stimuli.Value(1,i)        = {''};
% param.Stimuli.Value(2,i)        = {''};
% param.Stimuli.Value(3,i)        = {''};
% param.Stimuli.Value(4,i)        = {'2'}; %isi
% end

param.Stimuli.RowLabels    = stimRowLabs;

res = checkStructFields(param.Stimuli);
%% set up sequence;
param.Sequence.Section       = 'Application';
param.Sequence.Type          = 'intlist';
param.Sequence.DefaultValue  = '';
param.Sequence.LowRange      = '';
param.Sequence.HighRange     = '';
param.Sequence.Comment       = '';
param.Sequence.Value         = '';
param.Sequence.Value{1,1}    = '1'; % startup screen
param.Sequence.Value{2,1}    = '2'; % begin movie
for i=1:size(trial_seq,1)
    param.Sequence.Value{i+2,1}  = sprintf('%d',trial_seq(i,7));
    % param.Sequence.NumericValue(i+2,1)  = trial_seq(i,7);
end
param.Sequence.Value{end+1,1} = sprintf('%d',max(trial_seq(:,7))+1); % end

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
for i=1:numConfigs
    param.StimulationTriggers.Value{1,2*i-1} = sprintf('%s == %d',TriggerExp,stimMat_cathode(i,7)); %Expression
    param.StimulationTriggers.Value{2,2*i-1} = sprintf('%d',stimMat_cathode(i,8));%Config
    param.StimulationTriggers.Value{3,2*i-1} = sprintf('%d',stimMat_cathode(i,5));%Electrode
    param.StimulationTriggers.Value{1,2*i}   = sprintf('%s == %d',TriggerExp,stimMat_anode(i,7)); %Expression
    param.StimulationTriggers.Value{2,2*i}   = sprintf('%d',stimMat_anode(i,8));%Config
    param.StimulationTriggers.Value{3,2*i}   = sprintf('%d',stimMat_anode(i,5));%Electrode

    param.StimulationTriggers.ColumnLabels{2*i-1,1} = sprintf('Trigger %d',2*i-1);
    param.StimulationTriggers.ColumnLabels{2*i,1}   = sprintf('Trigger %d',2*i);
end
param.StimulationTriggers.RowLabels = {'Expression';'Config ID'; 'Electrode(s)'};





%% Setup Dynamic Mode
param.DynamicConfiguration.Section = 'CereStim';
param.DynamicConfiguration.Type = 'int';
param.DynamicConfiguration.DefaultValue = '0';
param.DynamicConfiguration.LowRange = '0';
param.DynamicConfiguration.HighRange = '1';
param.DynamicConfiguration.Comment= '';
param.DynamicConfiguration.Value = {'1'};



%% Write Trial Sweep Parm
fname = sprintf('BLAES_stimsweep_run_%s.prm',parmName_dec);
filename = fullfile(parmDir,fname);
parameter_lines = convert_bciprm( param );
fid = fopen(filename, 'w');

for i=1:length(parameter_lines)
    fprintf( fid, '%s', parameter_lines{i} );
    fprintf( fid, '\r\n' );
end
fclose(fid);

%% Stimulation Testing Parm Generation

% set up Stim Configs
stim_configs = getUniqueConfigs(stimMat_cathode,stimMat_anode);
testing_param.StimulationConfigurations.Section = 'CereStim';
testing_param.StimulationConfigurations.Type = 'matrix';
testing_param.StimulationConfigurations.DefaultValue = '';
testing_param.StimulationConfigurations.LowRange = '';
testing_param.StimulationConfigurations.HighRange = '';
testing_param.StimulationConfigurations.Comment = 'Configurations for BLAES PARAM SWEEP';
testing_param.StimulationConfigurations.Value = cell(length(rowLabs),size(stim_configs,1));

for i=1:size(stim_configs,1)
    testing_param.StimulationConfigurations.Value{1, i}  = sprintf('%d',stim_configs(i,1));
    testing_param.StimulationConfigurations.Value{2, i}  = sprintf('%d',stim_configs(i,2));
    testing_param.StimulationConfigurations.Value{3, i}  = sprintf('%d',stim_configs(i,4));
    testing_param.StimulationConfigurations.Value{4, i}  = sprintf('%d',stim_configs(i,4));
    testing_param.StimulationConfigurations.Value{5, i}  = sprintf('%d',stim_configs(i,3));
    testing_param.StimulationConfigurations.Value{6, i}  = sprintf('%d',stim_configs(i,3));
    testing_param.StimulationConfigurations.Value{7, i}  = sprintf('%d',stim_configs(i,6));
    testing_param.StimulationConfigurations.Value{8, i}  = sprintf('%d',53);
    testing_param.StimulationConfigurations.Value{9, i}  = sprintf('%d',1);
    testing_param.StimulationConfigurations.Value{10,i}  = sprintf('%d',8);
    testing_param.StimulationConfigurations.ColumnLabels{i,1} = sprintf('Configuration %d',i);
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
testing_param.Sequence.Value         = {'1'}; % startup screen


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
load("C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\keymapping.mat") % this loads the keymapping matfile
keyboardMap = keyMap;
% new_map = cell(2,numConfigs);
% for i=1:(size(new_map,2)/2)
% new_map{1,2*i-1} = keyboardMap{1,i};
% new_map{2,2*i-1} = keyboardMap{2,i};
% new_map{1,2*i} = keyboardMap{1,i+size(new_map,2)/2};
% new_map{2,2*i} = keyboardMap{2,i+size(new_map,2)/2};
% end
% keyboardMap = new_map;
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
% set(gcf, 'Position', get(0, 'Screensize')); %fullscreen fig generation
set(gcf, 'Position', [1 1 1980 1080])
for i=1:numConfigs
    testing_param.StimulationTriggers.Value{1,2*i-1} = sprintf('%s == %s',TriggerExp,keyboardMap{2,i}); %Expression
    testing_param.StimulationTriggers.Value{2,2*i-1} = sprintf('%d',stimMat_cathode(i,8));%Config
    testing_param.StimulationTriggers.Value{3,2*i-1} = sprintf('%d',stimMat_cathode(i,5));%Electrode
    subplot(3,numConfigs/3,i)
    f_g = stimMat_cathode(i,6);
    cath_g = stimMat_cathode(i,5);
    anode_g = stimMat_anode(i,5);
    amp_g = stimMat_cathode(i,4);
    pw_g = stimMat_cathode(i,3);
    np_g = stimMat_cathode(i,2);
    [t,y] = generate_theta_burst_waveform(f_g,amp_g);
    plot(1000*t,y,'Linewidth',2)
    xlim([0,200])
    ylim([-0.5,2.5])
    titleblock = sprintf('Combination #%d, Key %s \n [Cathode, Anode]: [%d, %d]\namp: %d mA, f: %d Hz\npw: %d us, np: %d',stimLabel(i),keyboardMap{2,i},cath_g,anode_g,amp_g,f_g,pw_g,np_g);
    title(titleblock)
    ylabel('amp (mA)')
    xlabel('time (ms)')
    stimDescription{i,9} = keyboardMap{1,i};
    stimDescription{i,10} = keyboardMap{2,i};
    testing_param.StimulationTriggers.Value{1,2*i}   = sprintf('%s == %s',TriggerExp,keyboardMap{2,i}); %Expression
    testing_param.StimulationTriggers.Value{2,2*i}   = sprintf('%d',stimMat_anode(i,8));%Config
    testing_param.StimulationTriggers.Value{3,2*i}   = sprintf('%d',stimMat_anode(i,5));%Electrode

    testing_param.StimulationTriggers.ColumnLabels{2*i-1,1} = sprintf('Trigger %d',2*i-1);
    testing_param.StimulationTriggers.ColumnLabels{2*i,1}   = sprintf('Trigger %d',2*i);
end
testing_param.StimulationTriggers.RowLabels = {'Expression';'Config ID'; 'Electrode(s)'};
locDescriptions = cell2struct(stimDescription,stimDescription_labels,2);
saveas(1,cerestim_map_path);
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

for i=1:length(parameter_lines)
    fprintf( fid, '%s', parameter_lines{i} );
    fprintf( fid, '\r\n' );
end
fclose(fid);


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

function bci2ktools(BCI2KPath)
olddir = pwd;
cd(BCI2KPath);
cd tools, cd matlab;
bci2000path -AddToMatlabPath tools/matlab;
bci2000path -AddToMatlabPath tools/mex;
bci2000path -AddToSystemPath tools/cmdline;   % required so that BCI2000CHAIN can call the command-line tools
cd(olddir); % change directory back to where we were before
clear olddir;
end

function x = rm_combinations(input, idx2remove)

x = input;
if numel(idx2remove) > 0
x(idx2remove,:) = [];
end
end
