%% Parse stimsweep stimcodes
% This function takes a bci2000 parameter file from a stimulus presentation
% experiment utilizing the cerestim, and creates a map of stimulus codes
% tied to their trigger expression. 
% ----------
% Parameters:
%   param: bci2000 param file
% ----------
% Returns
%   stimMap: struct
%       - Trigger Type: str 
%           expression of the trigger state
%       - Code: int 
%           value of the trigger state
%       - Static ID: str 
%           configuration ID as shown in the manual testing
%           PNG image
%       - Electrode: int 
%           stimulating electrode contact
%       - Electrode Type: str
%           anodic or cathodic leading stimulation
%       - Configuration ID: int 
%           Cerestim Stimulation Configuration ID
%       - Config: Struct
%           - Config: cell 
%               contains stimulation parameters
%           - Labels: cell 
%               contains stimulation parameter labels
                    
function stimMap = parse_stimsweep_stimcodes(param)
stimuliNames = param.Stimuli.ColumnLabels;
stimConfigs = param.StimulationConfigurations.Value;
stimTriggers = param.StimulationTriggers.Value;
stimMap = struct();
for i=1:size(stimTriggers,2)
    % disp(i)
    % if i==34
    %     disp(1)
    % end
    s = stimTriggers{1,i};
    regSpace = regexp(s,'\s*');
    regNum = regexp(s,'\d*');
    if length(regSpace) > 1
        regSpace = regSpace(1);
    end
    stimMap(i).TriggerType = s(1:regSpace-1);
    stimMap(i).Code = str2double(s(regNum:end));
    stimMap(i).StaticID = stimuliNames(stimMap(i).Code);
    stimMap(i).Electrode = str2double(stimTriggers{3,i});
    if stimConfigs{1,str2double(stimTriggers{2,i})} == '1'
        eType = 'Cathode';
    else
        eType = 'Anode';
    end
    stimMap(i).ElectrodeType = eType; 
    stimMap(i).ConfigurationID = str2double(stimTriggers{2,i});
    temp = struct();
    temp.Config = stimConfigs(:,str2double(stimTriggers{2,i}));
    temp.Labels = param.StimulationConfigurations.RowLabels;
    % stimMap(i).Config = stimConfigs{:,str2double(cathodeTriggers{2,i})};
    stimMap(i).Config = temp;
end
end
