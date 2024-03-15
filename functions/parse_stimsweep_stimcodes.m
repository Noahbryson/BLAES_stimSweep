%% Parse stimsweep stimcodes
function stimMap = parse_stimsweep_stimcodes(param)
stimConfigs = param.StimulationConfigurations.Value;
stimTriggers = param.StimulationTriggers.Value;
cathodeTriggers = stimTriggers(:,1:2:end);
anodeTriggers = stimTriggers(:,2:2:end);
stimMap = struct();
for i=1:size(cathodeTriggers,2)
    s = cathodeTriggers{1,i};
    reg = regexp(s,'\d*');
    s = s(reg:end);
    stimMap(i).Code = str2double(s);
    stimMap(i).Cathode = str2double(cathodeTriggers{3,i});
    stimMap(i).Anode = str2double(anodeTriggers{3,i});
    stimMap(i).ConfigurationID = str2double(cathodeTriggers{2,i});
    temp = struct();
    temp.Config = stimConfigs(:,str2double(cathodeTriggers{2,i}));
    temp.Labels = param.StimulationConfigurations.RowLabels;
    % stimMap(i).Config = stimConfigs{:,str2double(cathodeTriggers{2,i})};
    stimMap(i).Config = temp;
    
    
end
end
