%% test parse_stimsweep_stimcodes.m
BCI2KPath = 'C:\BCI2000\BCI2000'; % Set the path of the BCI2000 main directory here
bci2ktools(BCI2KPath);
root = 'C:\Paradigms';
parmDir  = fullfile(root,'\parms\BLAES\_BLAES_param_sweep'); % path to write parameter file to
checkDir(parmDir);
parms = dir(strcat(parmDir,'\','*.prm'));
param = read_bciprm(fullfile(parms(2).folder,parms(2).name));
map = parse_stimsweep_stimcodes(param);





function StimMap = parse_stimsweep_stimcodes(param)
stimcodes = param.Stimuli.Value;
stimConfigs = param.StimulationConfigurations.Value;
stimTriggers = param.StimulationTriggers.Value;
config2code = cell();

StimMap = 0;
end