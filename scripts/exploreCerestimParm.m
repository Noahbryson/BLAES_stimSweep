close all
clc
loadBCI2kTools();
[~, cereStates, cereParm] = load_bcidat('dat\ECOGS001R01.dat');
stimConfig = cereParm.StimulationConfigurations;
stimTrig = cereParm.StimulationTriggers;
vidParm = read_bciprm('parms\video.prm');
dynamicParm = read_bciprm('C:\Paradigms\parms\CereStim\dynamic_config_frag.prm');

root = 'C:\Paradigms';
stimuliDir = fullfile(root,'\tasks\BLAES\BLAES_param_sweep\stimuli');
checkDir(stimuliDir);
parmDir  = fullfile(root,'\parms\BLAES\_BLAES_param_sweep');
checkDir(parmDir);


videos = dir(strcat(stimuliDir,'\','*.mp4'));
video_num = 4;
parmName_dec = videos(video_num).name(1:20);
video = fullfile(videos(video_num).folder,videos(video_num).name);
v=VideoReader(video);
v1 = read(v,1000);
imshow(v1)