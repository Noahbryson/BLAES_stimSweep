loadBCI2kTools();
[~, cereStates, cereParm] = load_bcidat('C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\dat\ECOGS001R01.dat');
stimConfig = cereParm.StimulationConfigurations;
stimTrig = cereParm.StimulationTriggers;
vidParm = read_bciprm('C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\parms\video.prm');