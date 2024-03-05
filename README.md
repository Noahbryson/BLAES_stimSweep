# Repository for Parameter Generator and Parm Fragments for the BLAES stimSweep Task.
Version 1.0
------------
# BLAES_stimsweep_generator_CCEP.m
## Notes
### waveforms
  all stimulation currently occurs in one-second bursts, modulated by the carrier frequency. This will be updated in a later version
### timing
  1 Block takes (6s x number of configurations) when CCEPs are included
  1 Block takes (5s x number of configurations) when CCEPs are excluded
  NOTE: The start of a new BLOCK will have a 10s delay.

### exclusions
  excluding a 1 mA condition will cause this to break unfortunately,
  since the pseduorandomization uses recursive methods to ensure rules
  about amplitude and channels are not violated. 
## Parameters
- BCI2kPath: String
  - Path on your machine to the BCI2000 distribution. Required for .mex files needed to write parameter files
- cathodeChannels: array
  - cathode leading stim channels passed to stimulator, be sure to align the index of this with anodeChannels for correct pairs
- anodeChannels: array
  - anode leading stim channels passed to stimulator, be sure to align the index of this with cathodeChannels for correct pairs 
- stimAmps: array
  -  array of stimulation amplitudes in mA 
- pulseWidth: array
  - array of pulsewidths in us 
- frequencies: array
  - array of frequencies (Hz)
- num_pulses: array, optional
  - array of the number of pulses for frequency. Indicies must be matched with frequencies variable. Will be generated if empty
- carrier_freq: float, default = 8 Hz
  - carrier frequency of the stimulation waveforms
- num_blocks: int
  - number of blocks in a run. one trial per configuration in a block
- video_num: int
  - video number in the stimuli directory to write to the file (alphabetical order)
- conditions2remove: array[int]
  - conditons removed after testing, please enter the number as it appears on the testing image
  - currently, conditions must be removed in pairs at minimum.
- generateTest: bool
  - boolean to generate the testing .prm file and associated figure. do not close the figure before the script finishes
- rmHighestChargeCondition: bool
  - boolean to remove conditions with the highest stimultion charge (max(stimAmp) U max(pulseWidth))
- allowCCEPs: bool
  - boolean to allow leading and lagging CCEP stimulation around stimulation.
  - CCEPS are 1s before and after stimulation
- root: Path
  - base directory to where stimuli are currently stored and where parms should be stored. typically 'C:\Paradigms'
  - inside the root structure, make sure the following directories exist:
    - ..\tasks\BLAES\BLAES_param_sweep\stimuli
    - ..\parms\BLAES\_BLAES_param_sweep



 ## Dependencies
 1. MATLAB (developed on version 2023a)
 2. BCI2000 revision 7897 or later
 3. Video files in an appropriate directory to be loaded into MATLAB
 4. keys2numbers.mat file (can be found in this repostiory, will clone automatically, do not move)
 

# Directory Structure
- parms
  - BLAES_stimsweep_appWindowBLAES_stimsweep_appWindow.prm
    - .prm file configuring stimulus presentation application
    - this can be edited or swapped out to meet your needs, but is the default in the .bat files     
- functions
  - sliceVideo.m (Deprecated)
    - used to cicrumvent a BCI2000 bug
   
  - modulateSquareWave.m (Deprecated)
    - function to generate modulated waveforms for testing image
   
  -  generate_theta_burst_waveform.m
    -  function to generate modulated waveforms for testing image
 
  - frequency_feasibility.m (Deprecated)
    - just my scratch script
 
  - checkDir.m
    - verifies if a string is a directory.
 
# TODO
1. make CCEP configurations dynamic based on the number of channels.
  Currently CCEPs are statically set for two channels, this should not be
  hard to implement I am just tired of working on this.

2. add an option for only including CCEPs
on every other block. 

3. make an option to include user feedback so I can use the same script to
  generate Tao's stimulation titration experiment

4. make train frequency and duration configurable
