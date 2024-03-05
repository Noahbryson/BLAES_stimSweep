#! ../prog/BCI2000Shell
@cls & ..\prog\BCI2000Shell %0 %* #! && exit /b 0 || exit /b 1\n
#######################################################################################
## $Id: StimulusPresentation_SignalGenerator.bat 6172 2020-12-31 19:20:56Z mellinger $
## Description: BCI2000 startup Operator module script. For an Operator scripting
##   reference, see
##   http://doc.bci2000.org/index/User_Reference:Operator_Module_Scripting
##
## $BEGIN_BCI2000_LICENSE$
##
## This file is part of BCI2000, a platform for real-time bio-signal research.
## [ Copyright (C) 2000-2021: BCI2000 team and many external contributors ]
##
## BCI2000 is free software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## BCI2000 is distributed in the hope that it will be useful, but
##                         WITHOUT ANY WARRANTY
## - without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## $END_BCI2000_LICENSE$
#######################################################################################
Change directory $BCI2000LAUNCHDIR
Show window; Set title ${Extract file base $0}
Reset system
Startup system localhost
Start executable NihonKohdenSource --local --LogKeyboard=1 --LogWebcam=1 
Start executable DummySignalProcessing   --local
Start executable StimulusPresentation    --local --EnableCereStim=1
Wait for Connected
Load parameterfile "C:/Paradigms/parms/BLAES/_BLAES_param_sweep/BLAES_stimsweep_run_500 Days Survival An.prm"
Load parameterfile "C:/Paradigms/parms/BLAES/_BLAES_param_sweep/BLAES_stimsweep_appWindow.prm"
Load parameterfile "C:/Paradigms/parms/Source.Nihon/raw/current_source.prm"
Load parameterfile "C:/Paradigms/parms/Source/transmit_channel_one.prm"
Load parameterfile "C:/Paradigms/parms/AudioExtension/fragment_audio_recording_on.prm"

set parameter SubjectName ECOG
set parameter SubjectSession 001
set parameter DataDirectory ..\data\BJH\BLAES_Stimulation_sweep
set parameter VisualizeSourceTime 8s
set parameter DisplayStream 1
visualize  watch StimulusCode
visualize  watch CereStimUploadStatus
visualize  watch CereStimStimulation
visualize  watch CereStimStimulationSoftware