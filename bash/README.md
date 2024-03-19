# Bash Directory
This directory contains bash shell scripts to add builtin photodiode patches to videos. The description of each file is the output format.
- these scripts only work on Linux
	- adding windows functionality is likely not hard but I will not do it as I hate windows terminal
- Audio is preserved for these videos as well. 
- There is currently some loss of compression, and the encoding is software based since I am running this on WSL1 on windows. 
- Future updates will include hardware acceleration.

## dependencies
** FFmpeg **
- video processing software
	- https://ffmpeg.org/download.html
- A linux distro
- dos2unix (if on windows)
	- https://askubuntu.com/questions/1117623/how-to-install-dos2unix-on-a-ubuntu-app-on-a-windows-10-machine

## pathing
These files must be called from the root directory containing the videos, but the files themselves can be stored anywhere. 
- EX)
	- /Paradigms/tasks/param_sweep
		- video1.mp4
		- video2.mp4
	- /Path2code/bash
		- video_edit_x.sh
		- README.md (this file)
		
	to execute on a dir structure like this, open a linux terminal
		- cd /Paradigms/tasks/param_sweep
		- chmod +x Path2code/bash/video_edit_x.sh (only needed once to make this a linux executable)
		- dos2unix Path2code/bash/video_edit_x.sh (needed when there are windows line endings rather than unix)
		- ./Path2code/bash/video_edit_x.sh
			- this will execute the command and begin the batch processing
			- files will be stored in ../param_sweep/processed Directory in this example