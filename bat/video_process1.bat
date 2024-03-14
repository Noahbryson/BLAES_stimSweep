@echo off
setlocal enabledelayedexpansion

REM Set your source and destination directories
set "SOURCE_DIR=C:\Users\nbrys\Documents\Brunner\code\BLAES_stimSweep\stimuli"
set "DEST_DIR=%SOURCE_DIR%\editing"

REM Create the destination directory if it doesn't exist
if not exist "%DEST_DIR%" mkdir "%DEST_DIR%"

REM Process each .mp4 file in the source directory
for %%F in ("%SOURCE_DIR%\*.mp4") do (
    echo Processing: %%~nxF
    ffmpeg -i "%%F" -vf "drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=black@1:t=fill" -codec:a copy "%DEST_DIR%\%%~nxF"
)

echo Processing complete.
pause
