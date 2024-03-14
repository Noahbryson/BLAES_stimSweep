#!/bin/bash

# Directory containing the original videos
VIDEO_DIR="./stimuli"

# Subdirectory for processed videos
OUTPUT_DIR="${VIDEO_DIR}/editing"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .mp4 files in the video directory
for video in "${VIDEO_DIR}"/*.mp4; do
  # Extract filename without extension
  filename=$(basename "$video" .mp4)
  
  # Define output video path
  output="${OUTPUT_DIR}/${filename}.mp4"
  
  # Apply FFmpeg command to add a black box
  ffmpeg -i "$video" -vf "drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=black@1:t=fill" -codec:a copy "$output"
  
  echo "Processed: $output"
done

echo "Processing complete."
