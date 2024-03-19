#!/bin/bash

# Directory containing the original videos
VIDEO_DIR="./stimuli"

# Subdirectory for processed videos
OUTPUT_DIR="${VIDEO_DIR}/processed"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .mp4 files in the video directory
for video in "${VIDEO_DIR}"/*.mp4; do
    # Extract filename without extension
    filename=$(basename "$video" .mp4)
    
    # Define output video path
    mid="${OUTPUT_DIR}/${filename}_trimmed.mp4"
	output="${OUTPUT_DIR}/${filename}_processed.mp4"
    
    # Get the frame rate of the video
    fps=$(ffprobe -v 0 -of csv=p=0 -select_streams v:0 -show_entries stream=r_frame_rate "$video" | bc)
    
    # Calculate the duration in seconds when the white box should appear
    # Assuming you want the white box to appear for 1 second every 15 seconds
    duration=$(echo "scale=2; 1" | bc)
    interval=$(echo "scale=2; 15" | bc)
    
    # FFmpeg command to add a black box that's always there
    # and a white box that appears for the specified duration every interval
    ffmpeg -i "$video" -filter_complex "drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=black@1:t=fill, drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=white@1:t=fill:enable='between(mod(t,${interval}),0,${duration})'" -codec:a copy "$mid"
    ffmpeg -i "$video" -vf "scale=2560:-1, crop=2560:1440:0:(ih-1440)/2" -codec:a copy "$output"

    echo "Processed: $output"
done

echo "Processing complete."