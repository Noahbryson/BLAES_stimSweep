#!/bin/bash

VIDEO_DIR="./stimuli"
OUTPUT_DIR="${VIDEO_DIR}/processed"

mkdir -p "$OUTPUT_DIR"

for video in "${VIDEO_DIR}"/*.mp4; do
    filename=$(basename "$video" .mp4)
    output="${OUTPUT_DIR}/${filename}_processed.mp4"

    interval=15
    duration=1
    
    ffmpeg -i "$video" -filter_complex \
    "[0:v]scale=2560:-1,crop=2560:1440:0:(ih-1440)/2, \
     drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=black@1:t=fill, \
     drawbox=x=0:y=ih-ih*0.1:w=iw*0.1:h=ih*0.1:color=white@1:t=fill:enable='between(mod(t,${interval}),0,${duration})'" \
    -codec:a copy "$output"

    echo "Processed: $output"
done

echo "Processing complete."
