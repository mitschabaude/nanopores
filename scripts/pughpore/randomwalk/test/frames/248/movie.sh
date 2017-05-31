#!/bin/bash

#ffmpeg -framerate 270 -start_number 00000027 -i %08d.png -r 30 -an -s 1920x1080 bind1.mp4

#ffmpeg -loop 1 -f image2 -framerate 270 -start_number 00000027 -i %08d.png -r 30 -an -s 1920x1080 bind1.mp4
#ffmpeg -loop 1 -f image2 -i 00002354.png -r 30 -t 3 -an -s 1920x1080 bind2.mp4
#ffmpeg -framerate 270 -start_number 00002355 -i %08d.png -r 30 -an -s 1920x1080 bind3.mp4
ffmpeg -f concat -i videos.txt -c copy bind_final.mp4
