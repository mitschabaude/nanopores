#!/bin/bash

ffmpeg -framerate 290 -start_number 00000007 -i %08d.png -r 30 -an -s 1920x1080 movie.mp4
