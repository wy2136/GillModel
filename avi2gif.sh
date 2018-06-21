#!/usr/bin/env bash
ffmpeg -i $1.avi -vf scale=640:-1 -r 10 -f image2pipe -vcodec ppm - | convert -delay 1 -loop 0 - $1.gif
