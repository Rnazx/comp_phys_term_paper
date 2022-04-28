#!pip install ffmpeg-python

import numpy as np

import ffmpeg
stream = ffmpeg.input('jpg_plots/*.jpg', pattern_type='glob', framerate=5)
stream = ffmpeg.output(stream, 'final.mp4')
ffmpeg.run(stream)