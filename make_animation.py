#!pip install ffmpeg

import cv2
import numpy as np

images = ['/home/rnazx/comp_phys_term_paper/jpg_plots/'+str(float(p))+'.jpg' for p in range(0,20,1)]
print(images)

frameSize = (1080, 1080)

out = cv2.VideoWriter('movie.mp4',cv2.VideoWriter_fourcc(*'DIVX'), 5, frameSize)

for filename in images:
    img = cv2.imread(filename)
    out.write(img)

out.release()
