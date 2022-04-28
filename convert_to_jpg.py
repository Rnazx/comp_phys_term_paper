from PIL import Image
import numpy as np
for p in range(200,510,10):
    p = np.float64(p)
    img = Image.open(str(p)+'cic_10Mpc_1024g.eps')
    rgb_img = img.convert('RGB')
    rgb_img.save('/home/rnazx/eps_plots/jpg_plots/'+str(p)+'.jpg',quality = 95)