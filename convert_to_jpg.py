from PIL import Image
import numpy as np
for p in range(0,20,1):
    p = np.float64(p)
    img = Image.open('/home/rnazx/comp_phys_term_paper/eps_plots/'+str(p)+'cic_10Mpc_1024g.eps')
    rgb_img = img.convert('RGB')
    rgb_img.save('jpg_plots/'+str(p)+'.jpg',quality = 95)