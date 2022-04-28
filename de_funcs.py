import numpy as np
import astropy.units as u
import pandas as pd
##########################functions to find ngp and cic####################################

def ngp(ngrid, slice_width, nbox, posi, offset):
    posi[2] -= offset
    new_pos = np.swapaxes(posi,0,1)
    boxfac = np.float32(ngrid/nbox)
    density = np.zeros((ngrid+1, ngrid+1, np.int32(boxfac*slice_width)+1)) #defining box of size ngrid

    new_pos = boxfac*new_pos
    ints = np.floor(new_pos).astype(np.int32)  
    
    for i in range(len(ints)):
        density[ints[i][0], ints[i][1], ints[i][2]] += 1
                
    return density


def cic(ngrid,slice_width, nbox, posi, offset):
    posi[2] -= offset
    new_pos = np.swapaxes(posi,0,1)
    boxfac = np.float32(ngrid/nbox)
    density = np.zeros((ngrid+1, ngrid+1, np.int32(boxfac*slice_width)+1), dtype=np.float64) #defining box of size ngrid
    new_pos = boxfac*new_pos
    
    frac = np.modf(new_pos)[0]
    ints = np.floor(new_pos).astype(np.int32)
                                                                                                                                                                                                                                                                                                                                                   
    for i in range(-1,len(ints)-1):
        for j in range(2):
            if j == 1 : frac[i][0] = 1 - frac[i][0]
            for k in range(2):
                if k == 1 : frac[i][1] = 1 - frac[i][1]
                for l in range(2):
                    if l == 1 : frac[i][2] = 1 - frac[i][2]
                    density[ints[i+j][0], ints[i+k][1], ints[i+l][2]] += (frac[i][0])*(frac[i][1])*(frac[i][2])
                    
            
#What it actually does
#         density[ints[i+1][0], ints[i][1], ints[i][2]] += (frac[i][0])*(1 - frac[i][1])*(1 - frac[i][2])
#         density[ints[i][0], ints[i+1][1], ints[i][2]] += (1 - frac[i][0])*(frac[i][1])*(1 - frac[i][2])
#         density[ints[i][0], ints[i][1], ints[i+1][2]] += (1 - frac[i][0])*(1 - frac[i][1])*(frac[i][2])
#         density[ints[i+1][0], ints[i+1][1], ints[i][2]] += (frac[i][0])*(frac[i][1])*(1 - frac[i][2])
#         density[ints[i][0], ints[i+1][1], ints[i+1][2]] += (1 - frac[i][0])*(frac[i][1])*(frac[i][2])
#         density[ints[i+1][0], ints[i][1], ints[i+1][2]] += (frac[i][0])*(1 - frac[i][1])*(frac[i][2])
#         density[ints[i][0], ints[i][1], ints[i][2]] += (1 - frac[i][0])*(1 - frac[i][1])*(1 - frac[i][2])
#         density[ints[i+1][0], ints[i+1][1], ints[i+1][2]] += (frac[i][0])*(frac[i][1])*(frac[i][2])
        
        

    return density

###########################################################################################################


