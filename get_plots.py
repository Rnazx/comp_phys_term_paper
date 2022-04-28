
# coding: utf-8

#Importing necessary libraries
import matplotlib.pyplot as plt
plt.switch_backend('agg') #necessary for plotting in the cluster since no display server


import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import copy
from de_funcs import ngp, cic
from read_files import Gadget2_snapshot_hdr
import multiprocessing



#########################################Reading the Header#################################################




snapID = np.int32(22)
nfiles = np.int32(32)
dir_name= ("/mnt/data1/sandeep/Kalinga_run_Analysis/512Mpc_1024/kalinga_snapdir_seed_1690811/")

#/mnt/data1/sandeep/New_Gadget2_run/200Mpc_256/snapdir_seed_1690811/snap_06
print('Reading',dir_name)
if nfiles<=1:
    inp_file = dir_name + "snapshot_%03d"%(snapID)
else:
    inp_file = dir_name + "snapshot_%03d.%d"%(snapID, 0)

snapshot = Gadget2_snapshot_hdr(inp_file, DataFrame=True)
snapshot_header =  snapshot.read_gadget_header()

########################################################################################################



#########################################Making the Slice#####################################################




scale_fac = 1./(1.+snapshot.redshift)
unit_kpc_Mpc = np.float32(0.001)
vel_conv_factor = np.float32(1./np.sqrt(scale_fac))

posx = np.zeros(1, dtype=np.float64)
posy = np.zeros(1, dtype=np.float64)
posz = np.zeros(1, dtype=np.float64)

velx = np.zeros(1, dtype=np.float64)
vely = np.zeros(1, dtype=np.float64)
velz = np.zeros(1, dtype=np.float64)

for i in range(nfiles):

    if nfiles<=1:
        inp_file = dir_name + "snapshot_%03d"%(snapID)
    else:
        inp_file = dir_name + "snapshot_%03d.%d"%(snapID, i)
    print(inp_file)

    snapshot = Gadget2_snapshot_hdr(inp_file, DataFrame=False)
    Pos  = snapshot.read_Pos() # is u.kpc*snapshot.HubbleParam**(-1)

    Pos[:,0] = np.float64(Pos[:,0]*unit_kpc_Mpc)
    Pos[:,1] = np.float64(Pos[:,1]*unit_kpc_Mpc)
    Pos[:,2] = np.float64(Pos[:,2]*unit_kpc_Mpc)

    posx = np.concatenate((posx, Pos[:,0]))
    posy = np.concatenate((posy, Pos[:,1]))
    posz = np.concatenate((posz, Pos[:,2]))


    Vel  = snapshot.read_Vel() # to convert the Gadget velocity to comoving velocity
    Vel[:,0] *= np.float64(vel_conv_factor)
    Vel[:,1] *= np.float64(vel_conv_factor)
    Vel[:,2] *= np.float64(vel_conv_factor)

    velx = np.concatenate((velx, Vel[:,0]))        
    vely = np.concatenate((vely, Vel[:,1]))        
    velz = np.concatenate((velz, Vel[:,2]))        
    del Pos, Vel 


#####################################################################################################



###########################Creating a 2D plot from starting point p #############################################
Slice_width = np.float64(10)

def get_plot(p):

        p = np.float64(p)
        print('Starting from z = '+str(p)+ ' until z = '+str(p+Slice_width))
        index_z = (posz>=p + 0.0)*(posz<=p+Slice_width)

        positions = np.array([(posx[index_z])[1:], (posy[index_z])[1:], (posz[index_z])[1:]])


        pos = copy.deepcopy(positions)  
        pos1 = copy.deepcopy(positions) 


        ngrid = np.int32(1024)
        nbox = np.float64(512)
        lengrid = np.float64(nbox/ngrid)
        slice_width = np.float64(10)

        numgrid = lengrid*slice_width
        npart = np.int32(len(pos[0]))

        #finding ngp and cic 3-D arrays

        den = ngp(ngrid,slice_width,nbox,pos,p)
        denc = cic(ngrid,slice_width,nbox,pos1,p)

        #density check for ngp
        if int(den.sum()) == npart: print('NGP is consistent with the mass density distribution')
        print('Sum of the densities is ',den.sum())
        print('Sum of the masses is ', len(pos[0]))

        #density check for cic
        if int(denc.sum()) == npart : print('CIC is consistent with the mass density distribution')
        print('Sum of the densities is ',denc.sum())
        print('Sum of the masses is ', len(pos[0]))



        #slicing the density
        grid_num = 0
        print(den.shape)
        den_2d = den[ :,:,grid_num]
        denc_2d = denc[ :,:,grid_num]
        for i in range(1,int(numgrid)+1):
                den_2d += den[ :,:,grid_num + i]
                denc_2d += denc[ :,:,grid_num + i]


        min_den = den_2d.min()
        max_den = den_2d.max()

        min_denc = denc_2d.min()
        max_denc = denc_2d.max()

        #print(min_den, max_den, min_denc, max_denc)

        fig = plt.figure(figsize = (15,15))
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(denc_2d, 
                cmap=cm.afmhot, norm=LogNorm(vmin =0.8, vmax= 1000))

        fig.savefig('/mnt/home1/rionx/plotting_ngp_cic/eps_plots/'+str(p)+'cic_10Mpc_'+str(ngrid)+'g.eps', format='eps')#saving the plot as eps since we get better quality

        #plt.show()
        #saving the slice as an .npy file
        #np.save('pos_512.npy', positions, fix_imports = True)
####################################################################################################################

##############################parallel processing to get the plots faster##############################################

if __name__ == '__main__':
    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes=64)
    inputs = range(0,20,1)
    pool.map(get_plot, inputs)
    pool.close()



