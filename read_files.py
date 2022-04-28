#============ class for Gadget2 snapshot header ====================
import numpy as np
import astropy.units as u
import pandas as pd
import os
import sys

class Gadget2_snapshot_hdr:

    def __init__(self, filename, DataFrame=False):

        if os.path.exists(filename):
            inp_file = filename
        else:
            print "file not found:"
            sys.exit()
        self.filename = filename   # keepig filename as self variable 
        self.Data_frame=DataFrame
        f = open(inp_file,'rb')    # rb is for read binary
        blocksize = np.fromfile(f,dtype=np.int32,count=1)
        if blocksize[0]!=256:
            raise ValueError("incorrect file format encountered when reading header of")

	#==== Gadget2 struct_io_header =====
	# number of particles of each type in this file 
        self.npart = np.fromfile(f,dtype=np.int32,count=6)

	# mass of particles of each type. If 0, then the masses are explicitly
	# stored in the mass-block of the snapshot file, otherwise they are omitted 
        self.massarr   = np.fromfile(f,dtype=np.float64,count=6)

	# time of snapshot file 
        self.time      = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# redshift of snapshot file 
        self.redshift  = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# flags whether the simulation was including star formation 
        self.flag_sfr  = (np.fromfile(f,dtype=np.int32,count=1))[0]

	# flags whether feedback was included (obsolete) 
        self.flag_feedback = (np.fromfile(f,dtype=np.int32,count=1))[0]

	# total number of particles of each type in this snapshot. This can be
	# different from npart if one is dealing with a multi-file snapshot. 
        self.npartTotal    = np.fromfile(f,dtype=np.uint32,count=6)

	# flags whether cooling was included  */
        self.flag_cooling  = (np.fromfile(f,dtype=np.int32,count=1))[0]

	# number of files in multi-file snapshot 
        self.num_files     = (np.fromfile(f,dtype=np.int32,count=1))[0]

	# box-size of simulation in case periodic boundaries were used 
        self.BoxSize       = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# matter density in units of critical density 
        self.Omega0        = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# cosmological constant parameter 
        self.OmegaLambda   = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# Hubble parameter in units of 100 km/sec/Mpc 
        self.HubbleParam   = (np.fromfile(f,dtype=np.float64,count=1))[0]

	# flags whether the file contains formation times of star particles 
        self.flag_stellarage = (np.fromfile(f,dtype=np.float32,count=1))[0]

	# flags whether the file contains metallicity values for gas and star particles 
        self.flag_metals =  (np.fromfile(f,dtype=np.float32,count=1))[0]

     	# High word of the total number of particles of each type 
        self.npartTotalHighWord = np.fromfile(f,dtype=np.uint32,count=6)

        # flags that IC-file contains entropy instead of u 
        self.flag_entropy_instead_u = (np.fromfile(f,dtype=np.float32,count=1))[0]

        # fills to 256 Bytes 
        self.unused = np.fromfile(f,dtype=np.byte, count=60)
        
        blk_check = np.fromfile(f,dtype=np.int32,count=1)
        if blocksize[0]!=blk_check[0]:
            raise ValueError("I/O:ERRORBLOCKS NOT MATCHING")
        del blocksize, blk_check 
        f.close()
        
   #---------------------------------------------------------------------------------

    def read_gadget_header(self):
        
	"""
	Simply prints out the header of the file
	""" 

        print 'npar=',self.npart
        print 'nall=',self.npartTotal
        print 'a=',self.time
        print 'z=',self.redshift
        print 'masses=',self.massarr*1e10,'Msun/h'
        print 'boxsize=',self.BoxSize,'kpc/h'
        print 'filenum=',self.num_files
        print 'cooling=',self.flag_cooling
        print 'Omega_m,Omega_l=',self.Omega0,self.OmegaLambda
        print 'h=',self.HubbleParam,'\n'
        print 'H0=', self.HubbleParam*100.* u.km*u.s**(-1)*u.Mpc**(-1)

        rhocrit=2.77536627e11 #h**2 M_sun/Mpc**3
        rhocrit=rhocrit/1e9 #h**2M_sun/kpc**3
        
        Omega_CDM=self.npartTotal[1]*self.massarr[1]*1e10/(self.BoxSize**3*rhocrit)
        print 'DM mass=%.5e  Omega_DM = %.5f'          %(self.massarr[1]*1e10, Omega_CDM)

    #------------- Read Particle Positions in comoving units kpc h^-1 ------------------
    
    def read_Pos(self, Double=False):
        
	"""
	Param
	Double: If True then positions are written in 
		double else they are written in 
	
	returns
	If Data_frame is flase then returns the numpy 
	array of shape (3, Npart) else Panda's data frame
	of the same shape.  
	"""

        with open(self.filename,'rb') as f:   # rb is for read binary
                offset = 4+256+4

                f.seek(offset, os.SEEK_CUR) # https://www.geeksforgeeks.org/python-seek-function/(???????)
                blocksize = np.fromfile(f,dtype=np.int32,count=1) 
                #https://numpy.org/doc/stable/reference/generated/numpy.fromfile.html
                
                if not Double:
                        dt = np.dtype((np.float32,3))  # Positions are in float not double
                else:     
                        dt = np.dtype((np.float64,3))  # Positions are in float not double

                Pos = np.fromfile(f,dtype=dt,count=self.npart[1]) # https://numpy.org/doc/stable/reference/arrays.dtypes.html
                blk_check = np.fromfile(f,dtype=np.int32,count=1)

                if blocksize[0]!=blk_check[0]:
                        raise ValueError("I/O:ERRORBLOCKS NOT MATCHING")
                
                del blk_check, blocksize
 
        
        if self.Data_frame == True:
            df = pd.DataFrame(Pos)
            return df
        else:
            return Pos

    #-------- Read Particle velocities in Gadget internal units as km s^-1 ---------

    def read_Vel(self):
            
    
	"""
	returns
	If Data_frame is flase then returns the numpy 
	array of shape (3, Npart) else Panda's data frame
	of the same shape.  
	"""
       
        f = open(self.filename,'rb')  # rb is for read binary
        
        # Calculate the offset from the beginning of the 
	# file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)

        offset = 4+256+8
        
        #Skip all the particle Position (four bytes per position * three coordinates* number of particles)
        offset += 4 * 3 * self.npart[1]
        f.seek(offset+4, os.SEEK_CUR)
        
        blocksize = np.fromfile(f,dtype=np.int32,count=1)
        dt = np.dtype((np.float32,3))  # velocities are in float not double
        Vel = np.fromfile(f,dtype=dt,count=self.npart[1])

        blk_check = np.fromfile(f,dtype=np.int32,count=1)
        
        if blocksize[0]!=blk_check[0]:
            raise ValueError("I/O:ERRORBLOCKS NOT MATCHING")
                    
        del blk_check, blocksize
        f.close()    
        
        if self.Data_frame == True:
            df = pd.DataFrame(Vel)
            return df
        else:
            return Vel
        

    #----------------- Read Particle ID --------------------------

    def read_ID(self):
    
        f = open(self.filename,'rb')    # rb is for read binary
        
        # Calculate the offset from the beginning of the 
	# file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
        
        offset = 4+256+8
        #Skip all the particle Position
        offset += 4 * 3 * self.npart[1]
        offset+=8
        #Skip all the particle velocities
        offset += 4 * 3 * self.npart[1]

        f.seek(offset+4, os.SEEK_CUR)
        
        blocksize = np.fromfile(f,dtype=np.int32,count=1)

        Part_ID = np.fromfile(f,dtype=np.uint32, count=self.npart[1])
        blk_check = np.fromfile(f,dtype=np.int32,count=1)
        
        if blocksize[0]!=blk_check[0]:
            raise ValueError("I/O:ERRORBLOCKS NOT MATCHING")
            
        del blk_check, blocksize
        f.close()    

        return Part_ID