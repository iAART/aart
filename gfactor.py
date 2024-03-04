from aart_func import *
from params import * 

print("Computing the redshift factor at each point in the image plane \n")

fnbands=path+"LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
mask0=h5f['mask0'][:]
N0=int(h5f["N0"][0])

supergrid1=h5f['grid1'][:]
mask1=h5f['mask1'][:]
N1=int(h5f["N1"][0])
	
supergrid2=h5f['grid2'][:]
mask2=h5f['mask2'][:]
N2=int(h5f["N2"][0])

h5f.close()

fnbands=path+"Rays_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

rs0=h5f['rs0'][:]
sign0=h5f['sign0'][:]

rs1=h5f['rs1'][:]
sign1=h5f['sign1'][:]


rs2=h5f['rs2'][:]
sign2=h5f['sign2'][:]

h5f.close()


i_g0 = obsint.gfactorf(supergrid0,mask0,sign0,spin_case,isco,rs0,thetao)
i_g1 = obsint.gfactorf(supergrid1,mask1,sign1,spin_case,isco,rs1,thetao)
i_g2 = obsint.gfactorf(supergrid2,mask2,sign2,spin_case,isco,rs2,thetao)

i_g0 = (i_g0).reshape(N0,N0).T
i_g1 = (i_g1).reshape(N1,N1).T
i_g2 = (i_g2).reshape(N2,N2).T

filename=path+"gfactors_a_%s_i_%s.h5"%(spin_case,i_case)

h5f = h5py.File(filename, 'w')
h5f.create_dataset('gs0', data=i_g0)
h5f.create_dataset('gs1', data=i_g1)
h5f.create_dataset('gs2', data=i_g2)

h5f.close()
