from aart_func import *
from params import * 

print("Polarization")

fnbands=path+"LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
mask0=h5f['mask0'][:]
N0=int(h5f["N0"][0])

h5f.close()

if bvapp!=1:

	fnbands=path+"Rays_a_%s_i_%s.h5"%(spin_case,i_case)

	print("Reading file: ",fnbands)

	h5f = h5py.File(fnbands,'r')

	rs0=h5f['rs0'][:]
	sign0=h5f['sign0'][:]
	h5f.close()

	polarizationk.kappa(supergrid0,mask0,N0,rs0,sign0,spin_case,thetao)

else:
	
	fnrays=path+"Rays_bv_a_%s_i_%s.h5"%(spin_case,i_case)
	
	print("Reading file: ",fnrays)
	h5f = h5py.File(fnrays,'r')
	rs0_bv=h5f['rs0_bv'][:]
	sign0_bv=h5f['sign0_bv'][:]
	h5f.close()

	polarizationk.kappa_bv(supergrid0,mask0,N0,rs0_bv,sign0_bv,spin_case,thetao)