from aart_func import *
from params import * 

print("Ray-tracing")

fnbands=path+"LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
mask0=h5f['mask0'][:]

if bvapp!=1:

	supergrid1=h5f['grid1'][:]
	mask1=h5f['mask1'][:]
	supergrid2=h5f['grid2'][:]
	mask2=h5f['mask2'][:]
	h5f.close()

	rt.rt(supergrid0,mask0,supergrid1,mask1,supergrid2,mask2)
	print("A total of",supergrid0.shape[0]+supergrid1.shape[0]+supergrid2.shape[0],"photons were ray-traced")

else:

	h5f.close()

	rt.rs_bv(supergrid0,mask0)
	print("A total of",supergrid0.shape[0],"photons were approximately ray-traced.")
