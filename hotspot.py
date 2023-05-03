from aart_func import *
from params import * 
from multiprocessing import get_context

print("Intensity for a Flare")

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

fnbands=path+"Rays_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

rs0=h5f['rs0'][:]
sign0=h5f['sign0'][:]
t0=h5f['t0'][:]
phi0=h5f['phi0'][:]

rs1=h5f['rs1'][:]
sign1=h5f['sign1'][:]
t1=h5f['t1'][:]
phi1=h5f['phi1'][:]

rs2=h5f['rs2'][:]
sign2=h5f['sign2'][:]
t2=h5f['t2'][:]
phi2=h5f['phi2'][:]

h5f.close()

I0s = []
I1s = []
I2s = []

def mp_worker(tsnap):
	bghts0 = obsint.flare_model(supergrid0,mask0,sign0,spin_case,rs0,phi0,t0,thetao,rwidth,tsnap)
	bghts1 = obsint.flare_model(supergrid1,mask1,sign1,spin_case,rs1,phi1,t1,thetao,rwidth,tsnap)
	bghts2 = obsint.flare_model(supergrid2,mask2,sign2,spin_case,rs2,phi2,t2,thetao,rwidth,tsnap)

	i_I0 = (bghts0).reshape(N0,N0).T
	i_I1 = (bghts1).reshape(N1,N1).T
	i_I2 = (bghts2).reshape(N2,N2).T

	print("Calculating an image at time t=%s (M)"%np.round(tsnap,5))
	return(i_I0,i_I1,i_I2)

p = get_context("fork").Pool(nthreads) #using n threads
    
if __name__ == '__main__':
	I0s,I1s,I2s = zip(*p.map(mp_worker, np.linspace(i_tM+i_frame,f_tM,snapshots)))

p.close()

filename=path+"Intensities_a_%s_i_%s_rhs_%s_velhs_%s_rwidth_%s.h5"%(spin_case,i_case,radhs,velhs,rwidth)
h5f = h5py.File(filename, 'w')

h5f.create_dataset('bghts0', data=np.array(I0s))
h5f.create_dataset('bghts1', data=np.array(I1s))
h5f.create_dataset('bghts2', data=np.array(I2s))

h5f.close()

print("File ",filename," created.")
	


