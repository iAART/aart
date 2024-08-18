from aart_func import *
from params import * 
from multiprocessing import get_context

print("Movies")

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

print("Reading inoisy file: ",i_fname)

hf = h5py.File(i_fname, 'r')

data = np.array(hf['data/data_env'])
#inoisy has periodic boudaries, so we need to copy wrap the data with one frame
data=np.concatenate((data,data[0,:,:][np.newaxis,:,:]),axis=0)

nt =data.shape[0] #inoisy time resolution
ni = data.shape[1] #inoisy x resolution
nj = data.shape[2] #inoisy y resolution

try: 
	xtstart = np.array(hf['params/x0start'])[0]
	xtend = np.array(hf['params/x0end'])[0]

	x1start = np.array(hf['params/x1start'])[0]
	x2start = np.array(hf['params/x2start'])[0]

	x1end = np.array(hf['params/x1end'])[0]
	x2end = np.array(hf['params/x2end'])[0]

except:
	xtstart = np.array(hf['params/x0start'])
	xtend = np.array(hf['params/x0end'])

	x1start = np.array(hf['params/x1start'])
	x2start = np.array(hf['params/x2start'])

	x1end = np.array(hf['params/x1end'])
	x2end = np.array(hf['params/x2end'])


x1 = np.linspace(x1start, x1end, ni) 
x2 = np.linspace(x2start, x2end, nj)

times = np.linspace(xtstart, xtend, nt)

h5py.File.close(hf)

fact=-(D_obs+2*np.log(D_obs))

t0-=fact
t1-=fact
t2-=fact

fact2=xtend/2-np.nanmax(t0)

t0+=fact2
t1+=fact2
t2+=fact2

print("AART starts!")

i_dt = xtend/nt
timeconversion=i_dt*MMkg*Gc/cc**3/(3600*24) # [days]

interpolated3_R=RegularGridInterpolator((times,x1,x2),data,fill_value=0,bounds_error=False,method='linear')

I0s = []
I1s = []
I2s = []

def mp_worker(tsnap):
	i_bghts0 = obsint.slow_light(supergrid0,mask0,sign0,spin_case,isco,rs0,phi0,np.mod(t0+tsnap,xtend), interpolated3_R,thetao)
	i_bghts1 = obsint.slow_light(supergrid1,mask1,sign1,spin_case,isco,rs1,phi1,np.mod(t1+tsnap,xtend), interpolated3_R,thetao)
	i_bghts2 = obsint.slow_light(supergrid2,mask2,sign2,spin_case,isco,rs2,phi2,np.mod(t2+tsnap,xtend), interpolated3_R,thetao)

	i_I0 = (i_bghts0).reshape(N0,N0).T
	i_I1 = (i_bghts1).reshape(N1,N1).T
	i_I2 = (i_bghts2).reshape(N2,N2).T

	print("Calculating an image at time t=%s (M)"%np.round(tsnap,5))
	return(i_I0,i_I1,i_I2)

p = get_context("fork").Pool(nthreads) #using n threads
    
if __name__ == '__main__':
	I0s,I1s,I2s = zip(*p.map(mp_worker, np.linspace(i_tM+i_frame,f_tM,snapshots)))

p.close()

filename=path+"Images_a_%s_i_%s.h5"%(spin_case,i_case)

h5f = h5py.File(filename, 'w')
h5f.create_dataset('bghts0', data=np.array(I0s))
h5f.create_dataset('bghts1', data=np.array(I1s))
h5f.create_dataset('bghts2', data=np.array(I2s))
h5f.create_dataset('tc', data=np.array([timeconversion]))
h5f.create_dataset('limits', data=np.array([limits]))

h5f.close()

print("Images ",filename," created.\n")