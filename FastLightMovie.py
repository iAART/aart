from aart_func import *
from params import * 
from multiprocessing import get_context

print("Computing the fast-light Movie")

fnbands="./Results/LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

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

fnbands="./Results/Rays_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

rs0=h5f['rs0'][:]
sign0=h5f['sign0'][:]
phi0=h5f['phi0'][:]

rs1=h5f['rs1'][:]
sign1=h5f['sign1'][:]
phi1=h5f['phi1'][:]

rs2=h5f['rs2'][:]
sign2=h5f['sign2'][:]
phi2=h5f['phi2'][:]

h5f.close()

print("Reading inoisy file: ",i_fname)

hf = h5py.File(i_fname, 'r')

try:
    data = np.array(hf['data/data_env'])
except:
    data = np.array(hf['data/data_raw'])

# inoisy has periodic boundaries, so we need to wrap the data with one frame
data=np.concatenate((data,data[0,:,:][np.newaxis,:,:]),axis=0)

nt = data.shape[0]
ni = data.shape[1]
nj = data.shape[2]

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

print("Fast-Light calculation starts!")

i_dt = xtend/nt
timeconversion=i_dt*MMkg*Gc/cc**3/(3600*24) # [days]

maxintensity=np.nanmax(data)

interpolated3_R=RegularGridInterpolator((times,x1,x2),data,fill_value=0,bounds_error=False,method='linear')

# Define a COMMON time grid for both fast-light and slow-light
# This ensures frame-by-frame temporal alignment
t_frames = np.linspace(i_tM, f_tM, snapshots, endpoint=False)


# Worker: receives PHYSICAL OBSERVER TIME directly
def MovieWorker(tobs):

    # Fast-light: same emission time for all pixels
    i_bghts0 = obsint.fast_light(
        supergrid0, mask0, sign0, spin_case, isco,
        rs0, phi0, tobs, interpolated3_R, thetao
    )

    i_bghts1 = obsint.fast_light(
        supergrid1, mask1, sign1, spin_case, isco,
        rs1, phi1, tobs, interpolated3_R, thetao
    )

    i_bghts2 = obsint.fast_light(
        supergrid2, mask2, sign2, spin_case, isco,
        rs2, phi2, tobs, interpolated3_R, thetao
    )

    # Reshape images
    i_I0 = i_bghts0.reshape(N0, N0).T
    i_I1 = i_bghts1.reshape(N1, N1).T
    i_I2 = i_bghts2.reshape(N2, N2).T

    print(f"Calculating an image at time t={np.round(tobs, 5)} (M)")

    return i_I0, i_I1, i_I2


def main():
    p = get_context("spawn").Pool(nthreads)
    I0s, I1s, I2s = zip(*p.map(MovieWorker, t_frames))
    p.close()

    filename="./Results/FastLight_Images_a_%s_i_%s_%s.h5"%(spin_case,i_case,i_fname[:-3])

    h5f = h5py.File(filename, 'w')
    h5f.create_dataset('bghts0', data=np.array(I0s))
    h5f.create_dataset('bghts1', data=np.array(I1s))
    h5f.create_dataset('bghts2', data=np.array(I2s))

    print("Images ", filename, " created.\n")
    h5f.close()

if __name__ == '__main__':
    main()