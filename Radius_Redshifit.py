from aart_func import *
from params import *
from scipy.interpolate import griddata

fnbands="./Results/LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

#The grid points for each lensing band
supergrid0=h5f['grid0'][:]
N0=int(h5f["N0"][0])
lim0=int(h5f["lim0"][0])

h5f.close()

fnrays="./Results/Rays_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnrays)

h5f = h5py.File(fnrays,'r')

rs0=h5f['rs0'][:]

h5f.close()	

fnrays="./Results/gfactors_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnrays)

h5f = h5py.File(fnrays,'r')

gs0=h5f['gs0'][:]

h5f.close()

radius=isco

fig, ax = plt.subplots(figsize=[5,5],dpi=400)

CSg=ax.imshow(gs0,extent=[-lim0,lim0,-lim0,lim0],origin="lower");

CSr=ax.contour(rs0.reshape(N0,N0).T,levels=[radius],extent=[-lim0,lim0,-lim0,lim0],origin="lower");

nPointCont=400

xs_HR=[]
ys_HR=[]
for i in range(len(CSr.allsegs[0])):
    v=CSr.collections[0].get_paths()[i].vertices
    for j in range(len(v[:,0])):
        xs_HR.append(v[:,0][j])
        ys_HR.append(v[:,1][j])
    
xs_HR=np.array(xs_HR)[::int(CSr.allsegs[0][0].shape[0]/nPointCont)]
ys_HR=np.array(ys_HR)[::int(CSr.allsegs[0][0].shape[0]/nPointCont)]
radiuscoords=np.array([xs_HR,ys_HR]).T

redshifts = griddata(supergrid0, gs0.reshape(N0*N0), radiuscoords, method='linear')

maxredshift=np.max(redshifts)

np.save("result.npy",maxredshift)

