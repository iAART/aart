from aart_func import *
from params import *
#from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator

def round_up_to_even(f):
    return int(np.ceil(f / 2.) * 2)

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

fig, ax = plt.subplots()

CSg=ax.imshow(gs0,extent=[-lim0,lim0,-lim0,lim0],origin="lower");

CSr=ax.contour(rs0.reshape(N0,N0).T,levels=[radius],extent=[-lim0,lim0,-lim0,lim0],origin="lower");

#nPointCont=400

xs_HR=[]
ys_HR=[]
for i in range(len(CSr.allsegs[0])):
    v=CSr.collections[0].get_paths()[i].vertices
    for j in range(len(v[:,0])):
        xs_HR.append(v[:,0][j])
        ys_HR.append(v[:,1][j])
    
#xs_HR=np.array(xs_HR)[::int(CSr.allsegs[0][0].shape[0]/nPointCont)]
#ys_HR=np.array(ys_HR)[::int(CSr.allsegs[0][0].shape[0]/nPointCont)]
xs_HR=np.array(xs_HR)
ys_HR=np.array(ys_HR)

radiuscoords=np.array([xs_HR,ys_HR]).T

gsy,gsx=np.where(gs0==np.max(gs0))

x = np.linspace(-lim0, lim0, round_up_to_even(2*lim0/dx0))

ginterp = RegularGridInterpolator((x,x),gs0.T,method='linear')
rinterp = RegularGridInterpolator((x,x),rs0.reshape(N0,N0),method='linear')

redshifts =ginterp(radiuscoords)

maxredshift=np.max(redshifts)
minredshift=np.min(redshifts)

xfiner = np.linspace(x[gsx]-0.4, x[gsx]+0.4, round_up_to_even(2/(dx0/10)))
yfiner = np.linspace(x[gsy]-0.4, x[gsy]+0.4, round_up_to_even(2/(dx0/10)))

newpoints=np.column_stack((xfiner, yfiner))

radii=rinterp(newpoints)
gss=ginterp(newpoints)

maxgind=np.where(gss==np.max(gss))[0]

maxr=radii[maxgind][0]
maxg=gss[maxgind][0]

np.save("gmax.npy",maxredshift)
np.save("gmin.npy",minredshift)
np.save("rmax.npy",maxr)
np.save("grmax.npy",maxg)
