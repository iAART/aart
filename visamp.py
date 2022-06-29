from aart_func import *
from params import * 

print("Visamp")

fnbands=path+"LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
N0=int(h5f["N0"][0])
supergrid1=h5f['grid1'][:]
N1=int(h5f["N1"][0])
supergrid2=h5f['grid2'][:]
N2=int(h5f["N2"][0])
h5f.close()

fint=path+"Intensity_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fint)

h5f = h5py.File(fint,'r')

bghts0=h5f['bghts0'][:]
bghts1=h5f['bghts1'][:]
bghts2=h5f['bghts2'][:]

h5f.close()

vamp.radon_cut(radonangles,bghts0,bghts1,bghts2,supergrid0,supergrid1,supergrid2)