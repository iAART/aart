from aart_func import *

print("\nThanks for using AART")
#print("Copyright (C) 2025, A. Cardenas-Avendano, H. Zhu & A. Lupsasca\n")

#BH's Spin
spin_case=0.94
#Observer's inclination  
i_case=17

# Distance to the BH in meters (default: M87)
dBH=5.214795112e23  
# BH mass-to-distance ratio (default: 1/psi= 6.2e9 Kg)
psi=1.07473555940836 
#Observer's distance in units of M
D_obs=1e5

#Velocity Profile for the gas

#Sub-Kepleniarity param
sub_kep=1.0;
#Radial velocity param
betar=1.0;
#Angular velocity param
betaphi=1.0;

# If equal to 1, the radon cuts profiles will be stored   
radonfile=0
# If equal to 1, the Beloborodov approximation will also be computed
bvapp=0
 
#For the Image resolution in the Bardeen coordinates
#Limits for the image [M]. It should coincide with the inoisy if used.
#If equal to 1, the sizes of the grids will be equal and an image can be computed
#by summing the contributions    
p_image=1
limits=25
#Resolution for the n=0 image [M]
dx0=0.02
#Resolution for the n=1 image [M]
dx1=0.02
#Resolution for the n=2 image [M]
dx2=0.02

# Projection angles for the radon transformation
radonangles=[0,90]

# Image treatment 
fudge=1.5 #Fudge factor (For n>0)

# Sample Equatorial Profile
i_fname="inoisy.h5"
#i_fname="hotspot.h5"

# Stationary assumes a single inoisy frame. "stationary" or "dynamical" 
disk="dynamical" 

# inoisy initial time frame for single images
i_frame=0

# Initial and final times in units of M
i_tM=0  
#Makes sense when is less than the inosy temporal length 
f_tM=12
#Number of snapshots in that range    
snapshots=12

isco = rms(spin_case)
# SU's parameters for the envelope 
# Just used for the analytical profiles
gammap=-3/2
mup=1-sqrt(1-spin_case**2)
sigmap=1/2 

#Magnetic field parametrs
cr=0.0
cphi=1.0

#Hotspot
#Radius of the hotspot
radhs=8#2.2
#Radius of the hotspot
velhs=0.05#0.2
rwidth = 0.5 #0.1

# Useful for disk visualizations or when studying truncated disks.
# 0: Neglected
# 1: Radii computed up to thar radius
# 2: Adds 5M to r_cutoff for interpolation purposes. 
imag_cut=0
# Cutoff radius   
r_cutoff=20

#The power of the redshift factor
gfactor=3

# Max baseline in G\lambda
maxbaseline=500 

# Number of points in the critical curve 
npointsS=100    

#For the parallel generation of images (movies)
nthreads=4

#Observer's inclination in radians
thetao=i_case*np.pi/180
#Disk's inclination  
#Current version is just implemented for equatorial models   
i_disk=90    
thetad=i_disk*np.pi/180

Gc=6.67e-11 # G constant [m^3 kg^-1 s^-2]
cc= 2.99792458e8 # c constant [m/s]
Msc=1.988435e30 # Solar Mass [Kg]

MMkg= 6.2e9*psi*Msc # [Kg]
MM=MMkg *Gc/cc**2 # Mass of the BH in meters, i.e., for M87(psi*6.2*10^9) psi ("Best fit") Solar Masses 

# Size of the real image in meters
sizeim_Real=(limits)*MM 
#1 microarcsec in radians
muas_to_rad = np.pi/648000 *1e-6 
fov_Real=np.arctan(sizeim_Real/(dBH))/muas_to_rad #muas
#print("FOV= ",np.round(2*fov,2),"muas")

#Path where the results will be stored
path = './Results/'

# Create a directory for the results
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
    print("A directory (Results) was created to store the results")

'''
MIT license
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.
'''