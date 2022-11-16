from aart_func import *
from params import * 

def imagetreat(image,radonangle,limsn,lims0):
    # Computes the radon cut of an image and scales it
    NN = image.shape[0]
    fov=fov_Real*(limsn/lims0)
    fov_rad=fov*1e-6*1./3600.*np.pi/180.
    dfovreal=fov_rad/NN
    radon_scaled = radon(image, theta=[radonangle]).flatten()    
    xaxis=np.linspace(-limsn,limsn,num=NN)
    return dfovreal*radon_scaled, xaxis

def radon_cut(radonangles,I0,I1,I2,supergrid0,supergrid1,supergrid2,Ncut=0):

    for i in range(len(radonangles)):

        radonangle = radonangles[i]
        radon0=imagetreat(I0,radonangle,supergrid0[-1,0],supergrid0[-1,0])
        radon1=imagetreat(fudge*I1,radonangle,supergrid1[-1,0],supergrid0[-1,0])
        radon2=imagetreat(fudge*I2,radonangle,supergrid2[-1,0],supergrid0[-1,0])

        R0 = interpolate.interp1d(radon0[1],radon0[0],fill_value=0, bounds_error=False,kind="linear")
        R1 = interpolate.interp1d(radon1[1],radon1[0],fill_value=0, bounds_error=False,kind="linear")
        R2 = interpolate.interp1d(radon2[1],radon2[0],fill_value=0, bounds_error=False,kind="linear")

        dx=np.min([dx0,dx1,dx2])
        xvalues =np.round(np.arange(-limits,limits+dx, dx),4)

        R=R0(xvalues)+R1(xvalues)+R2(xvalues)

        # Compute 1D FFT of the projection, shift and take modulus
        padding=16
        radonff = fft(R,padding*xvalues.shape[0]) #1D FFT of the projection
        radonshift=fftshift(radonff) # recenter FFT
        radonvisamp=np.abs(radonshift) # this is the visamp

        xaxis1=np.linspace(-fov_Real,fov_Real,num=xvalues.shape[0]) # in muas
        deltax1=xaxis1[1]-xaxis1[0]
        xfourier1=fftshift(fftfreq(padding*xvalues.shape[0],d=deltax1)) #re centered frequencies in muas^-1
        xfourier1/= 1e-6 * 1./3600. * np.pi/180. # in rad^-1
        xfourier1 /= 1e9 # in Glambda

        indice1=np.where((xfourier1>=0.) & (xfourier1<maxbaseline))[0] # select only the positive freqs the FFT is symmetrical anyway
        visamp=radonvisamp[indice1[0]:(indice1[-1]+1)] # select FFT for these positive freqs

        #Normalize the visamp
        norm=visamp[0]
        visamp/=norm

        #print("V(0)= ",norm)
        
        filename=path+"Visamp_%s_a_%s_i_%s_%s.h5"%(radonangle,spin_case,i_case,Ncut)
        
        h5f = h5py.File(filename, 'w')

        h5f.create_dataset('visamp', data=visamp)

        if radonfile==1:
            h5f.create_dataset('radon', data=R)

        if Ncut==0:
            h5f.create_dataset('freqs', data=xfourier1[indice1])
            if radonfile==1:
                h5f.create_dataset('x_radon', data=xvalues)

        h5f.close()

        print("File ",filename," created.")