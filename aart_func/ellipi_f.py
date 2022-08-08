from aart_func import *
from params import * 

def ellipi(n,phi,m): 
    csc=1/np.sin(phi)**2
    return 1/3*n*elliprj(csc-1,csc-m,csc,csc-n)+ellipf(phi,m)

def ellippi_lim(n,phi,m):
    #take arrays n, m, phi can either be a float or an array
    #restricted to work on -np.pi/2 < phi < np.pi, 0<n<1, m<1.
    mask = np.ones(n.shape)
    mask[n>1/(np.sin(phi)**2+1e-16)] = 0
    mask[m>1] = 0
    mask[phi>np.pi/2] = 0
    mask[phi<-np.pi/2] = 0
    result = np.zeros(n.shape)
    result=ellipi(mask*n.real,mask*phi.real,mask*m.real)
    mask[mask==0] = np.nan
    return(result*mask)

#extending the range of phi to (-pi/2,pi)
def ellippi1(n,phi,m):
    phi = np.ones(n.shape)*phi
    result = np.zeros(n.shape)
    mask = np.ones(n.shape,dtype=bool)
    mask[phi>np.pi/2] = False
    result[mask] = ellippi_lim(n[mask],phi[mask],m[mask])
    mask1 = np.invert(mask)
    result[mask1] = 2*ellippi_lim(n[mask1],np.pi/2,m[mask1]) - ellippi_lim(n[mask1],np.pi-phi[mask1],m[mask1])
    return(result)

#extending the range of n to larger than 1
def ellippi2(n,phi,m):
    return(-ellippi1(m/n, phi, m) + ellipf(phi, m)
           +1/(2* np.sqrt((n-m)*(n-1)/n))* np.log(np.abs((np.sqrt(1 - m*np.sin(phi)**2) + 
            np.sqrt((n-m)*(n-1)/n)*np.tan(phi))/(np.sqrt(1 - m*np.sin(phi)**2) - 
                                                          np.sqrt((n-m)*(n-1)/n)*np.tan(phi)))))

def ellippi(n,phi,m):
    phi = np.ones(n.shape)*phi
    result = np.zeros(n.shape)
    result[n<1] = ellippi1(n[n<1],phi[n<1],m[n<1])
    result[n>1] = ellippi2(n[n>1],phi[n>1],m[n>1])
    return(result)

