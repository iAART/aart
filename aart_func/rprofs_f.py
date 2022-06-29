from aart_func import *
from params import * 

# Johnson’s SU distribution (e.g., Eq.12 2008.03879)
def profile(r,a,gammap,mup,sigmap):
    return (1-np.tanh(speed_p*(r-cutoff_p)))/2*np.exp(-(1/2)*(gammap+ np.arcsinh((r-mup)/sigmap))**2)/sqrt((r-mup)**2+sigmap**2)