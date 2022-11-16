from aart_func import *

def cbrt(x):
    '''
    Cubic root
    :param x: number to compute cbrt
    '''
    if x.imag==0:
        return np.cbrt(x)
    else:
        return x**(1/3)

def rms(a):
    '''
    ISCO value
    (Eq. B16 P1)
    :param a: BH spin
    '''
    Z1=1 + (1 - a**2)**(1/3) *((1 + a)**(1/3) + (1 - a)**(1/3))
    Z2=sqrt(3*a**2 + Z1**2)
    return (3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2*Z2)))
