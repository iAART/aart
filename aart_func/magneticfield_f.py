from aart_func import *
from params import * 

#I am not proud of this module. I wrote this in less than an hour. For sure there are
#better ways to write the equations.
#THIS JUST WORKS FOR KEPLERIAN CIRCULAR ORBITS!  

#### Cos(\theta_B)^2, where \theta_B is the between the field b^\mu and photon momentum k^\mu
def Br_f(r,Delta,a):
    return cr/r**2 

def Bth_f(r,Delta,a):
    return 0 

def Bphi_f(r,Delta,a):
    return -cphi/Delta

def anglemagout(r,Delta,a,signr,signth,eta,lamb,rms):
    return ((2*a*r**1.5 + (-3 + r)*r**2)*(-(Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(r) + r**2)) + (Bth_f(r,Delta,a)*signth*(2*a + (-3 + r)*np.sqrt(r))*r**1.5*np.sqrt(eta))/(a + r**1.5) + Bphi_f(r,Delta,a)*(a + (-2 + r)*np.sqrt(r))*lamb + (Br_f(r,Delta,a)*signr*(2*a + (-3 + r)*np.sqrt(r))*r**1.5*np.sqrt(-((a**2 + (-2 + r)*r)*(eta + (a - lamb)**2)) + (a**2 + r**2 - a*lamb)**2))/((a**2 + (-2 + r)*r)*(a + r**1.5)))**2)/(((Bth_f(r,Delta,a)**2*(2*a + (-3 + r)*np.sqrt(r))**2*r**5)/(a + r**1.5)**2 + (Br_f(r,Delta,a)**2*(2*a + (-3 + r)*np.sqrt(r))**2*r**5)/((a**2 + (-2 + r)*r)*(a + r**1.5)**2) + Bphi_f(r,Delta,a)**2*(a + (-2 + r)*np.sqrt(r))*(a**2 + (-2 + r)*r)*(a + r**1.5) - Bphi_f(r,Delta,a)**2*(a**2 + (-2 + r)*r)*(a**2 - 2*a*np.sqrt(r) + r**2))*(a + r**1.5 - lamb)**2)

def anglemagin(r,Delta,a,signr,signth,eta,lamb,rms):
    return (3*(a**2 + (-2 + r)*r)**2*((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r) - (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2) + (Bth_f(r,Delta,a)*signth*r**2*(a**2 + (-2 + r)*r)*np.sqrt(eta))/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms)))) + (r**2*(a**2 + (-2 + r)*r)*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*rms))*(a**2*r + (-2 + r)*rms**2 + 2*a*np.sqrt(rms)*(-r + rms))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + rms)*np.sqrt(rms))))*lamb)/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms)))) + (signr*r**2*(Br_f(r,Delta,a) + (np.sqrt(2/3)*(-1 + rms/r)**1.5*((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r) - (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/np.sqrt(rms))*np.sqrt(-((a**2 + (-2 + r)*r)*(eta + (a - lamb)**2)) + (a**2 + r**2 - a*lamb)**2))/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms)))))**2)/(((3*Bth_f(r,Delta,a)**2*r**4*(a**2 + (-2 + r)*r)**2*(a + (-2 + rms)*np.sqrt(rms))**2*rms)/((-2 + 3*rms)*(a**3*r + r**3*(-2 + rms)*np.sqrt(rms) + a**2*np.sqrt(rms)*(r*(-2 + rms) + 2*rms) + a*(r**3 - 2*rms**2))**2) + (3*r**6*(a**2 + (-2 + r)*r)*rms*(Br_f(r,Delta,a) + (np.sqrt(2/3)*(-1 + rms/r)**1.5*((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r) - (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/np.sqrt(rms))**2)/((-2 + 3*rms)*(a**4 + a**2*(-2 + r)*r - (a**2 + r**2)**2 + (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms)))**2) + (-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2))*((-1 + 2/r)*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)) - (2*a*r*(a**2 + (-2 + r)*r)*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*rms))*(a**2*r + (-2 + r)*rms**2 + 2*a*np.sqrt(rms)*(-r + rms))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + rms)*np.sqrt(rms)))))/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms))))) + (r**2*(a**2 + (-2 + r)*r)*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*rms))*(a**2*r + (-2 + r)*rms**2 + 2*a*np.sqrt(rms)*(-r + rms))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + rms)*np.sqrt(rms))))*((-2*a*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/r + (r**2*(a**2 + (-2 + r)*r)*(r**2 + (a**2*(2 + r))/r)*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*rms))*(a**2*r + (-2 + r)*rms**2 + 2*a*np.sqrt(rms)*(-r + rms))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*(-r + rms)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(rms) + rms**2))/np.sqrt(2*a*rms**1.5 + (-3 + rms)*rms**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + rms)*np.sqrt(rms)))))/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms))))))/(np.sqrt(1 - 2/(3.*rms))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms)))))*((np.sqrt(3 - 2/rms)*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(rms) + rms**2))/(a + (-2 + rms)*np.sqrt(rms))))/r**2 - (np.sqrt(3 - 2/rms)*(a**2*r + (-2 + r)*rms**2 + 2*a*np.sqrt(rms)*(-r + rms))*lamb)/(r*(a + (-2 + rms)*np.sqrt(rms))) + np.sqrt(2)*signr*(-1 + rms/r)**1.5*np.sqrt((-((a**2 + (-2 + r)*r)*(eta + (a - lamb)**2)) + (a**2 + r**2 - a*lamb)**2)/rms))**2)

def cos2angB_f(grid,mask,N,rs,redshift_sign,a,thetao):

    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    redshift_sign = redshift_sign[mask]

    Delta=rs**2-2*rs+a**2

    #Conserved quantities
    lamb = -alpha*sin(thetao)
    eta = (alpha**2-a**2)*cos(thetao)**2+beta**2

    cos2angB = np.zeros(rs.shape[0])

    cos2angB[rs>=isco]= anglemagout(rs[rs>=isco],Delta[rs>=isco],a,redshift_sign[rs>=isco],np.sign(beta[rs>=isco]),eta[rs>=isco],lamb[rs>=isco],isco)
    cos2angB[rs<isco]= anglemagin(rs[rs<isco],Delta[rs<isco],a,redshift_sign[rs<isco],np.sign(beta[rs<isco]),eta[rs<isco],lamb[rs<isco],isco)

    r_p = 1+sqrt(1-a**2)

    cos2angB[rs<=r_p] = np.nan

    cos2thB = np.zeros(mask.shape)
    cos2thB[mask] = cos2angB

    return cos2thB