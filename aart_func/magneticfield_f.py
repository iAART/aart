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

def anglemagout(r,Delta,a,signr,signth,eta,lamb,isco):
    return ((2*a*r**2.5 + (-3 + r)*r**3)*(Bth_f(r,Delta,a)*signth*(2*a + (-3 + r)*np.sqrt(r))*r**4.5*(a**2 + (-2 + r)*r)*np.sqrt(eta) + Bphi_f(r,Delta,a)*(2 - 2*a - r)*r*(a + r**1.5)*(a**2 - 2*a*np.sqrt(r) + r**2)*(r**3 + a**2*(2 + r) - 2*a*lamb) - Bphi_f(r,Delta,a)*(a + (-2 + r)*np.sqrt(r))*r*(a + r**1.5)*(2*a - r**3 - a**2*(2 + r))*(2*a + (-2 + r)*lamb) + Br_f(r,Delta,a)*signr*(2*a + (-3 + r)*np.sqrt(r))*r**4.5*np.sqrt(-((a**2 + (-2 + r)*r)*(eta + (a - lamb)**2)) + (a**2 + r**2 + a*lamb)**2))**2)/(r**6*(a**2 + (-2 + r)*r)**2*(a + r**1.5)**2*((Bth_f(r,Delta,a)**2*(2*a + (-3 + r)*np.sqrt(r))**2*r**6)/(a + r**1.5)**2 + (Br_f(r,Delta,a)**2*(2*a + (-3 + r)*np.sqrt(r))**2*r**6)/((a**2 + (-2 + r)*r)*(a + r**1.5)**2) - Bphi_f(r,Delta,a)**2*(-2 + 2*a + r)*(a**2 - 2*a*np.sqrt(r) + r**2)**2 + Bphi_f(r,Delta,a)**2*(a + (-2 + r)*np.sqrt(r))**2*(-2*a + r**3 + a**2*(2 + r)))*(a + r**1.5 - lamb)**2)

def anglemagin(r,Delta,a,signr,signth,eta,lamb,isco):
    return ((Bth_f(r,Delta,a)*signth*r**2*(a**2 + (-2 + r)*r)*np.sqrt(eta))/(np.sqrt(1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))) + ((a**2 + (-2 + r)*r)*(-2*a + r**3 + a**2*(2 + r))*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*isco))*(a**2*r + (-2 + r)*isco**2 + 2*a*np.sqrt(isco)*(-r + isco))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*r**1.5*(-1 + isco/r)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(isco) + isco**2))/np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + isco)*np.sqrt(isco))))*(-a + lamb + (a*(a**2 + r**2 - a*lamb))/(a**2 - 2*r + r**2)))/(r*np.sqrt(1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))) - ((-2 + 2*a + r)*(3*a**4*Bphi_f(r,Delta,a) - 6*a**3*Bphi_f(r,Delta,a)*np.sqrt(isco) - 6*a*Bphi_f(r,Delta,a)*(-2 + r)*r*np.sqrt(isco) + 3*Bphi_f(r,Delta,a)*(-2 + r)*r*isco**2 + 3*a**2*Bphi_f(r,Delta,a)*(-2*r + r**2 + isco**2) + np.sqrt(6)*Br_f(r,Delta,a)*np.sqrt(r)*(r - isco)*np.sqrt(-1 + isco/r)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2))*(-(a*(a - lamb)) + ((a**2 + r**2)*(a**2 + r**2 - a*lamb))/(a**2 - 2*r + r**2)))/(3.*r**3*(a**2 + (-2 + r)*r)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)) + (signr*r**2*(Br_f(r,Delta,a) - (2*Br_f(r,Delta,a)*(r - isco)**3)/(3.*r**1.5*(a**2 + (-2 + r)*r)*np.sqrt(isco)) - (np.sqrt(2/3)*Bphi_f(r,Delta,a)*(-1 + isco/r)**1.5*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(np.sqrt(isco)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)))*np.sqrt((a**2 + r**2 + a*lamb)**2 - (a**2 - 2*r + r**2)*(eta + (-a + lamb)**2)))/(np.sqrt(1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))))**2/(((Bth_f(r,Delta,a)**2*r**6*(a**2 + (-2 + r)*r)**2)/((1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))**2) + (r**6*(a**2 + (-2 + r)*r)*(Br_f(r,Delta,a) - (2*Br_f(r,Delta,a)*(r - isco)**3)/(3.*r**1.5*(a**2 + (-2 + r)*r)*np.sqrt(isco)) - (np.sqrt(2/3)*Bphi_f(r,Delta,a)*(-1 + isco/r)**1.5*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(np.sqrt(isco)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)))**2)/((1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))**2) - ((-2 + 2*a + r)*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*r**1.5*(-1 + isco/r)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(isco) + isco**2))/np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2))*(3*a**4*Bphi_f(r,Delta,a) - 6*a**3*Bphi_f(r,Delta,a)*np.sqrt(isco) - 6*a*Bphi_f(r,Delta,a)*(-2 + r)*r*np.sqrt(isco) + 3*Bphi_f(r,Delta,a)*(-2 + r)*r*isco**2 + 3*a**2*Bphi_f(r,Delta,a)*(-2*r + r**2 + isco**2) + np.sqrt(6)*Br_f(r,Delta,a)*np.sqrt(r)*(r - isco)*np.sqrt(-1 + isco/r)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)))/(3.*r*(a**2 + (-2 + r)*r)*np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)) + (r**3*(a**2 + (-2 + r)*r)**2*(-2*a + r**3 + a**2*(2 + r))*(Bphi_f(r,Delta,a) + (np.sqrt(1 - 2/(3.*isco))*(a**2*r + (-2 + r)*isco**2 + 2*a*np.sqrt(isco)*(-r + isco))*(-((np.sqrt(2/3)*Br_f(r,Delta,a)*r**1.5*(-1 + isco/r)**1.5)/(a**2 + (-2 + r)*r)) + (Bphi_f(r,Delta,a)*(a**2 - 2*a*np.sqrt(isco) + isco**2))/np.sqrt(2*a*isco**1.5 + (-3 + isco)*isco**2)))/(r*(a**2 + (-2 + r)*r)*(a + (-2 + isco)*np.sqrt(isco))))**2)/((1 - 2/(3.*isco))*(-(a**2*(a**2 + (-2 + r)*r)) + (a**2 + r**2)**2 - (2*a*r*(a**2 - 2*a*np.sqrt(isco) + isco**2))/(a + (-2 + isco)*np.sqrt(isco)))**2))*(((a**2 - 2*a*np.sqrt(isco) + isco**2)*(-a + lamb + (a*(a**2 + r**2 - a*lamb))/(a**2 - 2*r + r**2)))/(r**2*np.sqrt(2*a*isco**1.5 - 3*isco**2 + isco**3)) - (np.sqrt(1 - 2/(3.*isco))*(-(a*(a - lamb)) + ((a**2 + r**2)*(a**2 + r**2 - a*lamb))/(a**2 - 2*r + r**2)))/r**2 - (np.sqrt(2/3)*signr*(-1 + isco/r)**1.5*np.sqrt((a**2 + r**2 + a*lamb)**2 - (a**2 - 2*r + r**2)*(eta + (-a + lamb)**2)))/(np.sqrt(r)*(a**2 - 2*r + r**2)))**2)

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
    
    filename=path+"MagneticAngle_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('cos2angB',data=cos2thB.reshape(N,N).T)
    h5f.close()

    print("File ",filename," created.")