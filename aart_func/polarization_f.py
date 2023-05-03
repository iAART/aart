from aart_func import *
from params import * 

#Magnetic field components
#Radial profile as an example

def Bt_f(r,th,Delta,ut,uphi,a):
    return 0

def Br_f(r,th,Delta,ut,uphi,a):
    return sqrt(Delta)/r 

def Bth_f(r,th,Delta,ut,uphi,a):
    return 0 

def Bphi_f(r,th,Delta,ut,uphi,a):
    return 0

def kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,ur,uphi,Bt,Br,Bth,Bphi):
    """
    Calculates the Penrose–Walker constant (Eq. 6 P1)
    """    
    Omega=uphi/ut
    iota=-ur/ut

    k1=(a**2*(Bth*(-4 + r)*r*lamb + Bphi*sqrt(eta)*(2*(-1 + r)*r + sqR*iota*redshift_sign) + Bth*(sqR**2 - 2*r*(r**2 + r**3 + lamb**2))*Omega + sqrt(eta)*(-2*Bt*(-1 + r)*r + Br*sqR*redshift_sign)*Omega) + r**2*((Bphi*sqrt(eta) + Bth*lamb)*((-2 + r)*r + sqR*iota*redshift_sign) + (Bth*(-r**4 + sqR**2) + sqrt(eta)*(-(Bt*(-2 + r)*r) + Br*sqR*redshift_sign))*Omega) + a**4*(-(Bth*r*(2 + r)*Omega) + sqrt(eta)*(Bphi - Bt*Omega)) + a**3*(sqrt(eta)*lamb*(-Bphi + Bt*Omega) + Bth*r*(2 + (4 + r)*lamb*Omega)) + a*(-(sqrt(eta)*(sqR*(Br + Bt*iota)*redshift_sign + (-2 + r)*r*lamb*(Bphi - Bt*Omega))) + Bth*(-sqR**2 + r*(2*lamb**2 + r*(2*r - lamb**2 - sqR*iota*redshift_sign + r*(2 + r)*lamb*Omega)))))/(iota*(2*a*Bth*r + a**2*Bphi*sqrt(eta) + (-2 + r)*r*(Bphi*sqrt(eta) + Bth*lamb)) + (a**2 + (-2 + r)*r)*(Br*sqrt(eta) + Bth*sqR*redshift_sign)*Omega)
    k2=(a*r*(4*Br*lamb + 4*Bt*iota*lamb - r*(Br + (Bt + Bth*sqrt(eta))*iota)*lamb + Bphi*iota*(r**3 - r*eta + 2*(eta + lamb**2)) + Bphi*(-2 + r)*sqR*redshift_sign + Br*(r**3 - r*eta + 2*(eta + lamb**2))*Omega - (-2 + r)*sqR*(Bt + Bth*sqrt(eta))*redshift_sign*Omega) + r*(-2*Bt*eta*iota + Bt*r*eta*iota - Bphi*r**3*iota*lamb - 2*Bt*iota*lamb**2 + Bt*r*iota*lamb**2 + Br*(-2 + r)*(eta + lamb**2) + 2*Bphi*sqR*lamb*redshift_sign - Bphi*r*sqR*lamb*redshift_sign + Bth*sqrt(eta)*(r**3*iota + (-2 + r)*sqR*redshift_sign) - Br*r**3*lamb*Omega + Bt*(-2 + r)*sqR*lamb*redshift_sign*Omega) + a**3*(Bphi*(r*(2 + r) - eta)*iota + Bphi*sqR*redshift_sign + (Br*r*(2 + r) - Br*eta - sqR*(Bt + Bth*sqrt(eta))*redshift_sign)*Omega) + a**2*(r*iota*(Bth*r*sqrt(eta) - Bphi*(4 + r)*lamb) + sqR*(Bth*sqrt(eta) - Bphi*lamb)*redshift_sign + Bt*(-2*r*iota + eta*iota + sqR*lamb*redshift_sign*Omega) + Br*(eta - r*(2 + (4 + r)*lamb*Omega))))/(iota*(2*a*Bth*r + a**2*Bphi*sqrt(eta) + (-2 + r)*r*(Bphi*sqrt(eta) + Bth*lamb)) + (a**2 + (-2 + r)*r)*(Br*sqrt(eta) + Bth*sqR*redshift_sign)*Omega)

    return np.array((k1+1j*k2)/r,dtype=np.complex_)

def KDisk(r,thetad,a,lamb,eta,Delta,redshift_sign,sqR):
    """
    Cunningham's four velocity outside the ISCO (Eq. B17-B18 P1)
    """
    xi=sqrt(r**3-3*r**2+2*a*r**(3/2))
    ut=(r**(3/2)+a)/xi
    uphi=1/xi

    return kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,ut*0.0,uphi,Bt_f(r,thetad,Delta,ut,uphi,a),Br_f(r,thetad,Delta,ut,uphi,a),Bth_f(r,thetad,Delta,ut,uphi,a),Bphi_f(r,thetad,Delta,ut,uphi,a))

def KGas(r,thetad,a,lamb,eta,Delta,redshift_sign,sqR):
    """    
    Cunningham's four velocity inside the ISCO (Eq. B17-B18 P1)
    """
    
    #Eqns (P1 2)
    Delta=r**2 - 2*r + a**2

    #Eqns (P3 B13)
    lambe=((isco**2 - 2*a*sqrt(isco) + a**2))/(isco**(3/2) - 2*sqrt(isco) + a)
    #Eqns (P3 B12)
    H=(2*r - a*lambe)/Delta

    #Eqns (P3 B14)
    gamma=sqrt(1 - 2/3 *1/isco)

    #Eqns (P3 B9-B11)
    ut=gamma*(1 + 2/r *(1 + H))
    uphi=gamma/r**2*(lambe + a*H)
    ur=-np.sqrt(2/3/isco)*(isco/r - 1)**(3/2)

    return kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,ur,uphi,Bt_f(r,thetad,Delta,ut,uphi,a),Br_f(r,thetad,Delta,ut,uphi,a),Bth_f(r,thetad,Delta,ut,uphi,a),Bphi_f(r,thetad,Delta,ut,uphi,a))

def kappa(grid,mask,N,rs,redshift_sign,a,thetao):
    """
    Computes the linear polarization
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]

    #Conserved quantities
    lamb = -alpha*sin(thetao)
    eta = (alpha**2-a**2)*cos(thetao)**2+beta**2

    polk = np.zeros(rs.shape[0],dtype=np.complex_)

    redshift_sign = redshift_sign[mask]

    Delta=rs**2-2*rs+a**2

    sqR=sqrt((rs**2+a**2-a*lamb)**2-Delta*(eta+(lamb-a)**2))

    polk[rs>=isco]= KDisk(rs[rs>=isco],thetad,a,lamb[rs>=isco],eta[rs>=isco],Delta[rs>=isco],redshift_sign[rs>=isco],sqR[rs>=isco])
    polk[rs<isco] = KGas(rs[rs<isco],thetad,a,lamb[rs<isco],eta[rs<isco],Delta[rs<isco],redshift_sign[rs<isco],sqR[rs<isco])

    r_p = 1+sqrt(1-a**2)
    polk[rs<=r_p] = 0
    
    k1=np.real(polk)
    k2=np.imag(polk)

    #Electric vector polarization angle EVPA (Eq. 5 P1)
    nu=-(alpha+a*sin(thetao))

    EVPA_d=sqrt((k1**2+k2**2)*(beta**2+nu**2))
    EVPA_i=(beta*k2-nu*k1)
    EVPA_j=(beta*k1+nu*k2)

    EVPA_i[rs<=2] = np.nan
    EVPA_j[rs<=2] = np.nan

    mask_d=EVPA_d>0
    EVPA_i[mask_d]/=EVPA_d[mask_d]
    EVPA_j[mask_d]/=EVPA_d[mask_d]

    PK = np.zeros(mask.shape,dtype=np.complex_)
    EVPA_x = np.zeros(mask.shape)
    EVPA_y = np.zeros(mask.shape)

    EVPA_x[mask] = EVPA_i
    EVPA_y[mask] = EVPA_j
    PK[mask]     = polk

    filename=path+"Polarization_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('PK',     data=PK.reshape(N,N).T)
    h5f.create_dataset('EVPA_x', data=EVPA_x.reshape(N,N).T)
    h5f.create_dataset('EVPA_y', data=EVPA_y.reshape(N,N).T)

    h5f.close()

    print("File ",filename," created.")

def kappa_bv(grid,mask,N,rs,redshift_sign,a,thetao):
    """
    Computes the linear polarization using the Beloborodov approximation
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]

    #Conserved quantities
    lamb = -alpha*sin(thetao)
    eta = (alpha**2-a**2)*cos(thetao)**2+beta**2

    polk = np.zeros(rs.shape[0],dtype=np.complex_)

    redshift_sign = redshift_sign[mask]

    Delta=rs**2-2*rs+a**2

    sqR=sqrt((rs**2+a**2-a*lamb)**2-Delta*(eta+(lamb-a)**2))

    polk[rs>=isco]= KDisk(rs[rs>=isco],thetad,a,lamb[rs>=isco],eta[rs>=isco],Delta[rs>=isco],redshift_sign[rs>=isco],sqR[rs>=isco])
    polk[rs<isco] = KGas(rs[rs<isco],thetad,a,lamb[rs<isco],eta[rs<isco],Delta[rs<isco],redshift_sign[rs<isco],sqR[rs<isco])

    r_p = 1+sqrt(1-a**2)
    polk[rs<=r_p] = 0
    
    k1=np.real(polk)
    k2=np.imag(polk)

    #Electric vector polarization angle EVPA (Eq. 5 P1)
    nu=-(alpha+a*sin(thetao))

    EVPA_d=sqrt((k1**2+k2**2)*(beta**2+nu**2))
    EVPA_i=(beta*k2-nu*k1)
    EVPA_j=(beta*k1+nu*k2)

    EVPA_i[rs<=2] = np.nan
    EVPA_j[rs<=2] = np.nan

    mask_d=EVPA_d>0
    EVPA_i[mask_d]/=EVPA_d[mask_d]
    EVPA_j[mask_d]/=EVPA_d[mask_d]

    PK = np.zeros(mask.shape,dtype=np.complex_)
    EVPA_x = np.zeros(mask.shape)
    EVPA_y = np.zeros(mask.shape)

    EVPA_x[mask] = EVPA_i
    EVPA_y[mask] = EVPA_j
    PK[mask]     = polk

    filename=path+"Polarization_bv_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('PK',     data=PK.reshape(N,N).T)
    h5f.create_dataset('EVPA_x', data=EVPA_x.reshape(N,N).T)
    h5f.create_dataset('EVPA_y', data=EVPA_y.reshape(N,N).T)

    h5f.close()

    print("File ",filename," created.")