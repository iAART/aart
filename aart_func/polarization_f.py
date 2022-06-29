from aart_func import *
from params import * 

#Radial profile as an example

def Bt_f(r,th,Delta,ut,uphi,a):
    return 0

def Br_f(r,th,Delta,ut,uphi,a):
    return sqrt(Delta)/r 

def Bth_f(r,th,Delta,ut,uphi,a):
    return 0 

def Bphi_f(r,th,Delta,ut,uphi,a):
    return 0

def kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,uphi,Bt,Br,Bth,Bphi):
    """
    """

    X=(r**2*(r**2+a**2)+2*a*r*(a-lamb)-ut/uphi*(r*(r-2)*lamb+2*a*r))/Delta

    return np.array(1/r*((r**2+a**2-ut/uphi*a)*(redshift_sign*sqR)/Delta -1j*(a-ut/uphi)*sqrt(eta)+((Bphi*ut/uphi-Bt)*((r**2+a**2-a*lamb)*sqrt(eta)-1j*(lamb-a)*redshift_sign*sqR)-X*((r**2+a**2-a*lamb)*Bth+1j*(lamb-a)*Br))/(redshift_sign*sqR*Bth + sqrt(eta)*Br)),dtype=np.complex_)

def KDisk(r,thetad,a,lamb,eta,Delta,redshift_sign,sqR):

    xi=sqrt(r**3-3*r**2+2*a*r**(3/2))
    ut=(r**(3/2)+a)/xi
    uphi=1/xi

    return kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,uphi,Bt_f(r,thetad,Delta,ut,uphi,a),Br_f(r,thetad,Delta,ut,uphi,a),Bth_f(r,thetad,Delta,ut,uphi,a),Bphi_f(r,thetad,Delta,ut,uphi,a))

def KGas(r,thetad,a,lamb,eta,Delta,redshift_sign,sqR):
    
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
    return kappaeq(r,a,lamb,eta,Delta,redshift_sign,sqR,ut,uphi,Bt_f(r,thetad,Delta,ut,uphi,a),Br_f(r,thetad,Delta,ut,uphi,a),Bth_f(r,thetad,Delta,ut,uphi,a),Bphi_f(r,thetad,Delta,ut,uphi,a))

def kappa(grid,mask,N,rs,redshift_sign,a,thetao):
    """
    Computes the Penrose-Walker constant
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]

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
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]

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