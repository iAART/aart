from aart_func import *
from params import * 

#Angular momentum and carter constant
def conserved_quantities(alpha,beta,theta_o,a):
    """
    Computes the conserved quantities for a photon landing at (alpha,beta) on the image plane
    :param alpha: x coordinate on the image plane (perpendicular to the projected spin axis)
    :param beta: y coordinate (parallel to the projected spin axis)
    :param tehta_o: inclination angle of the observer from the black hole spin axis
    :param a: spin of the black hole (in unit M=1)
    :returns: angular momentum and carter constant
    """
    # Eqs. (58-59 P2)
    lam = -alpha*np.sin(theta_o)
    eta = (alpha**2-a**2)*np.cos(theta_o)**2+beta**2
    return(lam,eta)

def cuberoot(a):
    '''
    Computes cube root for an array of real number or the principle root for complex number
    :param a: array of numbers
    :returns: cuberoot of a
    '''
    a_real_ind = np.ones(a.shape)
    a_real_ind[np.abs(a.imag)>1e-12] = 0
    croa1 = np.cbrt((a*a_real_ind).real)
    croa2 = (a-a*a_real_ind)**(1/3)
    return(croa1+croa2)

#calculate radial turning points
def radial_turning_points(alpha,beta,lam,eta,a):
    '''
    Computes radial turning points for a photon orbit
    :param alpha: x coordinate on the image plane (perpendicular to the projected spin axis)
    :param beta: y coordinate (parallel to the projected spin axis)
    :param lam: angular momentum
    :param eta: carter constant
    :param a: spin of black hole
    :returns: the four roots for the radial turning points r1,r2,r3,r4. 
    '''
    # Eqs. (A1-A5 P2)
    A = a**2-eta-lam**2
    B = 2*(eta+(lam-a)**2)
    C = -a**2*eta
    
    P = -A**2/12-C
    Q = -A/3*((A/6)**2-C)-B**2/8
    
    # Eqs. (A7 P2)
    pp = -Q/2+np.sqrt((P/3)**3+(Q/2)**2+0*1j)
    mm = -Q/2-np.sqrt((P/3)**3+(Q/2)**2+0*1j)
    w_p = cuberoot(pp)
    w_m = cuberoot(mm)

    # Eqs. (A6 P2)
    z = np.sqrt((w_p+w_m)/2 - A/6)
    
    # Eqs. (A8 P2)
    r1 = - z - np.sqrt(-A/2-z**2+B/4/z)
    r2 = - z + np.sqrt(-A/2-z**2+B/4/z)
    r3 = + z - np.sqrt(-A/2-z**2-B/4/z)
    r4 = + z + np.sqrt(-A/2-z**2-B/4/z)
    return(r1,r2,r3,r4)

#we can now calculate the angular turning points
def angular_turning_points(alpha,beta,lam,eta,a):
    '''
    Computes angular turning points for a photon orbit
    :param alpha: x coordinate on the image plane (perpendicular to the projected spin axis)
    :param beta: y coordinate (parallel to the projected spin axis)
    :param lam: angular momentum
    :param eta: carter constant
    :param a: spin of black hole
    :returns: u_p, u_m for evaluating the angular integrals;
              and theta_p, theta_m the angular turning points 
    '''
    # Eqs. (11 P2)
    Delta_theta = (1-(eta+lam**2)/a**2)/2
    u_p = Delta_theta+np.sqrt(Delta_theta**2+eta/a**2)
    u_m = Delta_theta-np.sqrt(Delta_theta**2+eta/a**2)
    
    # Eqs. (10 P2)
    theta_p = np.arccos(-np.sqrt(u_p))
    theta_m = np.arccos(np.sqrt(u_p))
    
    return(u_p,u_m,theta_p,theta_m)


def angular_integrals(mbar,beta,u_p,u_m,theta_p,theta_m,pm_o,theta_o,a):
    '''
    Computes angular path integrals along the photon trajectories (P2 section II.A)
    :param mbar: the number of angular turning points encountered along the trajectory
    :param beta: y coordinate (parallel to the projected spin axis)
    :param u_p,u_m: to construct elliptical parameter for the integrals
    :param theta_p,theta_m: the angular turning points 
    :param pm_o: the sign of p_theta (theta momentum) at the observer 
    :param tehta_o: inclination angle of the observer from the black hole spin axis
    :param a: spin of black hole
    :returns: G_theta, G_phi, G_t angular path integrals
    '''
    
    k = u_p/u_m
    
    # Eqs. (16 P2)
    K = ellipk(k)
    E_prime = (ellipe(k)-K)/2/k
    PI = ellippi(u_p,np.pi/2,k)
    
    arg = (np.arcsin(np.cos(theta_o)/np.sqrt(u_p)))
    
    # Eqs. (12-14 P2)
    F_o = ellipf(arg,k)
    E_prime_o = (ellipeinc(arg,k)-F_o)/2/k
    PI_o = ellippi(u_p,arg,k)
    #source terms are zero as we are assuming emission only from the equitorial plane 
    
    #m = m+1 for the disk on the back side of the black hole
    H_beta = np.zeros(beta.shape)
    H_beta[beta>=0] = 1
    m = mbar+H_beta

    # Eqs. (20-22 P2), where the source terms are zero as the emission is from the equitorial plane
    G_theta = (2*m*K-pm_o*F_o)/a/np.sqrt(-u_m)
    G_phi = (2*m*PI-pm_o*PI_o)/a/np.sqrt(-u_m)
    G_t = -2*u_p*(2*m*E_prime-pm_o*E_prime_o)/a/np.sqrt(-u_m)
    
    return(G_theta,G_phi,G_t)

def source_radius2_inf(r1,r2,r3,r4,G_theta):
    '''
    Computes radius for the equitorial source of a photon with Type II trajectory
    (outside the critical curve, one turning point, scattering) in Boyer-Lindquist coordinates
    :param r1-4: radial turning points
    :param G_theta: the angular path integral, G_theta=I_r=Mino time

    :returns:  radius of the equitorial source
    '''
    # Eqs. (31-32 P1)
    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r21 = (r2-r1)

    k2 = r32*r41/r31/r42

    # Eqs. (30 P1)
    F2 = ellipf(np.arcsin(np.sqrt(r31/r41).real),(k2).real)
    sn_square = np.square(ellipj(1/2*np.sqrt(r31*r42).real*G_theta-F2, (k2).real)[0])
    rs2 = np.nan_to_num((r4*r31-r3*r41*sn_square)/(r31-r41*sn_square))
    return(rs2)

def source_radius2(r,r1,r2,r3,r4,G_theta):
    '''
    Computes radius for the equitorial source of a photon with Type II trajectory
    (outside the critical curve, one turning point, scattering) in Boyer-Lindquist coordinates
    :param r: Observer's radius
    :param r1-4: radial turning points
    :param G_theta: the angular path integral, G_theta=I_r=Mino time

    :returns:  radius of the equitorial source
    '''
    # Eqs. (31-32 P1)
    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r21 = (r2-r1)

    k2 = r32*r41/r31/r42

    x2 = np.sqrt((r-r4)*r31/(r-r3)/r41).real
    arg = np.arcsin(x2)

    F2 = ellipf(arg.real,k2.real)
    sn_square = np.square(ellipj(1/2*np.sqrt(r31*r42).real*G_theta-F2, (k2).real)[0])
    rs2 = np.nan_to_num((r4*r31-r3*r41*sn_square)/(r31-r41*sn_square))
    return(rs2)

def source_radius3_inf(r1,r2,r3,r4,G_theta):
    '''
    Computes radius for the equitorial source of a photon with Type III trajectory
    (inside the critical curve, generated at the horizon, no turning points) in Boyer-Lindquist coordinates
    assuming that the observer is at infinity
    :param r1-4: radial turning points
    :param G_theta: the angular path integral, G_theta=I_r=Mino time
    :returns:  radius of the equitorial source
    '''

    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r21 = (r2-r1)

    #Eqs (B57 P3)
    A = np.sqrt(r32*r42).real
    B = np.sqrt(r31*r41).real
    #Eq (B59 P3)
    k3 = ((A+B)**2 - r21**2)/(4*A*B)

    F3 = ellipf(np.arccos((A-B)/(A+B)),(k3).real)

    # Eqs. (B74-75 P2)
    cn = ellipj(np.sqrt(A*B)*G_theta-F3,(k3).real)[1]
    #Eq (B75 P3)
    rs3 = np.nan_to_num(((B*r2-A*r1)+(B*r2+A*r1)*cn)/((B-A)+(B+A)*cn))
    return(rs3)

def source_radius3(r,r1,r2,r3,r4,G_theta):
    '''
    Computes radius for the equitorial source of a photon with Type III trajectory
    (inside the critical curve, generated at the horizon, no turning points) in Boyer-Lindquist coordinates
    :param r: Observer's radius
    :param r1-4: radial turning points
    :param G_theta: the angular path integral, G_theta=I_r=Mino time

    :returns:  radius of the equitorial source
    '''
    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r21 = (r2-r1)

    #Eqs (B57 P2)
    A = np.sqrt(r32*r42).real
    B = np.sqrt(r31*r41).real
    #Eq (B59 P2)
    k3 = ((A+B)**2 - r21**2)/(4*A*B)

    F3 = ellipf(np.arccos((A*(r-r1)-B*(r-r2))/(A*(r-r1)+B*(r-r2))).real,(k3).real)

    # Eqs. (B74-75 P1)
    cn = ellipj(np.sqrt(A*B)*G_theta-F3,(k3).real)[1]

    rs3 = np.nan_to_num(((A*r1-B*r2)-(A*r1+B*r2)*cn)/((A-B)-(A+B)*cn))
    return(rs3)

def radial_potential(r,a,lam,eta):
    """
    Evaluates the radial effective potential, roots of which are the turing points. 
    :params r: radius
    :params a: spin of the black hole (in units M=1)
    :params lam: angular momentum
    :params eta: Carter constant
    :return: value of the radial potential
    """
    #Eqs (2,5 P2)
    Delta = r**2-2*r+a**2
    return((r**2+a**2-a*lam)**2-Delta*(eta+(lam-a)**2))

#(indefinite) auxiliary functions for antiderivatives of the three radial integral in radial case 2.
def radial_case2_antiderivative(r,r1,r2,r3,r4,a,lam,eta):
    '''
    Computes auxiliary integrals for the antiderivatives of the radial path integrals in Type II radial trajectory (P3 Appendix B.2)
    :param r: radius of the equitorial photon source
    :param r1-r4: radial turning points
    :param a: spin of black hole
    :param lam: angular momentum
    :param eta: carter constant 
    :returns I0,I1,I2,Ip,Im: values of the auxiliary integrals, combinations of which yeilds the antiderivatives
    '''
    
    #Eqs (3 P2)
    r_m = 1-np.sqrt(1-a**2)
    r_p = 1+np.sqrt(1-a**2)
    
    #differences in radial turning points, and inner outer horizon
    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r43 = (r4-r3)
    r21 = (r2-r1)
    rp3 = (r_p-r3)
    rm3 = (r_m-r3)
    rp4 = (r_p-r4)
    rm4 = (r_m-r4)
    
    #Eqs (B13 P3)
    k = r32*r41/r31/r42
    
    #Eqs (B35 P3)
    x2 = np.sqrt((r-r4)*r31/(r-r3)/r41)
    arg = np.arcsin(x2)

    #Eqs (B40-43 P3)
    F2 = 2/np.sqrt(r31*r42)*ellipf(arg.real,k.real)
    E2 = np.sqrt(r31*r42)*ellipeinc(arg.real,k.real)
    PI21 = 2/np.sqrt(r31*r42)*ellippi((r41/r31).real,arg.real,k.real)
    PI2p = 2/np.sqrt(r31*r42)*r43/rp3/rp4*ellippi((rp3*r41/rp4/r31).real,arg.real,k.real)
    PI2m = 2/np.sqrt(r31*r42)*r43/rm3/rm4*ellippi((rm3*r41/rm4/r31).real,arg.real,k.real)
    
    
    #Eqs (B36-9 P3)
    I0 = F2
    I1 = r3*F2 + r43*PI21
    I2 = sqrt(radial_potential(r,a,lam,eta))/(r-r3)-(r1*r4+r2*r3)/2*F2-E2
    Ip = - PI2p - F2/rp3
    Im = - PI2m - F2/rm3
    return(I0,I1,I2,Ip,Im)

#radial integrals in case 2.
def radial_case2(rs,ro,r1,r2,r3,r4,a,beta,lam,eta,redshift_sign):
    """
    Evaluates definite radial integrals for case 2 photon trajectory
    :param rs: radius of the equitorial photon source
    :param ro: radius of the observer
    :param r1-r4: radial turning points
    :param a: spin of black hole
    :param beta: y axis of the observer plane
    :param lam: angular momentum
    :param eta: carter constant 
    :redshift_sign: sign of the redshift associated with a photon
    :returns I0,I1,I2,Ip,Im: values of the definite auxiliary integrals.  
    """
    #Eqs (B25 P3)
    #where w = v_r evaluated at the source
    w = redshift_sign
    w[w==1] = 0
    w[w==-1] = 1

    I0s,I1s,I2s,Ips,Ims = radial_case2_antiderivative(rs,r1,r2,r3,r4,a,lam,eta)
    
    #turning points terms
    I0t,I1t,I2t,Ipt,Imt = radial_case2_antiderivative(r4,r1,r2,r3,r4,a,lam,eta)
    
    #observer terms
    I0o,I1o,I2o,Ipo,Imo = radial_case2_antiderivative(ro,r1,r2,r3,r4,a,lam,eta)

    #Eqs (B25 P3)
    #analagous to that in function "angular_integrals"
    I0 = I0o-I0s + 2*w*(I0s-I0t)
    I1 = I1o-I1s + 2*w*(I1s-I1t)
    I2 = I2o-I2s + 2*w*(I2s-I2t)
    Ip = Ipo-Ips + 2*w*(Ips-Ipt)
    Im = Imo-Ims + 2*w*(Ims-Imt)
    
    return(I0,I1,I2,Ip,Im)

#(indefinite) antiderivatives of the three radial integral in radial case 3.
def radial_case3_antiderivative(r,r1,r2,r3,r4,a):
    '''
    Computes auxiliary integrals for the antiderivatives of the radial path integrals in Type III radial trajectory (P3 Appendix B.3)
    :param r: radius of the equitorial photon source
    :param r1-r4: radial turning points
    :param a: spin of black hole
    :returns I0,I1,I2,Ip,Im: values of the auxiliary integrals, combinations of which yeilds the antiderivatives
    '''
    #Eqs (3 P2)
    r_m = 1-np.sqrt(1-a**2)
    r_p = 1+np.sqrt(1-a**2)
    
    #differences in radial turning points, and inner outer horizon
    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    r43 = (r4-r3)
    r21 = (r2-r1)
    rp3 = (r_p-r3)
    rm3 = (r_m-r3)
    rp4 = (r_p-r4)
    rm4 = (r_m-r4)
    rp2 = (r_p-r2)
    rm2 = (r_m-r2)
    rp1 = (r_p-r1)
    rm1 = (r_m-r1)
    
    #Eqns (B57 P3)
    A = np.sqrt(r32*r42)
    B = np.sqrt(r31*r41)

    #Eqns (B58-59 P3)
    k = ((A+B)**2-r21**2)/(4*A*B)
    alpha0 = (B+A)/(B-A)
    x3 = (A*(r-r1)-B*(r-r2))/(A*(r-r1)+B*(r-r2))
    arg = np.arccos(x3)

    #Eqns (B66 P3)
    alphap = (B*rp2+A*rp1)/(B*rp2-A*rp1)
    alpham = (B*rm2+A*rm1)/(B*rm2-A*rm1)

    #Eqns (B61-62, B65 P3)
    def R1(alpha):
        p1 = np.sqrt((alpha**2-1)/(k+(1-k)*alpha**2))
        f1 = p1/2*np.log(np.abs(
            (p1*np.sqrt(1-k*np.sin(arg)**2)+np.sin(arg))/
            (p1*np.sqrt(1-k*np.sin(arg)**2)-np.sin(arg))))#natural log
        R1 = 1/(1-alpha**2)*(ellippi((alpha**2/(alpha**2-1)).real,arg.real,k.real)-alpha*f1)
        return(R1)
    
    F = ellipf(arg.real,k.real)
    E = ellipeinc(arg.real,k.real)
    
    #Eqns (B63-64 P3)
    def R2(alpha):
        R2 = 1/(alpha**2-1)*(F-alpha**2/(k+(1-k)*alpha**2)*(E-alpha*np.sin(arg)* \
                np.sqrt(1-k*np.sin(arg)**2)/(1+alpha*np.cos(arg)))) + \
                1/(k+(1-k)*alpha**2)*(2*k-alpha**2/(alpha**2-1))*R1(alpha)
        return(R2)

    #Eqns (B71 P3)
    F3 = 1/np.sqrt(A*B)*F

    #Eqns (B72 P3)
    PI31 = ((2*r21*np.sqrt(A*B)/(B**2-A**2))*R1(alpha0)).real
    PI32 = ((2*r21*np.sqrt(A*B)/(B**2-A**2))**2*R2(alpha0)).real

    #Eqns (B67-70 P3)
    I0 = F3
    I1 = ((B*r2+A*r1)/(B+A))*F3 + PI31
    I2 = ((B*r2+A*r1)/(B+A))**2*F3 + 2*((B*r2+A*r1)/(B+A))*PI31 + np.sqrt(A*B)*PI32
    Ip = -1/(B*rp2+A*rp1)*((B+A)*F3+2*r21*np.sqrt(A*B)/(B*rp2-A*rp1)*R1(alphap))
    Im = -1/(B*rm2+A*rm1)*((B+A)*F3+2*r21*np.sqrt(A*B)/(B*rm2-A*rm1)*R1(alpham))
    return(I0,I1,I2,Ip,Im,PI31)

#radial integrals in case 3.
def radial_case3(rs,ro,r1,r2,r3,r4,a):
    """
    Evaluates definite radial integrals for case III photon trajectory
    :param rs: radius of the equitorial photon source
    :param ro: radius of the observer
    :param r1-r4: radial turning points
    :param a: spin of black hole
    :returns I0,I1,I2,Ip,Im: values of the definite auxiliary integrals.  
    """
    #source terms
    I0s,I1s,I2s,Ips,Ims, PI31s = radial_case3_antiderivative(rs,r1,r2,r3,r4,a)

    #observer terms
    I0o,I1o,I2o,Ipo,Imo, PI31o = radial_case3_antiderivative(ro,r1,r2,r3,r4,a)

    #no turning points for case 3, hence the definite integral is just the difference
    #between the observer and the source.
    I0 = I0o-I0s
    I1 = I1o-I1s
    I2 = I2o-I2s
    Ip = Ipo-Ips
    Im = Imo-Ims
    PI31 = PI31o-PI31s

    return(I0,I1,I2,Ip,Im,PI31)

#radial integrals in both cases for every point on the observer plane.
def radial_integrals(rs,ro,r1,r2,r3,r4,a,beta, mask2, mask3, lam,eta,redshift_sign):
    """
    Evaluates definite radial integrals for case II and III photon trajectory, which are seperated by the critical curve
    :param rs: radius of the equitorial photon source
    :param ro: radius of the observer
    :param r1-r4: radial turning points
    :param a: spin of black hole
    :param beta: y axis of the observer plane
    :param mask2: 1 inside the region of case II photon (outside the critical curve), 0 otherwise
    :param mask3: 1 inside the region of case III photon (inside the critical curve, outside the horizon), 0 otherwise
    :param lam: angular momentum
    :param eta: carter constant 
    :redshift_sign: sign of the redshift associated with a photon
    :returns I_r,I_phi,I_t: values of the definite radial integrals.  
    """

    #Evaluates integrals for both regions

    I03,I13,I23,Ip3,Im3, PI31 = radial_case3(rs[mask3],ro,r1[mask3],r2[mask3],r3[mask3],r4[mask3],a)
    
    I02,I12,I22,Ip2,Im2 = radial_case2(rs[mask2],ro,r1[mask2],r2[mask2],r3[mask2],r4[mask2],a,beta[mask2],lam[mask2],eta[mask2],redshift_sign[mask2])
    

    I0 = np.zeros(mask2.shape)
    I1 = np.zeros(mask2.shape)
    I2 = np.zeros(mask2.shape)
    Ip = np.zeros(mask2.shape)
    Im = np.zeros(mask2.shape)
    PI1 = np.zeros(mask2.shape)

    #Stich up the solutions
    I0[mask2] = I02.real
    I0[mask3] = I03.real
    I1[mask2] = I12.real
    I1[mask3] = I13.real
    I2[mask2] = I22.real
    I2[mask3] = I23.real
    Ip[mask2] = Ip2.real
    Ip[mask3] = Ip3.real
    Im[mask2] = Im2.real
    Im[mask3] = Im3.real
    PI1[mask3] = PI31
    r_m = 1-np.sqrt(1-a**2)
    r_p = 1+np.sqrt(1-a**2)
    
    #Eqns (B1-3 P3)
    I_r = I0
    I_phi = 2*a/(r_p-r_m)*((r_p-a*lam/2)*Ip-(r_m-a*lam/2)*Im)
    I_t = (2)**2/(r_p-r_m)*(r_p*(r_p-a*lam/2)*Ip-r_m*(r_m-a*lam/2)*Im) + (2)**2*I0 + (2)*I1 + I2
    return(I_r,I_phi,I_t)

def delta_phi(I_phi,G_phi,lam):
    """
    Computes the net change of angle of the ray's trajectory due to gravitational lensing
    :param I_phi: radial phi path integral
    :param G_phi: angular phi path integral
    :param lam: angular momentum
    :return: net change in angle (geometric effect + frame dragging)
    """
    #Eqns (11 P3)
    return (-(I_phi+lam*G_phi)-np.pi/2)

#net change of time due to gravitational lensing
def delta_t(I_t,G_t,a):
    """
    Computes the total time ellapsed from the source to the observer along the ray's trajectory
    :param I_t: radial t path integral
    :param G_t: angular t path integral
    :param a: black hole spin
    :return: total time ellapsed
    """
    #Eqns (12 P3)
    return (-(I_t+a**2*G_t))

#calculate the source radius, angular, and time change for the rays associated with each grid. 
def calculate_observables(grid,mask,theta_o,a,mbar,distance=1000):
    """
    A master function to calculate all observable effects based on all above functions
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param theta_o: observer inclination
    :param a: black hole spin
    :param mbar: lensing band index 0,1,2,...
    :param distance: distance of the observer in units of black hole mass, set to 1000M by default
    :return rs: source radius
    :return redshift_sign: sign of the redshift
    :return deltat: time ellapsed
    :return deltaphi: source angular position
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]

    lam,eta = conserved_quantities(alpha,beta,theta_o,a)
    pm_o = np.sign(beta)
    r1,r2,r3,r4 = radial_turning_points(alpha,beta,lam,eta,a)
    mask2 = np.ones(r1.shape,dtype=bool)
    #Condition for turning points
    mask2[np.abs(r4.imag)>1e-13] = False
    mask3 = np.invert(mask2)
    u_p,u_m,theta_p,theta_m = angular_turning_points(alpha,beta,lam,eta,a)
    G_theta,G_phi,G_t = angular_integrals(mbar,beta,u_p,u_m,theta_p,theta_m,pm_o,theta_o,a)

    r31 = (r3-r1)
    r32 = (r3-r2)
    r41 = (r4-r1)
    r42 = (r4-r2)
    k = r32*r41/r31/r42
    taumax = np.zeros(alpha.shape)
    #Turn
    Jmax = 2/sqrt(r31[mask2]*r42[mask2])*ellipf(np.arcsin(np.sqrt(r31[mask2]/r41[mask2])).real, k[mask2].real)
    taumax[mask2]=Jmax.real
    
    r_sign0 = np.ones(taumax[mask2].shape)
    r_sign0[G_theta[mask2]>taumax[mask2]] = -1
    
    r_sign = np.ones(r1.shape)
    r_sign[mask2] = r_sign0
    redshift_sign = np.ones(mask.shape)
    redshift_sign[mask] = r_sign

    rs2 = source_radius2(distance,r1[mask2],r2[mask2],r3[mask2],r4[mask2],G_theta[mask2])
    rs3 = source_radius3(distance,r1[mask3],r2[mask3],r3[mask3],r4[mask3],G_theta[mask3])

    #rs2 = source_radius2_inf(r1[mask2],r2[mask2],r3[mask2],r4[mask2],G_theta[mask2])
    #rs3 = source_radius3_inf(r1[mask3],r2[mask3],r3[mask3],r4[mask3],G_theta[mask3])

    rs = np.zeros(mask.shape)
    r_mask = np.zeros(rs[mask].shape)
    r_mask[mask2] = rs2.real
    r_mask[mask3] = rs3.real
    rs[mask] = r_mask
    rs=np.nan_to_num(rs)
    r_p = 1+np.sqrt(1-a**2)
    rs[rs<=r_p] = r_p

    I_r,I_phi,I_t = radial_integrals(r_mask,distance,r1,r2,r3,r4,a,beta, mask2, mask3,lam,eta,r_sign)
    deltat = np.zeros(mask.shape)*np.nan
    deltaphi = np.zeros(mask.shape)*np.nan
    
    deltat[mask] = delta_t(I_t,G_t,a)
    deltaphi[mask] = delta_phi(I_phi,G_phi,lam)
    
    maskkk = np.ones(rs.shape)
    
    maskkk[rs<=r_p] = np.nan

    if imag_cut==1:
        maskkk[rs>=r_cutoff] = np.nan
    elif imag_cut==2:
        maskkk[rs>=r_cutoff+5] = np.nan

    return(rs*maskkk,redshift_sign*maskkk,deltat*maskkk,deltaphi*maskkk)

def rs_bv(grid,mask):
    """
    Beloborodov's approximation based on 0201117
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]

    b2=alpha**2+beta**2

    sign_bv=np.sign(-(beta*np.tan(thetao))/np.sqrt(b2+beta**2*np.tan(thetao)**2)-(1)/(1-np.sqrt(b2/3)*np.cos(1/3*np.arccos(-(3*np.sqrt(3))/(np.sqrt(b2))))))

    # This sign is only a nontrivial function outside the circle 3\sqrt{3} in Bardeen's coordinates
    sign_bv[np.isnan(sign_bv)]=1.0

    redshift_sign = np.ones(mask.shape)
    redshift_sign[mask]=sign_bv

    Y=-(beta*np.tan(thetao))/np.sqrt(b2+beta**2*np.tan(thetao)**2)
    rs0_bv=np.sqrt(((1-Y)/(1+Y))**2+b2/(1-Y**2))-(1-Y)/(1+Y)

    rs = np.zeros(mask.shape)
    rs[mask] = rs0_bv
    r_p = 2
    rs[rs<=r_p] = r_p

    maskkk = np.ones(rs.shape)
    
    maskkk[rs<=r_p] = np.nan

    print("Computing the n=0 based on Beloborodov's approximation")

    filename=path+"Rays_bv_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')
    h5f.create_dataset('rs0_bv', data=rs*maskkk)
    h5f.create_dataset('sign0_bv', data=redshift_sign)

    h5f.close()

def rt(supergrid0,mask0,supergrid1,mask1,supergrid2,mask2):
    """
    A master function to compute all observables for the n=0, 1, and 2 lensing bands
    :param supergrid0: alpha and beta grid for n=0 lensing band
    :param mask0: mask out the n=0 lensing band, see lb_f.py for detail
    :param supergrid1: alpha and beta grid for n=1 lensing band
    :param mask1: mask out the n=1 lensing band, see lb_f.py for detail
    :param supergrid2: alpha and beta grid for n=2 lensing band
    :param mask2: mask out the n=2 lensing band, see lb_f.py for detail
    
    :return: if save==1: save the results of calculate_observables as .h5 files for each lensing band.
             if live==1: output those results
    Note: the spin, observer inclination and distance can be changed in the input file.
    """

    print("Analytical ray-tracing of the n=0 band points")
    rs0,sign0,t0,phi0 = calculate_observables(supergrid0,mask0,thetao,spin_case,0,distance=D_obs)
    
    filename=path+"Rays_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('rs0', data=rs0)
    h5f.create_dataset('sign0', data=sign0)
    h5f.create_dataset('t0', data=t0)
    h5f.create_dataset('phi0', data=phi0)

    del rs0, sign0, t0, phi0, supergrid0, mask0

    print("Analytical ray-tracing of the n=1 band points")
    rs1,sign1,t1,phi1 = calculate_observables(supergrid1,mask1,thetao,spin_case,1,distance=D_obs)

    h5f.create_dataset('rs1', data=rs1)
    h5f.create_dataset('sign1', data=sign1)
    h5f.create_dataset('t1', data=t1)
    h5f.create_dataset('phi1', data=phi1)

    del rs1, sign1, t1, phi1, supergrid1, mask1

    print("Analytical ray-tracing of the n=2 band points")
    rs2,sign2,t2,phi2 = calculate_observables(supergrid2,mask2,thetao,spin_case,2,distance=D_obs)

    h5f.create_dataset('rs2', data=rs2)
    h5f.create_dataset('sign2', data=sign2)
    h5f.create_dataset('t2', data=t2)
    h5f.create_dataset('phi2', data=phi2)

    h5f.close()

    print("File ",filename," created.\n")