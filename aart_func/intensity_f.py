from aart_func import *
from params import * 

def gDisk(r,a,lamb):
    """
    Calculates the redshift factor for a photon outside the inner-most stable circular orbit(isco) (assume circular orbit)
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum

    :return: the redshift factor associated with the ray
    """
    #Eqns (P4 B7)
    return sqrt(r**3 - 3*r**2 + 2*a*r**(3/2))/(r**(3/2) - (lamb- a))

def Rint(r,a,lamb,eta):
    """
    Evaluates the "radial potential", for calculating the redshift factor for infalling material
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: carter constant

    :return: radial potential evaluated at the source
    """
    #Eqns (P2 5)
    return (r**2 + a**2 - a*lamb)**2 - (r**2 - 2*r + a**2)*(eta + (lamb - a)**2)

def gGas(r,b,a,lamb,eta):
    """
    Calculates the redshift factor for a photon inside the isco (assume infalling orbit)
    :param r: radius of the source
    :param b: sign for the redshift
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: carter constant

    :return: the redshift factor associated with the ray
    """
    #calculate radius of the inner-most stable circular orbit
    isco=rms(a)

    #Eqns (P2 2)
    Delta=r**2 - 2*r + a**2

    #Eqns (P3 B13)
    lambe=((isco**2 - 2*a*sqrt(isco) + a**2))/(isco**(3/2) - 2*sqrt(isco) + a)
    #Eqns (P3 B12)
    H=(2*r - a*lambe)/Delta

    #Eqns (P3 B14)
    gamma=sqrt(1 - 2/3 *1/isco)

    #Eqns (P3 B9-B11)
    ut=gamma*(1 + (2)/r *(1 + H))
    uphi=gamma/r**2*(lambe + a*H)
    ur=-np.sqrt(2/3/isco)*(isco/r - 1)**(3/2)
    #Eqns (P3 B15)
    return 1/(ut - uphi*lamb - ur*Delta**(-1)*b*sqrt(Rint(r,a,lamb,eta)))

#calculate the observed brightness for a purely radial profile
def bright_radial(grid,mask,redshift_sign,a,rs,isco,thetao):
    """
    Calculate the brightness of a rotationally symmetric disk
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param rs: source radius
    :param isco: radius of the inner-most stable circular orbit
    :param thetao: observer inclination
    :hidden param profile: radial profile of the intensity

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]
    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,lamb[rs>=isco])**gfactor*ilp.profile(rs[rs>=isco],a,gammap,mup,sigmap)

    brightness[rs<isco]= gGas(rs[rs<isco],redshift_sign[rs<isco],a,lamb[rs<isco],eta[rs<isco])**gfactor*ilp.profile(rs[rs<isco],a,gammap,mup,sigmap)
    
    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)

#calculate the observed brightness for an arbitrary profile, passed in as the interpolation object
#but ignoring the time delay due to lensing
def fast_light(grid,mask,redshift_sign,a,isco,rs,th,interpolation,thetao):
    """
    Calculate the black hole image ignoring the time delay due to lensing or geometric effect
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param interpolation: 2 dimensional brightness function of the source, interpolation object
    :param thetao: observer inclination

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    #Eqs. 60 and 61 in 0706.0622
    x_aux=np.sqrt(rs**2+a**2)*np.cos(th)
    y_aux=np.sqrt(rs**2+a**2)*np.sin(th)
 
    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,lamb[rs>=isco])**gfactor*interpolation(np.vstack([x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],redshift_sign[rs<isco],a,lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([x_aux[rs<isco],y_aux[rs<isco]]).T)
    
    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)

#calculate the observed brightness for an arbitrary, evolving profile, passed in as the interpolation object
def slow_light(grid,mask,redshift_sign,a,isco,rs,th,ts,interpolation,thetao):
    """
    Calculate the black hole image including the time delay due to lensing and geometric effect
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param ts: time of emission at the source
    :param interpolation: a time series of 2 dimensional brightness function of the source, 3d interpolation object
    :param thetao: observer inclination

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    ts = ts[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]
    
    #Eqs. 60 and 61 in 0706.0622
    x_aux=np.sqrt(rs**2+a**2)*np.cos(th)
    y_aux=np.sqrt(rs**2+a**2)*np.sin(th)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,lamb[rs>=isco])**gfactor*interpolation(np.vstack([ts[rs>=isco],x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],redshift_sign[rs<isco],a,lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([ts[rs<isco],x_aux[rs<isco],y_aux[rs<isco]]).T)

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)

def br(supergrid0,mask0,N0,rs0,sign0,supergrid1,mask1,N1,rs1,sign1,supergrid2,mask2,N2,rs2,sign2):
    """
    Calculate and save the radial brightness profile
    """
    bghts0 = bright_radial(supergrid0,mask0,sign0,spin_case,rs0,isco,thetao)
    bghts1 = bright_radial(supergrid1,mask1,sign1,spin_case,rs1,isco,thetao)
    bghts2 = bright_radial(supergrid2,mask2,sign2,spin_case,rs2,isco,thetao)

    I0 = bghts0.reshape(N0,N0).T
    I1 = bghts1.reshape(N1,N1).T
    I2 = bghts2.reshape(N2,N2).T

    filename=path+"Intensity_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('bghts0', data=I0)
    h5f.create_dataset('bghts1', data=I1)
    h5f.create_dataset('bghts2', data=I2)

    h5f.close()

    print("File ",filename," created.")

def br_bv(supergrid0,mask0,N0,rs0,sign0):
    """
    Calculate and save the radial brightness profile
    """
    bghts0 = bright_radial(supergrid0,mask0,sign0,spin_case,rs0,isco,thetao)

    I0 = bghts0.reshape(N0,N0).T

    filename=path+"Intensity_bv_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('bghts0', data=I0)

    h5f.close()

    print("File ",filename," created.")