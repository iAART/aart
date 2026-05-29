from aart_func import *
from params import * 

def Delta(r,a):
    """
    Calculates the Kerr metric function \Delta(t)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return r**2-2*r+a**2

def PIF(r,a):
    """
    Calculates PI(r) (Eq. B6 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (r**2+a**2)**2-a**2*Delta(r,a)

def urbar(r,a):
    """
    Calculates the r (contravariant) component of the four velocity for radial infall
    (Eq. B34b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return -np.sqrt(2*r*(r**2+a**2))/(r**2)

def Omegabar(r,a):
    """
    Calculates the angular velocity of the radial infall
    (Eq. B32a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (2*a*r)/PIF(r,a)

def Omegahat(r,a,laux):
    """
    Calculates the angular velocity of the sub-Keplerian orbit
    (Eq. B39 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (a+(1-2/r)*(laux-a))/(PIF(r,a)/(r**2)-(2*a*laux)/r)

def uttilde(r, a,urT,OT):
    """
    Calculates the t (contravariant) component of the general four velocity
    (Eq. B52 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param urT: r (contravariant) component of the general four velocity
    :param OT: Angular velocity of the general four velocity
    """
    return np.sqrt((1 + urT**2*r**2/Delta(r,a))/(1-(r**2+a**2)*OT**2-(2/r)*(1-a*OT)**2))

def Ehat(r,a,laux):
    """
    Calculates the orbital energy of the sub-Keplerian flow
    (Eq. B44a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    """
    return np.sqrt(Delta(r,a)/(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2))

def nuhat(r,a,laux,Ehataux):
    """
    Calculates the radial velocity of the sub-Keplerian flow
    (Eq. B45 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    :param Ehataux: sub-Keplerian orbital energy
    """
    return r/Delta(r,a)*np.sqrt(np.abs(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2-Delta(r,a)/(Ehataux**2)))

def lhat(r,a):
    """
    Calculates the rspecific angular momentum of the sub-Keplerian flow
    (Eq. B44b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return sub_kep*(r**2+a**2-2*a*np.sqrt(r))/(np.sqrt(r)*(r-2)+a)

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

def gDisk(r,a,b,lamb,eta):
    """
    Calculates the redshift factor for a photon outside the inner-most stable circular orbit(isco) (assume circular orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param b: The +- sign of p^r
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: Carter constant

    :return: the redshift factor associated with the ray
    """

    OH=Omegahat(r,a,lhat(r,a))
    OT=OH+(1-betaphi)*(Omegabar(r,a)-OH)
    ur=(1-betar)*urbar(r,a)
    ut=uttilde(r,a,ur,OT)
    uphi=ut*OT
    
    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))

def gGas(r,a,b,lamb,eta):
    """
    Calculates the redshift factor for a photon inside the isco (assume infalling orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param b: sign for the redshift
    :param lamb: angular momentum
    :param eta: carter constant

    :return: the redshift factor associated with the ray
    """
    #Calculate radius of the inner-most stable circular orbit
    isco=rms(a)

    lms=lhat(isco,a)
    OH=Omegahat(r,a,lms)
    OT=OH+(1-betaphi)*(Omegabar(r,a)-OH)

    Ems=Ehat(isco,a,lms)
    urhat=-Delta(r,a)/(r**2)*nuhat(r, a, lms ,Ems)*Ems
    ur=urhat+(1-betar)*(urbar(r,a)-urhat)
    ut=uttilde(r,a,ur,OT)
    uphi=OT*ut

    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))

#TODO: This expression will just work for sure for the Keplerian velocity. 
# I need to check if it has to be modified for the general four-velocity
def CosAng(r,a,b,lamb,eta):
    """
    Calculates the cosine of the emission angle
    :param r: radius of the source
    :param a: spin of the black hole
    :param b: sign for the redshift
    :param lamb: angular momentum
    :param eta: Carter constant

    :return: the  cosine of the emission angle
    """
    #From eta, solve for Sqrt(p_\theta/p_t)
    kthkt=np.sqrt(eta)
    #Sqrt(g^{\theta\theta}) Evaluated at the equatorial plane
    thth=1/r
    return thth*gDisk(r,a,b,lamb,eta)*kthkt

#calculate the observed brightness for a purely radial profile
def bright_radial(grid,mask,redshift_sign,a,rs,isco,thetao):
    """
    Calculate the brightness of a rotationally symmetric disk
    (Eq. 50 P1)
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param rs: source radius
    :param isco: radius of the inner-most stable circular orbit
    :param thetao: observer inclination

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]

    rs = rs[mask]

    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)

    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*ilp.profile(rs[rs>=isco],a,gammap,mup,sigmap)
    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*ilp.profile(rs[rs<isco],a,gammap,mup,sigmap)
    
    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    
    return(I)


#calculate the observed brightness for an arbitrary profile, passed in as the interpolation object
#but ignoring the time delay due to lensing
def fast_light(grid,mask,redshift_sign,a,isco,rs,th,ts,interpolation,thetao):
    """
    Calculate the black hole image ignoring the time delay due to lensing or geometric effect
    (Eq. 116 P1)
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

    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    #Added by Daniel for comparisons
    K_1 =np.count_nonzero(rs>=isco)
    tcol_1 = np.full(K_1, ts)

    K_2 =np.count_nonzero(rs<isco)
    tcol_2 = np.full(K_2, ts)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*interpolation(np.vstack([tcol_1,x_aux[rs>=isco],y_aux[rs>=isco]]).T)

    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([tcol_2,x_aux[rs<isco],y_aux[rs<isco]]).T)
    
    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)

#calculate the observed brightness for an arbitrary, evolving profile, passed in as the interpolation object
def slow_light(grid,mask,redshift_sign,a,isco,rs,th,ts,interpolation,thetao):
    """
    Calculate the black hole image including the time delay due to lensing and geometric effect
    (Eq. 50 P1)

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
    
    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*interpolation(np.vstack([ts[rs>=isco],x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([ts[rs<isco],x_aux[rs<isco],y_aux[rs<isco]]).T)

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

def gfactorf(grid,mask,redshift_sign,a,isco,rs,thetao):
    """
    Calculate the redshift factor
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param thetao: observer inclination

    :return: redshift factor at each point.

    """
    
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    gfact = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]
    
    gfact[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])
    gfact[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])
    
    r_p = 1+np.sqrt(1-a**2)
    gfact[rs<=r_p] = 0
    
    gs = np.zeros(mask.shape)
    gs[mask] = gfact
    return(gs)

# orbit for the centroid with radhs=Radius of the hotspot and velhs = 0.01 (angular frequency)
# one may put an arbitrary orbit
def x0(t):
    return(radhs*np.cos(t*velhs))

def y0(t):
    return(radhs*np.sin(t*velhs))

def flare_model(grid,mask,redshift_sign,a,rs,th,ts,thetao,rwidth,delta_t):

    """
    Calculate the black hole image including the time delay due to lensing and geometric effect
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param mbar: lensing band index 0,1,2,...
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
    
    x_aux = rs*np.cos(th)
    y_aux = rs*np.sin(th)
    # x0 and y0 is now a function of t, where one can specify an arbitrary equitorial orbit
    brightness = np.exp(-(x_aux-x0(ts+delta_t))**2/rwidth**2-(y_aux-y0(ts+delta_t))**2/rwidth**2)
    
    brightness[rs>=isco]*= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor
    brightness[rs<isco]*= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(np.nan_to_num(I))

#BRISK LIGHT


def clip_ts_to_interval(ts, left_s, right_s):
    """
    Clip finite values of a 1D time array to an interval on a periodic domain.

    Case 1: left_s <= right_s
        Valid interval is [left_s, right_s].
        Values below are mapped to left_s, values above to right_s.

    Case 2: left_s > right_s
        The valid interval wraps around the period boundary:
            [left_s, T) U [0, right_s]
        Then the excluded gap is (right_s, left_s).
        Values inside that gap are mapped to the nearest boundary,
        using the midpoint of the gap.

    Non-finite values (NaN, inf, -inf) are preserved.

    :param ts array_like, 1D: Array of emission times.
    :param left_s float: Left boundary.
    :param right_s float: Right boundary.

    :Returns ts_out : ndarray Output array with clipped values.
    """
    ts = np.asarray(ts, dtype=float)

    if ts.ndim != 1:
        raise ValueError("ts must be a 1D array.")
    if ts.size == 0:
        raise ValueError("ts cannot be empty.")

    ts_out = ts.copy()
    finite_mask = np.isfinite(ts_out)

    if left_s <= right_s:
        # Standard interval: [left_s, right_s]
        ts_out[finite_mask] = np.clip(ts_out[finite_mask], left_s, right_s)

    else:
        # Wrapped interval: [left_s, T) U [0, right_s]
        # Gap to collapse: (right_s, left_s)
        center = 0.5 * (left_s + right_s)

        gap_mask = finite_mask & (ts_out > right_s) & (ts_out < left_s)

        left_half_mask = gap_mask & (ts_out > center)
        right_half_mask = gap_mask & (ts_out < center)

        ts_out[left_half_mask] = left_s
        ts_out[right_half_mask] = right_s

        # If a value is exactly equal to center, choose one side consistently.
        ts_out[gap_mask & (ts_out == center)] = right_s

    return ts_out


def modal_hdi_kde(data,p,*,bounds=None,trim_quantiles=(0.005, 0.995),fit_on="trimmed",bw_method=None,gridsize=4096):
    """
    KDE modal mass interval with robust handling of long tails / extreme values.

    Parameters
    ----------
    data : array-like
        1D data.
    p : float
        Desired probability mass around the mode, between 0 and 1.
        p=0 returns only the mode.
        p=1 returns the whole accepted support.
    bounds : tuple or None
        Explicit bounds `(lower, upper)` of the region you want to consider.
        If None, bounds are computed from `trim_quantiles`.
    trim_quantiles : tuple
        Quantiles used to define the accepted support when `bounds=None`.
        Example: (0.005, 0.995) ignores the lowest 0.5% and highest 0.5%.
        For a one-sided long upper tail, use e.g. (0.0, 0.99).
    fit_on : {"trimmed", "winsorized", "all"}
        How to fit the KDE:
        - "trimmed": remove points outside bounds before fitting. Recommended.
        - "winsorized": clip points to bounds before fitting.
        - "all": fit KDE to all data, but compute mass only inside bounds.
    bw_method : str, float, callable, optional
        Passed to scipy.stats.gaussian_kde.
    gridsize : int
        Number of grid points used for numerical approximation.

    Returns
    -------
    dict
        Contains the KDE, mode, interval, selected mass, grid, density, mask,
        accepted bounds, and fraction of data ignored.
    """

    x = np.asarray(data, dtype=float).ravel()
    x = x[np.isfinite(x)]

    if x.size < 2:
        raise ValueError("Need at least two finite data points.")

    if not 0 <= p <= 1:
        raise ValueError("p must satisfy 0 <= p <= 1.")

    # Define the accepted region.
    if bounds is None:
        q_low, q_high = trim_quantiles
        if not 0 <= q_low < q_high <= 1:
            raise ValueError("trim_quantiles must satisfy 0 <= low < high <= 1.")
        lo, hi = np.quantile(x, [q_low, q_high])
    else:
        lo, hi = map(float, bounds)

    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("Invalid bounds.")

    ignored_fraction = np.mean((x < lo) | (x > hi))

    # Choose which data are used to fit the KDE.
    if fit_on == "trimmed":
        x_fit = x[(x >= lo) & (x <= hi)]
    elif fit_on == "winsorized":
        x_fit = np.clip(x, lo, hi)
    elif fit_on == "all":
        x_fit = x
    else:
        raise ValueError("fit_on must be 'trimmed', 'winsorized', or 'all'.")

    if x_fit.size < 2:
        raise ValueError("Too few data points remain after trimming.")

    if np.ptp(x_fit) == 0:
        raise ValueError("Remaining data are constant; KDE is not well-defined.")

    kde = gaussian_kde(x_fit, bw_method=bw_method)

    # Evaluate KDE only on the accepted support.
    x_grid = np.linspace(lo, hi, gridsize)
    density_raw = kde(x_grid)

    # Renormalize density on accepted support.
    total_area = np.trapz(density_raw, x_grid)

    if total_area <= 0 or not np.isfinite(total_area):
        raise ValueError("Could not normalize KDE on the accepted support.")

    density = density_raw / total_area

    mode_idx = np.argmax(density)
    mode = x_grid[mode_idx]

    if p == 0:
        mask = np.zeros_like(x_grid, dtype=bool)
        mask[mode_idx] = True

        return {
            "kde": kde,
            "mode": mode,
            "interval": (mode, mode),
            "mass": 0.0,
            "bounds": (lo, hi),
            "ignored_fraction": ignored_fraction,
            "x_grid": x_grid,
            "density": density,
            "mask": mask,
            "density_threshold": density[mode_idx]}

    if p == 1:
        mask = np.ones_like(x_grid, dtype=bool)

        return {
            "kde": kde,
            "mode": mode,
            "interval": (lo, hi),
            "mass": 1.0,
            "bounds": (lo, hi),
            "ignored_fraction": ignored_fraction,
            "x_grid": x_grid,
            "density": density,
            "mask": mask,
            "density_threshold": 0.0,
        }

    def modal_component_at_threshold(threshold):
        """Connected region around the mode where density >= threshold."""
        above = density >= threshold

        left = mode_idx
        while left > 0 and above[left - 1]:
            left -= 1

        right = mode_idx
        while right < len(x_grid) - 1 and above[right + 1]:
            right += 1

        return left, right

    def component_mass(left, right):
        if right <= left:
            return 0.0
        return np.trapz(density[left:right + 1], x_grid[left:right + 1])

    # Binary search for the highest density threshold whose modal component
    # contains probability mass at least p.
    low = 0.0
    high = density[mode_idx]

    for _ in range(60):
        mid = 0.5 * (low + high)
        left, right = modal_component_at_threshold(mid)
        mass_mid = component_mass(left, right)

        if mass_mid >= p:
            low = mid
        else:
            high = mid

    threshold = low
    left, right = modal_component_at_threshold(threshold)

    interval = (x_grid[left], x_grid[right])
    mass = component_mass(left, right)

    mask = np.zeros_like(x_grid, dtype=bool)
    mask[left:right + 1] = True

    return {
        "kde": kde,
        "mode": mode,
        "interval": interval,
        "mass": mass,
        "bounds": (lo, hi),
        "ignored_fraction": ignored_fraction,
        "x_grid": x_grid,
        "density": density,
        "mask": mask,
        "density_threshold": threshold,
    }

def brisk_light(grid, mask, redshift_sign, a, isco, rs, th, ts,interpolation, thetao, left_s, right_s):
    """
    Calculate the black hole image including the time delay due to lensing and geometric effect but with a restriction in the source to the range [left_s, right_s]

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
    :param left_s: left boundary of the modal HDI
    :param right_s: right boundary of the modal HDI 

    :return: image of a lensed equitorial source with only radial dependence. 
    """

    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    ts = ts[mask]

    ts_reduce = clip_ts_to_interval(ts, left_s,right_s)
    
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]
    
    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*interpolation(np.vstack([ts_reduce[rs>=isco],x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([ts_reduce[rs<isco],x_aux[rs<isco],y_aux[rs<isco]]).T)

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0
    
    I = np.zeros(mask.shape) 
    I[mask] = brightness
    return(I)

def plot_kde_mass_result(res, data=None, bins=100, show_hist=True):
    """
    Plot KDE, shaded modal mass region, mode, and accepted bounds.

    Parameters
    ----------
    res : dict
        Output from robust_kde_mass_around_mode or kde_mass_around_mode.
    data : array-like, optional
        Original data. Used only to show a histogram.
    bins : int
        Number of histogram bins.
    show_hist : bool
        Whether to show histogram of the data.
    """

    x_grid = res["x_grid"]
    density = res["density"]
    mask = res["mask"]

    fig, ax = plt.subplots(figsize=(8, 4))

    # Optional histogram
    if show_hist and data is not None:
        data = np.asarray(data)
        data = data[np.isfinite(data)]

        # If robust bounds exist, plot only data inside accepted support
        if "bounds" in res:
            lo, hi = res["bounds"]
            data = data[(data >= lo) & (data <= hi)]

        ax.hist(
            data,
            bins=bins,
            density=True,
            alpha=0.25,
            label="Data histogram"
        )

    # KDE curve
    ax.plot(x_grid, density, label="KDE")

    # Shaded mass around the mode
    ax.fill_between(
        x_grid,
        density,
        where=mask,
        alpha=0.35,
        label="Mass around mode"
    )

    # Mode
    ax.axvline(
        res["mode"],
        linestyle="--",
        label=f"Mode = {res['mode']:.4g}"
    )

    # Selected interval
    left, right = res["interval"]
    ax.axvline(left, linestyle=":", label=f"Interval = [{left:.4g}, {right:.4g}]")
    ax.axvline(right, linestyle=":")

    # Accepted bounds, if using the robust version
    if "bounds" in res:
        blo, bhi = res["bounds"]
        ax.axvline(blo, linestyle="-.", alpha=0.7, label="Accepted bounds")
        ax.axvline(bhi, linestyle="-.", alpha=0.7)

    #llo,lhi=res['interval']
    #ax.set_xlim(-200,0)

    ax.set_xlabel("x")
    ax.set_ylabel("Density")
    ax.set_title("KDE mass around the mode")
    ax.legend()
    plt.tight_layout()
    plt.show()