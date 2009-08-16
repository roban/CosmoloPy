"""Routines related to the reionization of the IGM."""

import math

import numpy
import scipy
import scipy.integrate as si

import perturbation as cp
import distance as cd
import constants as cc
import density as cden

def delta_lambda_delta_dl(z, delta_dl, **cosmo):
    """Calculate Lyman-alpha wavelength shift given light-travel distance.

    Wavelengths are in Angstroms.

    Returns lambda(z), lambda(z - Deltaz), z, z - Deltaz

    """

    dl = cd.light_travel_distance(z, **cosmo)[0]
    dl0 = dl - delta_dl
    z0 = cd.redshift_d_light(dl0, **cosmo)
    
    wavelength = cc.lambda_Lya_0 * (1+z)
    wavelength0 = cc.lambda_Lya_0 * (1+z0)
    return wavelength, wavelength0, z, z0


def recomb_rate_coeff_HG(temp, species, case):
    """Recombination rate coefficients for HII, HeII and HeIII.

    Paramaters:
    ----------

    temp is the temperature in K

    species is 'H', 'He0', or 'He1'.

    case is 'A' or 'B'.

    Notes:
    -----
    
    From Hui and Gnedin (1997MNRAS.292...27H).

    Valid for He0 for temperatures between 5e3 and 5e5 K.

    """

    if not((species == 'H') or (species == 'He0') or (species == 'He1')):
        raise(exception.ValueError("Bad species: " + species))

    if case == 'A':
        case_N = 0
    elif case == 'B':
        case_N = 1
    else:
        raise(exception.ValueError("Case must be 'A' or 'B': " + case))

    # transition threshold temps in K
    T_TR = {'H'   : 157807.,
            'He0' : 285335.,
            'He1' : 631515.
            }        
    
    # cm**3 s**-1
    a = {'H'   : [1.269e-13, 2.753e-14],
         'He0' : [3e-14, 1.26e-14],
         'He1' : [2. * 1.269e-13, 2. * 2.753e-14]
         }
    
    p0 = {'H'   : [1.503, 1.500],
          'He0' : [0.654, 0.750],
          'He1' : [1.503, 1.500]
          }        

    p1 = {'H'   : [0.470, 0.407],
          'He0' : [0, 0],
          'He1' : [0.470, 0.407]
          }        

    p2 = {'H'   : [1.923,2.242],
          'He0' : [0, 0],
          'He1' : [1.923,2.242]
	  }        

    con = {'H'   : [0.522, 2.740],
         'He0' : [1, 1],
         'He1' : [0.522, 2.740],
         }        

    lam = 2. * T_TR[species] / temp
    return (a[species][case_N] * lam ** p0[species][case_N] /
            (1.0 + (lam / con[species][case_N]) ** p1[species][case_N]) 
            ** p2[species][case_N])

def ionization_from_collapse(z, coeff_ion, temp_min, passed_min_mass = False,
                             **cosmo):
    """Calculate the ionized fraction of the universe using
    perturbation theory.

    Parameters:
    ----------

    z: 
    
       Redshift values at which to evaluate the ionized fraction.

    coeff_ion:

       Coefficient giving the ratio between collapse fraction and
       ionized fraction (neglecting recombinations and assuming all
       photons are instantly absorbed).

    temp_min: 

       Either the minimum Virial temperature (in Kelvin) or minimum
       mass of halos (in solar masses) contributing to reionization.

    passed_temp_min: Boolean

       Set this to True if you pass a minimum mass, False (default) if
       you pass a minimum Virial temperature.

    **cosmo: dict

       Cosmological parameters.

    Notes:
    -----

    See Furlanetto et al. (2004ApJ...613....1F).

    """
    sd = cp.sig_del(temp_min, z, passed_min_mass=passed_min_mass, **cosmo)
    cf = cp.collapse_fraction(*sd)
    w = cf * coeff_ion
    return w

def quick_ion_col_function(coeff_ion, temp_min, passed_min_mass = False,
                           zmax = 20., zmin = 0., zstep = 0.1, **cosmo):
    """Return a function giving ionization_from_collapse as a function
    of redshift (based on interpolation).

    Calling the resulting function is much faster than evaluating
    ionization_from_collapse.
    """
    z = numpy.arange(zmin, zmax, zstep)
    w = ionization_from_collapse(z, coeff_ion, temp_min, 
                                 passed_min_mass, **cosmo)
    return scipy.interpolate.interp1d(z, w)

def clumping_factor_BKP(z):
    """Clumping factor as a function of redshift used by Bagla et al. 2009.

    See Bagla, Kulkarni & Padmanabhan (2009MNRAS.397..971B).
    """
    return numpy.sqrt(26.2917 * numpy.exp(-0.1822 * z + 0.003505 * z**2.))

def _udot_BKP(u, t, coeff_rec_func, redshift_func, ion_func):
    """du/dt where u = x - f_* f_esc,gamma N_gamma F
    
    Parameters:
    ----------

    u: integral of du/dt as defined below

    t: cosmic age in s

    redshift_func: function returning redshift given t
    
    ion_func: function returing ionized fraction neglecting recombinations

    coeff_rec_func: function returning clumping_factor^2 alpha_B n_H_0 (1+z)^3

    Notes:
    -----

    This is implemented as a reformulation of the ODE described in
    Bagla, Kulkarni & Padmanabhan (2009MNRAS.397..971B).

    The original ODE is:

    dx/dt = -alpha_B C^2 n_H x + f_* f_esc,gamma N_gamma dF/dt

    If we let u = x - w, where w = f_* f_esc,gamma N_gamma F(t) then

    du/dt = dx/dt - dw/dt

    which gives

    du/dt = -alpha_B C^2 n_H x = -alpha_B C^2 n_H (u + w)

    We have an analytical expression for w, so we can numerically
    integrate the ODE to give us u(t) or x(t) = u(t) + w(t).

    """
    z = redshift_func(t)
    w = ion_func(z)
    crf = coeff_rec_func(z)
    #ionization_from_collapse(z, coeff_ion, temp_min, 
    #                             passed_min_mass = passed_min_mass,
    #                             **cosmo)

    x = u + w
    x = x * (x <= 1.) + 1.0 * (x > 1.) 
    udot = -1. * crf * x
    #if (abs(round(z,1) - z) < 0.01):
    if (False):
        print ("z=%.3f; t=%.1g; c=%.2g; udot=%.2g; w,x,u = %.2g, %.2g, %.2g" % 
               (z, t, crf, udot, w, x, u))
    return udot

def integrate_ionization_BKP(z, coeff_ion,
                             temp_min = 1e4,
                             passed_min_mass = False,
                             temp_gas=1e4, 
                             alpha_B=None,
                             clump_fact_func = clumping_factor_BKP,
                             **cosmo):  

    """Integrate an ODE describing ionization and recombination rates.

    As described in Bagla, Kulkarni & Padmanabhan (hereafter BKP,
    2009MNRAS.397..971B).

    Parameters:
    ----------

    z: array 

       The redshift values at which to calculate the ionized
       fraction. This array should be in reverse numerical order. The
       first redshift specified should be early enough that the
       universe is still completely neutral.

    coeff_ion: 

       The coefficient converting the collapse fraction to ionized
       fraction, neglecting recombinations. Equivalent to the product
       (f_star * f_esc_gamma * N_gamma) in the BKP formalism.


    temp_min: 

       See docs for ionization_from_collapse. Either the minimum virial
       temperature or minimum mass of halos contributing to
       reionization.

    passed_temp_min: 

       See documentation for ionization_from_collapse.

    temp_gas: 

       Gas temperature used to calculate the recombination coefficient
       if alpha_b is not specified.

    alpha_B:

       Optional recombination coefficient in units of cm^3
       s^-1. In alpha_B=None, it is calculated from temp_gas.

    clump_fact_func: function

      Function returning the clumping factor when given a redshift.

   **cosmo: dict

      Dictionary specifying the cosmological parameters.

    Notes:
    -----
    
    BKP, as is fairly standard, assume that the ionized fraction is
    contained in fully ionized bubbles surrounded by a fully neutral
    IGM. The output is therefore the volume filling factor of ionized
    regions, not the ionized fraction of a uniformly ionized IGM.

    I have also made the standard assumption that all ionized photons
    are immediately absorbed, which allows the two differential
    equations to be combined into a single ODE. BKP say that this
    induces less than 5% error in the final value of the
    electron-scattering optical depth.

    I believe there is an extraneous factor of m_p in the second term
    of their equation 6, unless I am misunderstanding the definition
    of their variables. I have left it out.

    """

    # Determine recombination coefficient.
    if alpha_B is None:
        alpha_B_cm = recomb_rate_coeff_HG(temp_gas, 'H', 'B')
    else:
        alpha_B_cm = alpha_B
    alpha_B = alpha_B_cm / (cc.Mpc_cm**3.)
    print ("Recombination rate alpha_B = %.4g (Mpc^3 s^-1) = %.4g (cm^3 s^-1)" 
           % (alpha_B, alpha_B_cm))

    # Normalize power spectrum.
    if 'deltaSqr' not in cosmo:
        cosmo['deltaSqr'] = cp.norm_power(**cosmo)

    # Calculate useful densities.
    rho_crit, rho_0, n_He_0, n_H_0 = cden.baryon_densities(**cosmo)

    # Function used in the integration.
    # Units: (Mpc^3 s^-1) * Mpc^-3 = s^-1
    coeff_rec_func = lambda z: (clump_fact_func(z)**2. * 
                                alpha_B * 
                                n_H_0 * (1.+z)**3.)

    # Generate a function that converts redshift to age of the universe.
    redfunc, rerr1, rerr2 = \
        cd.quick_redshift_age_function(zmax = 1.1 * numpy.max(z), 
                                       zmin = -0.0,
                                       **cosmo)
    # Function used in the integration.
    ionfunc = quick_ion_col_function(coeff_ion, 
                                     temp_min, 
                                     passed_min_mass = passed_min_mass,
                                     zmax = 1.1 * numpy.max(z), 
                                     zmin = -0.0, 
                                     zstep = 0.1, **cosmo)

    # Convert specified redshifts to cosmic time (age of the universe).
    t, terr1, terr2 = cd.age(z, **cosmo)

    # Integrate to find u(z) = x(z) - w(z), where w is the ionization fraction 
    u = si.odeint(_udot_BKP, y0=0.0, t=t,
                  args=(coeff_rec_func, redfunc, ionfunc))
    u = u.flatten()

    w = ionization_from_collapse(z, coeff_ion, temp_min, 
                                 passed_min_mass = passed_min_mass,
                                 **cosmo)
    x = u + w
    x[x > 1.0] = 1.0
    return x, w, t
 
def integrate_optical_depth(x_ionH, x_ionHe, z, **cosmo):
    """Calculated optical depth due to electron scattering given a profile
    of the filling factor of ionized gas vs. redshift.

    Parameters:
    ----------
    
    x_ionH: array

       Ionized fraction of hydrogen as a function of z. Should be [0,1].

    x_ionHe: array 

       Set x_ionHE to X_HeII + 2 * X_HeIII, where X_HeII is the
       fraction of helium that is singly ionized, and X_HeII is the
       fraction of helium that is doubly ionized. See Notes below.
    
    z: array
       Redshift values at which the filling factor is specified.

    **cosmo: cosmological parameters
    
       uses: 'X_H' and/or 'Y_He', plus parameters needed for hubble_z

    Returns:
    -------

    tau: array
       The optical depth as a function of z.

    Notes:
    -----

    The precision of your result depends on the spacing of the input
    arrays. When in doubt, try doubling your z resolution and see if
    the optical depth values have converged.

    100% singly ionized helium means x_ionHe = 1.0, 100% doubly
    ionized helium means x_ionHe = 2.0

    If you want helium to be singly ionized at the same rate as
    hydrogen, set x_ionHe = x_ionH.

    If you want helium to be doubly ionized at the same rate as
    hydrogen is ionized, set x_ionHe = 2 * x_ionH.

    """

    rho_crit, rho_0, n_He_0, n_H_0 = cden.baryon_densities(**cosmo)

    # comoving Mpc^-1
    n_p = n_H_0 + 2. * n_He_0
    
    # comoving Mpc^-1
    n_e = n_H_0 * x_ionH + n_He_0 * x_ionHe

    # fraction of electrons that are free
    x = n_e / n_p

    H_0 = cc.H100_s * cosmo['h']

    # Mpc s^-1 * Mpc^2 * Mpc^-3 / s^-1 -> unitless
    tau_star = cc.c_light_Mpc_s * cc.sigma_T_Mpc * n_p

    # s^-1
    H_z = cd.hubble_z(z, **cosmo)

    # Mpc^3 s^-1 * Mpc^-3 / s^-1 -> unitless
    integrand = -1. * tau_star * x * ((1. + z)**2.) / H_z

    integral = si.cumtrapz(integrand, z)
    return numpy.abs(integral)

def optical_depth_instant(z_r, x_ionH=1.0, x_ionHe=2.0, 
                          return_tau_star=False, **cosmo):
    """Optical depth assuming instantaneous reionization and a flat
    universe.

    Calculates the optical depth due to Thompson scattering off free
    electrons in the IGM. 

    Parameters:
    ----------

    z_r: 
       Redshift of instantaneos reionization.

    x_ionH: 
       Ionized fraction of hydrogen after reionization.

    x_ionHe:
       value of X_HeII + 2 * X_HeIII after reionization (where X_HeII
       is the fraction of helium that is singly ionized, and X_HeII is
       the fraction of helium that is doubly ionized). Set to 2. for
       fully ionized helium.

    return_tau_star: Boolean
      whether or not to return the value of tau_star, as defined by
      Griffiths et al. (arxiv:astro-ph/9812125v3)

    **cosmo: cosmological parameters

    Returns:
    -------

    tau: array 
       optical depth to election

    tau_star: array or scalar

    Notes:
    -----

    See, e.g. Griffiths et al. (arxiv:astro-ph/9812125v3, note that
    the published version [ 1999MNRAS.308..854G] has typos)

    """

    if numpy.any(cden.get_omega_k_0(**cosmo) != 0):
        raise ValueError, "Not valid for non-flat (omega_k_0 !=0) cosmology."

    rho_crit, rho_0, n_He_0, n_H_0 = cden.baryon_densities(**cosmo)

    # comoving Mpc^-1
    n_p = n_H_0 + 2. * n_He_0
    
    # comoving Mpc^-1
    n_e = n_H_0 * x_ionH + n_He_0 * x_ionHe

    # fraction of electrons that are free
    x = n_e / n_p

    H_0 = cc.H100_s * cosmo['h']

    # Mpc s^-1 * Mpc^2 * Mpc^-3 / s^-1 -> unitless
    tau_star = cc.c_light_Mpc_s * cc.sigma_T_Mpc * n_p * x / H_0

    e_z_reion = cd.e_z(z_r, **cosmo)
    
    tau = 2. * tau_star * (e_z_reion - 1.0) / (3. * cosmo['omega_M_0'])

    if return_tau_star:
        return tau, tau_star
    else:
        return tau
    
