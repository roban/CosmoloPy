/* The following routines implement all of the fitting formulae in 
Eisenstein \& Hu (1997) */

/* These routines have been modified by Berian James to remove the
use of pointers. It is possible to pass pointers with numpy arrays
through SWIG, but it is more work that just reconfiguring the function
calls to be like those in power.c */

/* There are two sets of routines here.  The first set,

	TFfit_hmpc(), TFset_parameters(), and TFfit_onek(),

calculate the transfer function for an arbitrary CDM+baryon universe using
the fitting formula in Section 3 of the paper.  The second set,

	TFsound_horizon_fit(), TFk_peak(), TFnowiggles(), and TFzerobaryon(),

calculate other quantities given in Section 4 of the paper. */

#include <math.h>
#include <stdio.h>
void TFset_parameters(float omega0hh, float f_baryon, float Tcmb);

//float TFfit_onek(float k, float *tf_baryon, float *tf_cdm); 
float TFfit_onek(float k); 

void TFfit_hmpc(float omega0, float f_baryon, float hubble, float Tcmb,
	int numk, float *k, float *tf_full, float *tf_baryon, float *tf_cdm);

float TFsound_horizon_fit(float omega0, float f_baryon, float hubble);
float TFk_peak(float omega0, float f_baryon, float hubble);
float TFnowiggles(float omega0, float f_baryon, float hubble, 
		float Tcmb, float k_hmpc);
float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc);

/* ------------------------ DRIVER ROUTINE --------------------------- */
/* The following is an example of a driver routine you might use. */
/* Basically, the driver routine needs to call TFset_parameters() to
set all the scalar parameters, and then call TFfit_onek() for each 
wavenumber k you desire. */

/* While the routines use Mpc^-1 units internally, this driver has been
written to take an array of wavenumbers in units of h Mpc^-1.  On the
other hand, if you want to use Mpc^-1 externally, you can do this by
altering the variables you pass to the driver:  
	omega0 -> omega0*hubble*hubble, hubble -> 1.0 		*/

/* INPUT: omega0 -- the matter density (baryons+CDM) in units of critical 
	  f_baryon -- the ratio of baryon density to matter density 
	  hubble -- the Hubble constant, in units of 100 km/s/Mpc
	  Tcmb -- the CMB temperature in Kelvin. T<=0 uses the COBE value 2.728.
	  numk -- the length of the following zero-offset array
	  k[] -- the array of wavevectors k[0..numk-1]  */

/* INPUT/OUTPUT: There are three output arrays of transfer functions. 
All are zero-offset and, if used, must have storage [0..numk-1] declared
in the calling program.  However, if you substitute the NULL pointer for
one or more of the arrays, then that particular transfer function won't
be outputted. The transfer functions are:

	tf_full[] -- The full fitting formula, eq. (16), for the matter
			transfer function. 
	tf_baryon[] -- The baryonic piece of the full fitting formula, eq. 21.
	tf_cdm[] -- The CDM piece of the full fitting formula, eq. 17. */

/* Again, you can set these pointers to NULL in the function call if
you don't want a particular output. */

/* Various intermediate scalar quantities are stored in global variables, 
so that you might more easily access them.  However, this also means that
you would be better off not simply #include'ing this file in your programs,
but rather compiling it separately, calling only the driver, and using
extern declarations to access the intermediate quantities. */

/* ----------------------------- DRIVER ------------------------------- */

void TFfit_hmpc(float omega0, float f_baryon, float hubble, float Tcmb,
	int numk, float *k, float *tf_full, float *tf_baryon, float *tf_cdm)
/* Remember: k[0..numk-1] is in units of h Mpc^-1. */
{
    int j;
    float tf_thisk, baryon_piece, cdm_piece;

    TFset_parameters(omega0*hubble*hubble, f_baryon, Tcmb);

    for (j=0;j<numk;j++) {
      //tf_thisk = TFfit_onek(k[j]*hubble, &baryon_piece, &cdm_piece); 
      tf_thisk = TFfit_onek(k[j]*hubble); 
      if (tf_full!=NULL) tf_full[j] = tf_thisk;
      //if (tf_baryon!=NULL) tf_baryon[j] = baryon_piece;
      //if (tf_cdm!=NULL) tf_cdm[j] = cdm_piece;
    }
    return;
}

/* ------------------------ FITTING FORMULAE ROUTINES ----------------- */

/* There are two routines here.  TFset_parameters() sets all the scalar
parameters, while TFfit_onek() calculates the transfer function for a 
given wavenumber k.  TFfit_onek() may be called many times after a single
call to TFset_parameters() */

/* Global variables -- We've left many of the intermediate results as 
global variables in case you wish to access them, e.g. by declaring
them as extern variables in your main program. */
/* Note that all internal scales are in Mpc, without any Hubble constants! */

float	omhh,		/* Omega_matter*h^2 */
	obhh,		/* Omega_baryon*h^2 */
	theta_cmb,	/* Tcmb in units of 2.7 K */
	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
	k_equality,	/* Scale of equality, in Mpc^-1 */
	z_drag,		/* Redshift of drag epoch */
	R_drag,		/* Photon-baryon ratio at drag epoch */
	R_equality,	/* Photon-baryon ratio at equality epoch */
	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
	k_silk,		/* Silk damping scale, in Mpc^-1 */
	alpha_c,	/* CDM suppression */
	beta_c,		/* CDM log shift */
	alpha_b,	/* Baryon suppression */
	beta_b,		/* Baryon envelope shift */
	beta_node,	/* Sound horizon shift */
	k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
	sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
        alpha_gamma,	/* Gamma suppression in approximate TF */
        T_full; /* EDIT: Make transfer function value global */

/* Convenience from Numerical Recipes in C, 2nd edition */
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static float cubearg;
#define CUBE(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)
static float pow4arg;
#define POW4(a) ((pow4arg=(a)) == 0.0 ? 0.0 : pow4arg*pow4arg*pow4arg*pow4arg)
	/* Yes, I know the last one isn't optimal; it doesn't appear much */

void TFset_parameters(float omega0hh, float f_baryon, float Tcmb)
/* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
/* Input: omega0hh -- The density of CDM and baryons, in units of critical dens,
		multiplied by the square of the Hubble constant, in units
		of 100 km/s/Mpc */
/* 	  f_baryon -- The fraction of baryons to CDM */
/*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
			of the COBE value of  2.728 K. */
/* Output: Nothing, but set many global variables used in TFfit_onek(). 
You can access them yourself, if you want. */
/* Note: Units are always Mpc, never h^-1 Mpc. */
{
    float z_drag_b1, z_drag_b2;
    float alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;

    if (f_baryon<=0.0 || omega0hh<=0.0) {
	fprintf(stderr, "TFset_parameters(): Illegal input.\n");
	exit(1);
    }
    omhh = omega0hh;
    obhh = omhh*f_baryon;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;

    z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*
		(1+z_drag_b1*pow(obhh,z_drag_b2));
    
    R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
    R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);

    sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
	    log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));

    k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));

    alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    alpha_c = pow(alpha_c_a1,-f_baryon)*
		pow(alpha_c_a2,-CUBE(f_baryon));
    
    beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    beta_c_b2 = pow(0.395*omhh, -0.0266);
    beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));

    y = z_equality/(1+z_drag);
    alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;

    beta_node = 8.41*pow(omhh, 0.435);
    beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);

    k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));

    alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
		SQR(f_baryon);
    
    return;
}

//float TFfit_onek(float k, float *tf_baryon, float *tf_cdm)
float TFfit_onek(float k)
/* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
	  *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
				the input was not NULL. */
/* Output: Returns the value of the full transfer function fitting formula.
		This is the form given in Section 3 of Eisenstein & Hu (1997).
  	  *tf_baryon -- The baryonic contribution to the full fit.
	  *tf_cdm -- The CDM contribution to the full fit. */
/* Notes: Units are Mpc, not h^-1 Mpc. */
{
    float T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
    float q, xx, xx_tilde, q_eff;
    //float T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
    float T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon; // Define T_full as global
    float T_0_L0, T_0_C0, T_0, gamma_eff; 
    float T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;

    k = fabs(k);	/* Just define negative k as positive */
    //if (k==0.0) {
      //if (tf_baryon!=NULL) *tf_baryon = 1.0;
      //if (tf_cdm!=NULL) *tf_cdm = 1.0;
    //	return 1.0;
    //}

    q = k/13.41/k_equality;
    xx = k*sound_horizon;

    T_c_ln_beta = log(2.718282+1.8*beta_c*q);
    T_c_ln_nobeta = log(2.718282+1.8*q);
    T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
    T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));

    T_c_f = 1.0/(1.0+POW4(xx/5.4));
    T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
	    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));
    
    s_tilde = sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
    xx_tilde = k*s_tilde;

    T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
    T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2))+
		alpha_b/(1+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
    
    f_baryon = obhh/omhh;
    T_full = f_baryon*T_b + (1-f_baryon)*T_c;

    /* Now to store these transfer functions */
    //if (tf_baryon!=NULL) *tf_baryon = T_b;
    //if (tf_cdm!=NULL) *tf_cdm = T_c;
    return T_full;
}

/* ======================= Approximate forms =========================== */

float TFsound_horizon_fit(float omega0, float f_baryon, float hubble)
/* Input: omega0 -- CDM density, in units of critical density
	  f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
	  hubble -- Hubble constant, in units of 100 km/s/Mpc
/* Output: The approximate value of the sound horizon, in h^-1 Mpc. */
/* Note: If you prefer to have the answer in  units of Mpc, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    float omhh, sound_horizon_fit_mpc;
    omhh = omega0*hubble*hubble;
    sound_horizon_fit_mpc = 
	44.5*log(9.83/omhh)/sqrt(1+10.0*pow(omhh*f_baryon,0.75));
    return sound_horizon_fit_mpc*hubble;
}

float TFk_peak(float omega0, float f_baryon, float hubble)
/* Input: omega0 -- CDM density, in units of critical density
	  f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
	  hubble -- Hubble constant, in units of 100 km/s/Mpc
/* Output: The approximate location of the first baryonic peak, in h Mpc^-1 */
/* Note: If you prefer to have the answer in  units of Mpc^-1, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    float omhh, k_peak_mpc;
    omhh = omega0*hubble*hubble;
    k_peak_mpc = 2.5*3.14159*(1+0.217*omhh)/TFsound_horizon_fit(omhh,f_baryon,1.0);
    return k_peak_mpc/hubble;
}

float TFnowiggles(float omega0, float f_baryon, float hubble, 
		float Tcmb, float k_hmpc)
/* Input: omega0 -- CDM density, in units of critical density
	  f_baryon -- Baryon fraction, the ratio of baryon to CDM density.
	  hubble -- Hubble constant, in units of 100 km/s/Mpc
	  Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
			COBE FIRAS value of 2.728 K
	  k_hmpc -- Wavenumber in units of (h Mpc^-1). */
/* Output: The value of an approximate transfer function that captures the
non-oscillatory part of a partial baryon transfer function.  In other words,
the baryon oscillations are left out, but the suppression of power below
the sound horizon is included. See equations (30) and (31).  */
/* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    float k, omhh, theta_cmb, k_equality, q, xx, alpha_gamma, gamma_eff;
    float q_eff, T_nowiggles_L0, T_nowiggles_C0;

    k = k_hmpc*hubble;	/* Convert to Mpc^-1 */
    omhh = omega0*hubble*hubble;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;

    k_equality = 0.0746*omhh/SQR(theta_cmb);
    q = k/13.41/k_equality;
    xx = k*TFsound_horizon_fit(omhh, f_baryon, 1.0);

    alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
		SQR(f_baryon);
    gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+POW4(0.43*xx)));
    q_eff = q*omhh/gamma_eff;

    T_nowiggles_L0 = log(2.0*2.718282+1.8*q_eff);
    T_nowiggles_C0 = 14.2 + 731.0/(1+62.5*q_eff);
    return T_nowiggles_L0/(T_nowiggles_L0+T_nowiggles_C0*SQR(q_eff));
}

/* ======================= Zero Baryon Formula =========================== */

float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc)
/* Input: omega0 -- CDM density, in units of critical density
	  hubble -- Hubble constant, in units of 100 km/s/Mpc
	  Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
			COBE FIRAS value of 2.728 K
	  k_hmpc -- Wavenumber in units of (h Mpc^-1). */
/* Output: The value of the transfer function for a zero-baryon universe. */
/* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
and omega0 -> omega0*hubble^2. */ 
{
    float k, omhh, theta_cmb, k_equality, q, T_0_L0, T_0_C0;

    k = k_hmpc*hubble;	/* Convert to Mpc^-1 */
    omhh = omega0*hubble*hubble;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;

    k_equality = 0.0746*omhh/SQR(theta_cmb);
    q = k/13.41/k_equality;

    T_0_L0 = log(2.0*2.718282+1.8*q);
    T_0_C0 = 14.2 + 731.0/(1+62.5*q);
    return T_0_L0/(T_0_L0+T_0_C0*q*q);
}
