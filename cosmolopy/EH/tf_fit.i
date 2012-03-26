/* tf_fit.i */
%module tf_fit
%{
  /* Put header files here or function declarations like below */
  extern void TFset_parameters(float omega0hh, float f_baryon, float Tcmb);
  extern float TFfit_onek(float k);
  extern float T_full;
  %}

extern void TFset_parameters(float omega0hh, float f_baryon, float Tcmb);
extern float TFfit_onek(float k);
extern float T_full;
