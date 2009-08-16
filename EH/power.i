/* power.i */
%module power
%{
  /* Put header files here or function declarations like below */
  extern int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
		     int degen_hdm, float omega_lambda, float hubble, 
		     float redshift);
  extern float TFmdm_onek_mpc(float     kk);
  extern float TFmdm_onek_hmpc(float kk );
  extern float tf_cbnu;
  %}
 
extern int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
		   int degen_hdm, float omega_lambda, float hubble, 
		   float redshift);
extern float TFmdm_onek_mpc(float     kk);
extern float TFmdm_onek_hmpc(float kk );

extern float tf_cbnu;
