#include <iostream>
#include <cmath>

void F_d(double* u, double* u_0, double* u_1, double* u_2, double xi, int N_index);
void C_s(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void compute_alphaC(double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2);
void compute_uhalf(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x, int index_mode);
void compute_coef(double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x, int index_mode);
void compute_jump(double* u_jump, double* u_plus, double* u_minus, int N_index);
void compute_ave(double* u_ave, double* u_plus, double* u_minus, int N_index); void compute_uxi(double* u, double* u_1, double* u_2, double xi);
double Compute_alpha_max(double* alpha_c); double Compute_alpha_min(double* alpha_c);

// Convection part
void Fd_hat(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* alpha_c);
void convection(double* U_0, double* U_1, double* U_2, double* F_hat, double* GP_0, double* GP_1, double* GP_2);
void convection_d(double* con_d_0, double* con_d_1, double* con_d_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* alpha_c, double h_x);
void F_m(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void Fm_hat(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double* alpha_c);
void convection_m(double* con_m_0, double* con_m_1, double* con_m_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double* alpha_c, double h_x);
void F_alpha(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double xi);
void Falpha_hat(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double* alpha_c);
void convection_alpha(double* con_alpha_0, double* con_alpha_1, double* con_alpha_2, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double* alpha_c, double h_x);

// Nonconservetive part
void p_alpha(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double xi);
void p_alphaxi(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_1, double* Ealpha_2, double xi);
void F_Ne(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi);
void F_Ni(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi);
void F_Nr(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi);
void palpha_jump(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2);
void Nonconservetive_Ne(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode);
void Nonconservetive_Ni(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode);
void Nonconservetive_Nr(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode);

// Source part
void compute_Tr4(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void compute_Talpha(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi);
void F_Si(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double xi);
void F_Sr(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void Source_i(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double h_x);
void Source_r(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Er_0, double* Er_1, double* Er_2, double h_x);
void Source_e(double* U_0, double* U_1, double* U_2, double* Si_0, double* Si_1, double* Si_2, double* Sr_0, double* Sr_1, double* Sr_2);

// Diffusion part
void compute_Talphaxi(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi);
void compute_Tr4xi(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void compute_Talphaxi2(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi);
void compute_Tr4xi2(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi);
void Talpha_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha);
void Talphaxi_average(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha);
void Talphaxi2_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha);
void Tr4_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2);
void Tr4xi_average(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2);
void Tr4xi2_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2);
void Diffusion_alpha(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double kappa_alpha, double diff_alpha, double h_x);
void Diffusion_r(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double diff_alpha, double h_x);

void space_d(double* U_0, double* U_1, double* U_2, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x);
void space_m(double* U_0, double* U_1, double* U_2, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x);
void space_e(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode);
void space_i(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode);
void space_r(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode);

