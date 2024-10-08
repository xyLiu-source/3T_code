#include <iostream>
#include <cmath>

void OE_compute_u(double* u, double xi, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2);
void OE_compute_u_ave(double* u_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2);
void OE_compute_e_alpha(double* e_alpha, double xi, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2);
void OE_compute_e_alpha_ave(double* e_alpha_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2);
void OE_compute_Cs_ave(double* Cs_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2);
void OE_compute_beta_j(double* beta_j, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2);
void OE_compute_scalar_sigma(double* Linf, double* sigma_0, double* sigma_1, double* sigma_2, double* u_0, double* u_1, double* u_2);
void OE_compute_system_sigma(double* sigma_0, double* sigma_1, double* sigma_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2);
void OE_procedure(double* u_1, double* u_2, double* beta_j, double* sigma_0, double* sigma_1, double* sigma_2, double tau, double h_x);

