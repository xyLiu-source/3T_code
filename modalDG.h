#include <iostream>
#include <cmath>
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
using namespace std;

#define N_x 80
#define N_timestep 100
#define pi 3.14159265358979
#define Mode 2
#define CFL 0.2 / (2.0 * Mode + 1.0)
#define T_final 1.0
#define OE_index 1
#define Omega 1.0

// parameters of model
#define gamma_e 5.0/3.0
#define gamma_i 5.0/3.0
#define gamma_r 4.0/3.0
#define C_ve 1.0
#define C_vi 1.0
#define a 1.0
#define omega_ei 1.0
#define omega_er 1.0
#define kappa_e 1.0
#define kappa_i 1.0
#define kappa_r 1.0

#define Coef_d2 4.0
#define Coef_d1 2.0
#define diff_beta 1.0/12.0
#define gamma 0.435866521508459

// Quadrature points and weights
#define GQxi_1 -0.906179845938664
#define weight_1 0.236926885056189

#define GQxi_2 -0.538469310105683
#define weight_2 0.478628670499366

#define GQxi_3 0.0
#define weight_3 0.568888888888889

#define GQxi_4 0.538469310105683
#define weight_4 0.478628670499366

#define GQxi_5 0.906179845938664
#define weight_5 0.236926885056189


// Modal DG general framework
void get_Gaussian_points(double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* x, int N_index, double h_x);
void get_L2_projection(double* U_0, double* U_1, double* U_2, double* u_0, double* u_1, double* u_2, int N_index, double h_x);
void get_Gauss_quadrature(double* U_0, double* U_1, double* U_2, double* u_x1, double* u_x2, double* u_x3, double* u_x4, double* u_x5, int N_index, double h_x);
void get_Gauss_quadrature_prime(double* U_0, double* U_1, double* U_2, double* u_0, double* u_1, int N_index, double h_x);
void get_time_prime(double* U_0, double* U_1, double* U_2, double* S_0, double* S_1, double* S_2, int N_index, double h_x);
void Mode_select(double* u_0, double* u_1, double* u_2, int index_mode, int N_index);
void update_boundary(double* u_0, double* u_1, double* u_2);

// accuracy solution
double rho_ex(double x, double t); double u_ex(double x, double t);
double rhoe_e_ex(double x, double t); double rhoe_i_ex(double x, double t); double rhoe_r_ex(double x, double t);
double rho_prime_ex(double x, double t); double rho_prime2_ex(double x, double t);

// model solution
double m_ex(double x, double t); double Ee_ex(double x, double t); double Ei_ex(double x, double t); double Er_ex(double x, double t);
double S3_ex(double x, double t); double S4_ex(double x, double t); double S5_ex(double x, double t);

// Compute error and order
double* rho_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t);
double* m_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t);
double* Ee_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t);
double* Ei_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t);
double* Er_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t);