#include <iostream>
#include <cmath>
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
using namespace std;

#define N_x 40
#define N_y N_x
#define N_timestep 100
#define pi 3.14159265358979
#define Mode 2
#define CFL 0.1 / (2.0 * Mode + 1.0)
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
#define omega_ei 0.1
#define omega_er 0.1
#define kappa_e 0.5
#define kappa_i 0.5
#define kappa_r 0.5

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

// 2D Modal DG general framework
void get_Gaussian_points_x(double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* x, double h_x);
void get_Gaussian_points_y(double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double* y, double h_y);
void get_Gauss_quadrature_2D(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** u_11, double** u_21, double** u_31, double** u_41, double** u_51, double** u_12, double** u_22, double** u_32, double** u_42, double** u_52, double** u_13, double** u_23, double** u_33, double** u_43, double** u_53, double** u_14, double** u_24, double** u_34, double** u_44, double** u_54, double** u_15, double** u_25, double** u_35, double** u_45, double** u_55, double h_x, double h_y);
void get_L2_projection_2D(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double h_x, double h_y);

void createMatrix(int rows, int cols, double**& arr);
void deleteMatrix(int rows, double** arr);
void Matrix_to_vector(double* v, double** M, int M_rows, int M_cols);
void vector_to_Matrix(double** M, double* v, int M_rows, int M_cols);
void update_boundary(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5);
double compute_MAX(double** Matrix);

// accuracy solution
double rho_ex(double x, double y, double t); double u_ex(double x, double y, double t); double v_ex(double x, double y, double t);
double rhoe_e_ex(double x, double y, double t); double rhoe_i_ex(double x, double y, double t); double rhoe_r_ex(double x, double y, double t);
double rho_prime_ex(double x, double y, double t); double rho_prime2_ex(double x, double y, double t);


// model solution
double m1_ex(double x, double y, double t); double m2_ex(double x, double y, double t);
double Ee_ex(double x, double y, double t); double Ei_ex(double x, double y, double t); double Er_ex(double x, double y, double t);
double S3_ex(double x, double y, double t); double S4_ex(double x, double y, double t); double S5_ex(double x, double y, double t);


// Compute error and order
double* rho_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
double* m1_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
double* m2_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
double* Ee_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
double* Ei_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
double* Er_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t);
