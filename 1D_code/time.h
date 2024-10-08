#include <iostream>
#include <cmath>
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

void piece_Linear_coef_matrix(Eigen::SparseMatrix<double>& mat, int M_rows, int M_cols, int rows, int cols, double alpha_0, double beta_0, double h, double m_index, double coff);

void IMEX3_FS(Eigen::VectorXd& Uf_vec, double* Uf_0, double* Uf_1, double* Uf_2, Eigen::VectorXd& U_vec, Eigen::VectorXd& N_f, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, double* U_0, double* U_1, double* U_2, double tau, double a0, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver);
void IMEX3_SS(Eigen::VectorXd& Us_vec, double* Us_0, double* Us_1, double* Us_2, Eigen::VectorXd& N_s, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd N_f, double tau, double a0, double a1, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver);
void IMEX3_TS(Eigen::VectorXd& Ut_vec, double* Ut_0, double* Ut_1, double* Ut_2, Eigen::VectorXd& N_t, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd Us_vec, Eigen::VectorXd N_s, double tau, double a0, double b1, double b2, double a2, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver);

// time step
void compute_st(double* s_j, double* rho_0, double* m_0, double* Ee_0);
double compute_sj_max(double* s_j); double comput_rho_min(double* rho); double compute_dj(double rho_min);