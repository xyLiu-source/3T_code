#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include <string.h>
#include "modalDG.h"
#include "space.h"
#include "time.h"
#include "OE.h"
using namespace std;


int main()
{
	// start timing
	std::time_t start, end;
	start = std::time(NULL);

	// spatial discretization
	double x_min = 0.0, x_max = 2 * pi;
	double h_x = (x_max - x_min) / N_x;

	double* x = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		x[i] = x_min - 0.5 * h_x + i * h_x;
	}

	// temporal discretization
	double t_ini = 0.0;
	
	// initialization
	double t = t_ini;
	//L2-projection
	double* X_1 = new double[N_x + 2]; double* X_2 = new double[N_x + 2]; double* X_3 = new double[N_x + 2];
	double* X_4 = new double[N_x + 2]; double* X_5 = new double[N_x + 2];
	get_Gaussian_points(X_1, X_2, X_3, X_4, X_5, x, N_x + 2, h_x);

	double* rho_x1 = new double[N_x + 2]; double* rho_x2 = new double[N_x + 2]; double* rho_x3 = new double[N_x + 2]; double* rho_x4 = new double[N_x + 2]; double* rho_x5 = new double[N_x + 2];
	double* m_x1 = new double[N_x + 2]; double* m_x2 = new double[N_x + 2]; double* m_x3 = new double[N_x + 2]; double* m_x4 = new double[N_x + 2]; double* m_x5 = new double[N_x + 2];
	double* Ee_x1 = new double[N_x + 2]; double* Ee_x2 = new double[N_x + 2]; double* Ee_x3 = new double[N_x + 2]; double* Ee_x4 = new double[N_x + 2]; double* Ee_x5 = new double[N_x + 2];
	double* Ei_x1 = new double[N_x + 2]; double* Ei_x2 = new double[N_x + 2]; double* Ei_x3 = new double[N_x + 2]; double* Ei_x4 = new double[N_x + 2]; double* Ei_x5 = new double[N_x + 2];
	double* Er_x1 = new double[N_x + 2]; double* Er_x2 = new double[N_x + 2]; double* Er_x3 = new double[N_x + 2]; double* Er_x4 = new double[N_x + 2]; double* Er_x5 = new double[N_x + 2];

	for (int i = 0; i < N_x + 2; i++)
	{
		rho_x1[i] = rho_ex(X_1[i], t); rho_x2[i] = rho_ex(X_2[i], t); rho_x3[i] = rho_ex(X_3[i], t);
		rho_x4[i] = rho_ex(X_4[i], t); rho_x5[i] = rho_ex(X_5[i], t);
		m_x1[i] = m_ex(X_1[i], t); m_x2[i] = m_ex(X_2[i], t); m_x3[i] = m_ex(X_3[i], t);
		m_x4[i] = m_ex(X_4[i], t); m_x5[i] = m_ex(X_5[i], t);
		Ee_x1[i] = Ee_ex(X_1[i], t); Ee_x2[i] = Ee_ex(X_2[i], t); Ee_x3[i] = Ee_ex(X_3[i], t);
		Ee_x4[i] = Ee_ex(X_4[i], t); Ee_x5[i] = Ee_ex(X_5[i], t);
		Ei_x1[i] = Ei_ex(X_1[i], t); Ei_x2[i] = Ei_ex(X_2[i], t); Ei_x3[i] = Ei_ex(X_3[i], t);
		Ei_x4[i] = Ei_ex(X_4[i], t); Ei_x5[i] = Ei_ex(X_5[i], t);
		Er_x1[i] = Er_ex(X_1[i], t); Er_x2[i] = Er_ex(X_2[i], t); Er_x3[i] = Er_ex(X_3[i], t);
		Er_x4[i] = Er_ex(X_4[i], t); Er_x5[i] = Er_ex(X_5[i], t);
	}
	double* GQ_0 = new double[N_x + 2]; double* GQ_1 = new double[N_x + 2]; double* GQ_2 = new double[N_x + 2];
	double* rho_0 = new double[N_x + 2]; double* rho_1 = new double[N_x + 2]; double* rho_2 = new double[N_x + 2];
	double* m_0 = new double[N_x + 2]; double* m_1 = new double[N_x + 2]; double* m_2 = new double[N_x + 2];
	double* Ee_0 = new double[N_x + 2]; double* Ee_1 = new double[N_x + 2]; double* Ee_2 = new double[N_x + 2];
	double* Ei_0 = new double[N_x + 2]; double* Ei_1 = new double[N_x + 2]; double* Ei_2 = new double[N_x + 2];
	double* Er_0 = new double[N_x + 2]; double* Er_1 = new double[N_x + 2]; double* Er_2 = new double[N_x + 2];

	get_Gauss_quadrature(GQ_0, GQ_1, GQ_2, rho_x1, rho_x2, rho_x3, rho_x4, rho_x5, N_x + 2, h_x);
	get_L2_projection(rho_0, rho_1, rho_2, GQ_0, GQ_1, GQ_2, N_x + 2, h_x);

	get_Gauss_quadrature(GQ_0, GQ_1, GQ_2, m_x1, m_x2, m_x3, m_x4, m_x5, N_x + 2, h_x);
	get_L2_projection(m_0, m_1, m_2, GQ_0, GQ_1, GQ_2, N_x + 2, h_x);

	get_Gauss_quadrature(GQ_0, GQ_1, GQ_2, Ee_x1, Ee_x2, Ee_x3, Ee_x4, Ee_x5, N_x + 2, h_x);
	get_L2_projection(Ee_0, Ee_1, Ee_2, GQ_0, GQ_1, GQ_2, N_x + 2, h_x);

	get_Gauss_quadrature(GQ_0, GQ_1, GQ_2, Ei_x1, Ei_x2, Ei_x3, Ei_x4, Ei_x5, N_x + 2, h_x);
	get_L2_projection(Ei_0, Ei_1, Ei_2, GQ_0, GQ_1, GQ_2, N_x + 2, h_x);

	get_Gauss_quadrature(GQ_0, GQ_1, GQ_2, Er_x1, Er_x2, Er_x3, Er_x4, Er_x5, N_x + 2, h_x);
	get_L2_projection(Er_0, Er_1, Er_2, GQ_0, GQ_1, GQ_2, N_x + 2, h_x);

	delete[] rho_x1; delete[] rho_x2; delete[] rho_x3; delete[] rho_x4; delete[] rho_x5;
	delete[] m_x1; delete[] m_x2; delete[] m_x3; delete[] m_x4; delete[] m_x5;
	delete[] Ee_x1; delete[] Ee_x2; delete[] Ee_x3; delete[] Ee_x4; delete[] Ee_x5;
	delete[] Ei_x1; delete[] Ei_x2; delete[] Ei_x3; delete[] Ei_x4; delete[] Ei_x5;
	delete[] Er_x1; delete[] Er_x2; delete[] Er_x3; delete[] Er_x4; delete[] Er_x5;
	delete[] GQ_0; delete[] GQ_1; delete[] GQ_2;

	// determining time step
	double* s_j = new double[N_x];
	compute_st(s_j, rho_0, m_0, Ee_0);
	double s_max = compute_sj_max(s_j);
	double tau_s = CFL / s_max;

	double* beta_j = new double[N_x];
	OE_compute_beta_j(beta_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2);
	double beta_max = compute_sj_max(beta_j);
	double tau_a = CFL * h_x / beta_max;

	double tau = tau_a;
	if (tau_s < tau)
	{
		tau = tau_s;
	}

	cout << tau_a << endl;
	cout << tau_s << endl;
	cout << tau << endl;

	double t_fin = T_final;

	// mass matrix
	double alpha_0 = Coef_d2; double beta_0 = diff_beta;

	int mrows = N_x; int mcols = N_x;
	int Mrows = (Mode + 1) * mrows; int Mcols = (Mode + 1) * mcols;
	Eigen::SparseMatrix<double> M(Mrows, Mcols);
	piece_Linear_coef_matrix(M, Mrows, Mcols, mrows, mcols, alpha_0, beta_0, h_x, 0.0, 1.0);


	// IMEX3 
	double b1 = -3.0 / 2.0 * pow(gamma, 2) + 4.0 * gamma - 1.0 / 4.0;
	double b2 = 3.0 / 2.0 * pow(gamma, 2) - 5.0 * gamma + 5.0 / 4.0;
	double a1 = -0.35;
	double a2 = (1.0 / 3.0 - 2.0 * pow(gamma, 2) - 2.0 * b2 * a1 * gamma) / (gamma * (1.0 - gamma));

	double a0 = kappa_r;
	double coff = -1.0 * tau * gamma * a0;
	Eigen::SparseMatrix<double> A(Mrows, Mcols);
	piece_Linear_coef_matrix(A, Mrows, Mcols, mrows, mcols, alpha_0, beta_0, h_x, 1.0, coff);

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	double* E = new double[N_x + 2]; double* E_new = new double[N_x + 2];
	double rho_c = 0.0; double m_c = 0.0; double E_c = 0.0;
	// time stepping
	cout << "Method: Modal DG with P" << Mode << " element. Mesh: " << N_x << endl;
	cout << "Start computing the numerical solution ... " << endl;
	int Timestep_index = 0;
	double Tau = 10.0;

	while (t < t_fin)
	{
		if (t_fin - t < tau)
		{
			tau = t_fin - t;
			coff = -1.0 * tau * gamma * a0;
			piece_Linear_coef_matrix(A, Mrows, Mcols, mrows, mcols, alpha_0, beta_0, h_x, 1.0, coff);
			solver.compute(A);
		}

		Timestep_index = Timestep_index + 1;
		if (fmod(Timestep_index, N_timestep) == 0)
		{
			compute_st(s_j, rho_0, m_0, Ee_0);
			s_max = compute_sj_max(s_j);
			tau_s = CFL / s_max;
			OE_compute_beta_j(beta_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2);
			beta_max = compute_sj_max(beta_j);
			tau_a = CFL * h_x / beta_max;

			if (tau_s < tau_a)
			{
				Tau = tau_s;
			}
			else
			{
				Tau = tau_a;
			}

			if (Tau < tau)
			{
				tau = Tau;
			}

			coff = -1.0 * tau * gamma * a0;
			piece_Linear_coef_matrix(A, Mrows, Mcols, mrows, mcols, alpha_0, beta_0, h_x, 1.0, coff);
			solver.compute(A);
		}

		for (int i = 0; i < 9; i++)
		{
			if (t > T_final / 10.0 * (i + 1) - tau / 2.0 && t < T_final / 10.0 * (i + 1) + tau / 2.0)
			{
				int j = i + 1;
				cout << "t = " << T_final / 10.0 * j << endl;
			}
		}
		
		for (int i = 0; i < N_x + 2; i++)
		{
			E[i] = Ee_0[i] + Ei_0[i] + Er_0[i];
		}

		double* Sigma_0 = new double[N_x]; double* Sigma_1 = new double[N_x]; double* Sigma_2 = new double[N_x];
		double* Beta_j = new double[N_x];
		double* alpha_c = new double[N_x + 1];
		double* coef_L = new double[N_x]; double* coef_R = new double[N_x];

		double* S_0 = new double[N_x + 2]; double* S_1 = new double[N_x + 2]; double* S_2 = new double[N_x + 2];

		double t_1 = t; double t_2 = t + gamma * tau; 
		double t_3 = t + (1.0 + gamma) / 2.0 * tau; double t_4 = t + tau;
		

		// First step of IMEX3
		Mode_select(rho_0, rho_1, rho_2, Mode, N_x + 2); Mode_select(m_0, m_1, m_2, Mode, N_x + 2); Mode_select(Ee_0, Ee_1, Ee_2, Mode, N_x + 2);
		Mode_select(Ei_0, Ei_1, Ei_2, Mode, N_x + 2); Mode_select(Er_0, Er_1, Er_2, Mode, N_x + 2);

		double* rhof_0 = new double[N_x + 2]; double* rhof_1 = new double[N_x + 2]; double* rhof_2 = new double[N_x + 2];
		double* mf_0 = new double[N_x + 2]; double* mf_1 = new double[N_x + 2]; double* mf_2 = new double[N_x + 2];
		double* Eef_0 = new double[N_x + 2]; double* Eef_1 = new double[N_x + 2]; double* Eef_2 = new double[N_x + 2];
		double* Eif_0 = new double[N_x + 2]; double* Eif_1 = new double[N_x + 2]; double* Eif_2 = new double[N_x + 2];
		double* Erf_0 = new double[N_x + 2]; double* Erf_1 = new double[N_x + 2]; double* Erf_2 = new double[N_x + 2];

		compute_alphaC(alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2);
		compute_coef(coef_L, coef_R, rho_0, rho_1, rho_2, m_0, m_1, m_2, h_x, Mode);

		if (t > t_fin - 5 * tau)
		{
			double alphac_min = Compute_alpha_min(alpha_c);
			double alphac_max = Compute_alpha_max(alpha_c);
			cout << alphac_min << " , " << alphac_max << endl;
		}

		space_d(S_0, S_1, S_2, alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, h_x);
		double* rhofTP_0 = new double[N_x + 2]; double* rhofTP_1 = new double[N_x + 2]; double* rhofTP_2 = new double[N_x + 2];
		get_time_prime(rhofTP_0, rhofTP_1, rhofTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			rhof_0[i] = rho_0[i] + tau * gamma * rhofTP_0[i];
			rhof_1[i] = rho_1[i] + tau * gamma * rhofTP_1[i];
			rhof_2[i] = rho_2[i] + tau * gamma * rhofTP_2[i];
		}
		
		space_m(S_0, S_1, S_2, alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x);
		double* mfTP_0 = new double[N_x + 2]; double* mfTP_1 = new double[N_x + 2]; double* mfTP_2 = new double[N_x + 2];
		get_time_prime(mfTP_0, mfTP_1, mfTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			mf_0[i] = m_0[i] + tau * gamma * mfTP_0[i];
			mf_1[i] = m_1[i] + tau * gamma * mfTP_1[i];
			mf_2[i] = m_2[i] + tau * gamma * mfTP_2[i];
		}

		space_e(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, t_1, Mode);
		double* EefTP_0 = new double[N_x + 2]; double* EefTP_1 = new double[N_x + 2]; double* EefTP_2 = new double[N_x + 2];
		get_time_prime(EefTP_0, EefTP_1, EefTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Ee_vec(Mrows); Eigen::VectorXd Eef_vec(Mrows); Eigen::VectorXd EeN_f(Mrows);
		IMEX3_FS(Eef_vec, Eef_0, Eef_1, Eef_2, Ee_vec, EeN_f, M, mrows, Mrows, EefTP_0, EefTP_1, EefTP_2, Ee_0, Ee_1, Ee_2, tau, a0, solver);

		space_i(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, t_1, Mode);
		double* EifTP_0 = new double[N_x + 2]; double* EifTP_1 = new double[N_x + 2]; double* EifTP_2 = new double[N_x + 2];
		get_time_prime(EifTP_0, EifTP_1, EifTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Ei_vec(Mrows); Eigen::VectorXd Eif_vec(Mrows); Eigen::VectorXd EiN_f(Mrows);
		IMEX3_FS(Eif_vec, Eif_0, Eif_1, Eif_2, Ei_vec, EiN_f, M, mrows, Mrows, EifTP_0, EifTP_1, EifTP_2, Ei_0, Ei_1, Ei_2, tau, a0, solver);

		space_r(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, t_1, Mode);
		double* ErfTP_0 = new double[N_x + 2]; double* ErfTP_1 = new double[N_x + 2]; double* ErfTP_2 = new double[N_x + 2];
		get_time_prime(ErfTP_0, ErfTP_1, ErfTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Er_vec(Mrows); Eigen::VectorXd Erf_vec(Mrows); Eigen::VectorXd ErN_f(Mrows);
		IMEX3_FS(Erf_vec, Erf_0, Erf_1, Erf_2, Er_vec, ErN_f, M, mrows, Mrows, ErfTP_0, ErfTP_1, ErfTP_2, Er_0, Er_1, Er_2, tau, a0, solver);

		// update boundary
		update_boundary(rhof_0, rhof_1, rhof_2); update_boundary(mf_0, mf_1, mf_2);
		update_boundary(Eef_0, Eef_1, Eef_2); update_boundary(Eif_0, Eif_1, Eif_2); update_boundary(Erf_0, Erf_1, Erf_2);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(Sigma_0, Sigma_1, Sigma_2, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2);
			OE_compute_beta_j(Beta_j, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2);

			//double tau_f = t_2 - t_1;
			double tau_f = tau;

			OE_procedure(rhof_1, rhof_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_f, h_x);
			OE_procedure(mf_1, mf_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_f, h_x);
			OE_procedure(Eef_1, Eef_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_f, h_x);
			OE_procedure(Eif_1, Eif_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_f, h_x);
			OE_procedure(Erf_1, Erf_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_f, h_x);

			for (int j = 0; j < mrows; j++)
			{
				Eef_vec(j) = Eef_0[j + 1]; Eef_vec(mrows + j) = Eef_1[j + 1]; Eef_vec(2 * mrows + j) = Eef_2[j + 1];
				Eif_vec(j) = Eif_0[j + 1]; Eif_vec(mrows + j) = Eif_1[j + 1]; Eif_vec(2 * mrows + j) = Eif_2[j + 1];
				Erf_vec(j) = Erf_0[j + 1]; Erf_vec(mrows + j) = Erf_1[j + 1]; Erf_vec(2 * mrows + j) = Erf_2[j + 1];
			}

			// update boundary
			update_boundary(rhof_0, rhof_1, rhof_2); update_boundary(mf_0, mf_1, mf_2);
			update_boundary(Eef_0, Eef_1, Eef_2); update_boundary(Eif_0, Eif_1, Eif_2); update_boundary(Erf_0, Erf_1, Erf_2);
		}
		
		// Second step of IMEX3
		Mode_select(rhof_0, rhof_1, rhof_2, Mode, N_x + 2); Mode_select(mf_0, mf_1, mf_2, Mode, N_x + 2); Mode_select(Eef_0, Eef_1, Eef_2, Mode, N_x + 2);
		Mode_select(Eif_0, Eif_1, Eif_2, Mode, N_x + 2); Mode_select(Erf_0, Erf_1, Erf_2, Mode, N_x + 2);

		double* rhos_0 = new double[N_x + 2]; double* rhos_1 = new double[N_x + 2]; double* rhos_2 = new double[N_x + 2];
		double* ms_0 = new double[N_x + 2]; double* ms_1 = new double[N_x + 2]; double* ms_2 = new double[N_x + 2];
		double* Ees_0 = new double[N_x + 2]; double* Ees_1 = new double[N_x + 2]; double* Ees_2 = new double[N_x + 2];
		double* Eis_0 = new double[N_x + 2]; double* Eis_1 = new double[N_x + 2]; double* Eis_2 = new double[N_x + 2];
		double* Ers_0 = new double[N_x + 2]; double* Ers_1 = new double[N_x + 2]; double* Ers_2 = new double[N_x + 2];

		compute_alphaC(alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2);
		compute_coef(coef_L, coef_R, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, h_x, Mode);
		

		space_d(S_0, S_1, S_2, alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, h_x);
		double* rhosTP_0 = new double[N_x + 2]; double* rhosTP_1 = new double[N_x + 2]; double* rhosTP_2 = new double[N_x + 2];
		get_time_prime(rhosTP_0, rhosTP_1, rhosTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			rhos_0[i] = rho_0[i] + tau * ((1.0 + gamma) / 2.0 - a1) * rhofTP_0[i] + tau * a1 * rhosTP_0[i];
			rhos_1[i] = rho_1[i] + tau * ((1.0 + gamma) / 2.0 - a1) * rhofTP_1[i] + tau * a1 * rhosTP_1[i];
			rhos_2[i] = rho_2[i] + tau * ((1.0 + gamma) / 2.0 - a1) * rhofTP_2[i] + tau * a1 * rhosTP_2[i];
		}

		space_m(S_0, S_1, S_2, alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2, h_x);
		double* msTP_0 = new double[N_x + 2]; double* msTP_1 = new double[N_x + 2]; double* msTP_2 = new double[N_x + 2];
		get_time_prime(msTP_0, msTP_1, msTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			ms_0[i] = m_0[i] + tau * ((1.0 + gamma) / 2.0 - a1) * mfTP_0[i] + tau * a1 * msTP_0[i];
			ms_1[i] = m_1[i] + tau * ((1.0 + gamma) / 2.0 - a1) * mfTP_1[i] + tau * a1 * msTP_1[i];
			ms_2[i] = m_2[i] + tau * ((1.0 + gamma) / 2.0 - a1) * mfTP_2[i] + tau * a1 * msTP_2[i];
		}

		space_e(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2, h_x, t_2, Mode);
		double* EesTP_0 = new double[N_x + 2]; double* EesTP_1 = new double[N_x + 2]; double* EesTP_2 = new double[N_x + 2];
		get_time_prime(EesTP_0, EesTP_1, EesTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Ees_vec(Mrows); Eigen::VectorXd EeN_s(Mrows);
		IMEX3_SS(Ees_vec, Ees_0, Ees_1, Ees_2, EeN_s, M, mrows, Mrows, EesTP_0, EesTP_1, EesTP_2, Ee_vec, Eef_vec, EeN_f, tau, a0, a1, solver);

		space_i(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2, h_x, t_2, Mode);
		double* EisTP_0 = new double[N_x + 2]; double* EisTP_1 = new double[N_x + 2]; double* EisTP_2 = new double[N_x + 2];
		get_time_prime(EisTP_0, EisTP_1, EisTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Eis_vec(Mrows); Eigen::VectorXd EiN_s(Mrows);
		IMEX3_SS(Eis_vec, Eis_0, Eis_1, Eis_2, EiN_s, M, mrows, Mrows, EisTP_0, EisTP_1, EisTP_2, Ei_vec, Eif_vec, EiN_f, tau, a0, a1, solver);

		space_r(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhof_0, rhof_1, rhof_2, mf_0, mf_1, mf_2, Eef_0, Eef_1, Eef_2, Eif_0, Eif_1, Eif_2, Erf_0, Erf_1, Erf_2, h_x, t_2, Mode);
		double* ErsTP_0 = new double[N_x + 2]; double* ErsTP_1 = new double[N_x + 2]; double* ErsTP_2 = new double[N_x + 2];
		get_time_prime(ErsTP_0, ErsTP_1, ErsTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Ers_vec(Mrows); Eigen::VectorXd ErN_s(Mrows);
		IMEX3_SS(Ers_vec, Ers_0, Ers_1, Ers_2, ErN_s, M, mrows, Mrows, ErsTP_0, ErsTP_1, ErsTP_2, Er_vec, Erf_vec, ErN_f, tau, a0, a1, solver);

		// update boundary
		update_boundary(rhos_0, rhos_1, rhos_2); update_boundary(ms_0, ms_1, ms_2);
		update_boundary(Ees_0, Ees_1, Ees_2); update_boundary(Eis_0, Eis_1, Eis_2); update_boundary(Ers_0, Ers_1, Ers_2);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(Sigma_0, Sigma_1, Sigma_2, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2);
			OE_compute_beta_j(Beta_j, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2);

			//double tau_s = t_3 - t_1;
			double tau_s = tau;

			OE_procedure(rhos_1, rhos_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_s, h_x);
			OE_procedure(ms_1, ms_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_s, h_x);
			OE_procedure(Ees_1, Ees_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_s, h_x);
			OE_procedure(Eis_1, Eis_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_s, h_x);
			OE_procedure(Ers_1, Ers_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau_s, h_x);

			for (int j = 0; j < mrows; j++)
			{
				Ees_vec(j) = Ees_0[j + 1]; Ees_vec(mrows + j) = Ees_1[j + 1]; Ees_vec(2 * mrows + j) = Ees_2[j + 1];
				Eis_vec(j) = Eis_0[j + 1]; Eis_vec(mrows + j) = Eis_1[j + 1]; Eis_vec(2 * mrows + j) = Eis_2[j + 1];
				Ers_vec(j) = Ers_0[j + 1]; Ers_vec(mrows + j) = Ers_1[j + 1]; Ers_vec(2 * mrows + j) = Ers_2[j + 1];
			}

			// update boundary
			update_boundary(rhos_0, rhos_1, rhos_2); update_boundary(ms_0, ms_1, ms_2);
			update_boundary(Ees_0, Ees_1, Ees_2); update_boundary(Eis_0, Eis_1, Eis_2); update_boundary(Ers_0, Ers_1, Ers_2);
		}


		// Third step of IMEX3
		Mode_select(rhos_0, rhos_1, rhos_2, Mode, N_x + 2); Mode_select(ms_0, ms_1, ms_2, Mode, N_x + 2); Mode_select(Ees_0, Ees_1, Ees_2, Mode, N_x + 2);
		Mode_select(Eis_0, Eis_1, Eis_2, Mode, N_x + 2); Mode_select(Ers_0, Ers_1, Ers_2, Mode, N_x + 2);

		double* rhot_0 = new double[N_x + 2]; double* rhot_1 = new double[N_x + 2]; double* rhot_2 = new double[N_x + 2];
		double* mt_0 = new double[N_x + 2]; double* mt_1 = new double[N_x + 2]; double* mt_2 = new double[N_x + 2];
		double* Eet_0 = new double[N_x + 2]; double* Eet_1 = new double[N_x + 2]; double* Eet_2 = new double[N_x + 2];
		double* Eit_0 = new double[N_x + 2]; double* Eit_1 = new double[N_x + 2]; double* Eit_2 = new double[N_x + 2];
		double* Ert_0 = new double[N_x + 2]; double* Ert_1 = new double[N_x + 2]; double* Ert_2 = new double[N_x + 2];

		compute_alphaC(alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2);
		compute_coef(coef_L, coef_R, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, h_x, Mode);


		space_d(S_0, S_1, S_2, alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, h_x);
		double* rhotTP_0 = new double[N_x + 2]; double* rhotTP_1 = new double[N_x + 2]; double* rhotTP_2 = new double[N_x + 2];
		get_time_prime(rhotTP_0, rhotTP_1, rhotTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			rhot_0[i] = rho_0[i] + tau * (1.0 - a2) * rhosTP_0[i] + tau * a2 * rhotTP_0[i];
			rhot_1[i] = rho_1[i] + tau * (1.0 - a2) * rhosTP_1[i] + tau * a2 * rhotTP_1[i];
			rhot_2[i] = rho_2[i] + tau * (1.0 - a2) * rhosTP_2[i] + tau * a2 * rhotTP_2[i];
		}

		space_m(S_0, S_1, S_2, alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2, h_x);
		double* mtTP_0 = new double[N_x + 2]; double* mtTP_1 = new double[N_x + 2]; double* mtTP_2 = new double[N_x + 2];
		get_time_prime(mtTP_0, mtTP_1, mtTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			mt_0[i] = m_0[i] + tau * (1.0 - a2) * msTP_0[i] + tau * a2 * mtTP_0[i];
			mt_1[i] = m_1[i] + tau * (1.0 - a2) * msTP_1[i] + tau * a2 * mtTP_1[i];
			mt_2[i] = m_2[i] + tau * (1.0 - a2) * msTP_2[i] + tau * a2 * mtTP_2[i];
		}
		
		space_e(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2, h_x, t_3, Mode);
		double* EetTP_0 = new double[N_x + 2]; double* EetTP_1 = new double[N_x + 2]; double* EetTP_2 = new double[N_x + 2];
		get_time_prime(EetTP_0, EetTP_1, EetTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Eet_vec(Mrows); Eigen::VectorXd EeN_t(Mrows);
		IMEX3_TS(Eet_vec, Eet_0, Eet_1, Eet_2, EeN_t, M, mrows, Mrows, EetTP_0, EetTP_1, EetTP_2, Ee_vec, Eef_vec, Ees_vec, EeN_s, tau, a0, b1, b2, a2, solver);

		space_i(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2, h_x, t_3, Mode);
		double* EitTP_0 = new double[N_x + 2]; double* EitTP_1 = new double[N_x + 2]; double* EitTP_2 = new double[N_x + 2];
		get_time_prime(EitTP_0, EitTP_1, EitTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Eit_vec(Mrows); Eigen::VectorXd EiN_t(Mrows);
		IMEX3_TS(Eit_vec, Eit_0, Eit_1, Eit_2, EiN_t, M, mrows, Mrows, EitTP_0, EitTP_1, EitTP_2, Ei_vec, Eif_vec, Eis_vec, EiN_s, tau, a0, b1, b2, a2, solver);

		space_r(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhos_0, rhos_1, rhos_2, ms_0, ms_1, ms_2, Ees_0, Ees_1, Ees_2, Eis_0, Eis_1, Eis_2, Ers_0, Ers_1, Ers_2, h_x, t_3, Mode);
		double* ErtTP_0 = new double[N_x + 2]; double* ErtTP_1 = new double[N_x + 2]; double* ErtTP_2 = new double[N_x + 2];
		get_time_prime(ErtTP_0, ErtTP_1, ErtTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		Eigen::VectorXd Ert_vec(Mrows); Eigen::VectorXd ErN_t(Mrows);
		IMEX3_TS(Ert_vec, Ert_0, Ert_1, Ert_2, ErN_t, M, mrows, Mrows, ErtTP_0, ErtTP_1, ErtTP_2, Er_vec, Erf_vec, Ers_vec, ErN_s, tau, a0, b1, b2, a2, solver);

		// update boundary
		update_boundary(rhot_0, rhot_1, rhot_2); update_boundary(mt_0, mt_1, mt_2);
		update_boundary(Eet_0, Eet_1, Eet_2); update_boundary(Eit_0, Eit_1, Eit_2); update_boundary(Ert_0, Ert_1, Ert_2);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(Sigma_0, Sigma_1, Sigma_2, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2);
			OE_compute_beta_j(Beta_j, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2);

			OE_procedure(rhot_1, rhot_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(mt_1, mt_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Eet_1, Eet_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Eit_1, Eit_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Ert_1, Ert_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);

			// update boundary
			update_boundary(rhot_0, rhot_1, rhot_2); update_boundary(mt_0, mt_1, mt_2);
			update_boundary(Eet_0, Eet_1, Eet_2); update_boundary(Eit_0, Eit_1, Eit_2); update_boundary(Ert_0, Ert_1, Ert_2);
		}
		

		// Last step of IMEX3
		Mode_select(rhot_0, rhot_1, rhot_2, Mode, N_x + 2); Mode_select(mt_0, mt_1, mt_2, Mode, N_x + 2); Mode_select(Eet_0, Eet_1, Eet_2, Mode, N_x + 2);
		Mode_select(Eit_0, Eit_1, Eit_2, Mode, N_x + 2); Mode_select(Ert_0, Ert_1, Ert_2, Mode, N_x + 2);

		double* rhol_0 = new double[N_x + 2]; double* rhol_1 = new double[N_x + 2]; double* rhol_2 = new double[N_x + 2];
		double* ml_0 = new double[N_x + 2]; double* ml_1 = new double[N_x + 2]; double* ml_2 = new double[N_x + 2];
		double* Eel_0 = new double[N_x + 2]; double* Eel_1 = new double[N_x + 2]; double* Eel_2 = new double[N_x + 2];
		double* Eil_0 = new double[N_x + 2]; double* Eil_1 = new double[N_x + 2]; double* Eil_2 = new double[N_x + 2];
		double* Erl_0 = new double[N_x + 2]; double* Erl_1 = new double[N_x + 2]; double* Erl_2 = new double[N_x + 2];

		compute_alphaC(alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2);
		compute_coef(coef_L, coef_R, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, h_x, Mode);

		space_d(S_0, S_1, S_2, alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, h_x);
		double* rholTP_0 = new double[N_x + 2]; double* rholTP_1 = new double[N_x + 2]; double* rholTP_2 = new double[N_x + 2];
		get_time_prime(rholTP_0, rholTP_1, rholTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			rhol_0[i] = rho_0[i] + tau * b1 * rhosTP_0[i] + tau * b2 * rhotTP_0[i] + tau * gamma * rholTP_0[i];
			rhol_1[i] = rho_1[i] + tau * b1 * rhosTP_1[i] + tau * b2 * rhotTP_1[i] + tau * gamma * rholTP_1[i];
			rhol_2[i] = rho_2[i] + tau * b1 * rhosTP_2[i] + tau * b2 * rhotTP_2[i] + tau * gamma * rholTP_2[i];
		}

		space_m(S_0, S_1, S_2, alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2, h_x);
		double* mlTP_0 = new double[N_x + 2]; double* mlTP_1 = new double[N_x + 2]; double* mlTP_2 = new double[N_x + 2];
		get_time_prime(mlTP_0, mlTP_1, mlTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			ml_0[i] = m_0[i] + tau * b1 * msTP_0[i] + tau * b2 * mtTP_0[i] + tau * gamma * mlTP_0[i];
			ml_1[i] = m_1[i] + tau * b1 * msTP_1[i] + tau * b2 * mtTP_1[i] + tau * gamma * mlTP_1[i];
			ml_2[i] = m_2[i] + tau * b1 * msTP_2[i] + tau * b2 * mtTP_2[i] + tau * gamma * mlTP_2[i];
		}
		
		space_e(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2, h_x, t_4, Mode);
		double* EelTP_0 = new double[N_x + 2]; double* EelTP_1 = new double[N_x + 2]; double* EelTP_2 = new double[N_x + 2];
		get_time_prime(EelTP_0, EelTP_1, EelTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			Eel_0[i] = Ee_0[i] + tau * b1 * EesTP_0[i] + tau * b2 * EetTP_0[i] + tau * gamma * EelTP_0[i];
			Eel_1[i] = Ee_1[i] + tau * b1 * EesTP_1[i] + tau * b2 * EetTP_1[i] + tau * gamma * EelTP_1[i];
			Eel_2[i] = Ee_2[i] + tau * b1 * EesTP_2[i] + tau * b2 * EetTP_2[i] + tau * gamma * EelTP_2[i];
		}

		space_i(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2, h_x, t_4, Mode);
		double* EilTP_0 = new double[N_x + 2]; double* EilTP_1 = new double[N_x + 2]; double* EilTP_2 = new double[N_x + 2];
		get_time_prime(EilTP_0, EilTP_1, EilTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			Eil_0[i] = Ei_0[i] + tau * b1 * EisTP_0[i] + tau * b2 * EitTP_0[i] + tau * gamma * EilTP_0[i];
			Eil_1[i] = Ei_1[i] + tau * b1 * EisTP_1[i] + tau * b2 * EitTP_1[i] + tau * gamma * EilTP_1[i];
			Eil_2[i] = Ei_2[i] + tau * b1 * EisTP_2[i] + tau * b2 * EitTP_2[i] + tau * gamma * EilTP_2[i];
		}

		space_r(S_0, S_1, S_2, coef_L, coef_R, X_1, X_2, X_3, X_4, X_5, alpha_c, rhot_0, rhot_1, rhot_2, mt_0, mt_1, mt_2, Eet_0, Eet_1, Eet_2, Eit_0, Eit_1, Eit_2, Ert_0, Ert_1, Ert_2, h_x, t_4, Mode);
		double* ErlTP_0 = new double[N_x + 2]; double* ErlTP_1 = new double[N_x + 2]; double* ErlTP_2 = new double[N_x + 2];
		get_time_prime(ErlTP_0, ErlTP_1, ErlTP_2, S_0, S_1, S_2, N_x + 2, h_x);
		for (int i = 1; i < N_x + 1; i++)
		{
			Erl_0[i] = Er_0[i] + tau * b1 * ErsTP_0[i] + tau * b2 * ErtTP_0[i] + tau * gamma * ErlTP_0[i];
			Erl_1[i] = Er_1[i] + tau * b1 * ErsTP_1[i] + tau * b2 * ErtTP_1[i] + tau * gamma * ErlTP_1[i];
			Erl_2[i] = Er_2[i] + tau * b1 * ErsTP_2[i] + tau * b2 * ErtTP_2[i] + tau * gamma * ErlTP_2[i];
		}

		// update boundary
		update_boundary(rhol_0, rhol_1, rhol_2); update_boundary(ml_0, ml_1, ml_2);
		update_boundary(Eel_0, Eel_1, Eel_2); update_boundary(Eil_0, Eil_1, Eil_2); update_boundary(Erl_0, Erl_1, Erl_2);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(Sigma_0, Sigma_1, Sigma_2, rhol_0, rhol_1, rhol_2, ml_0, ml_1, ml_2, Eel_0, Eel_1, Eel_2, Eil_0, Eil_1, Eil_2, Erl_0, Erl_1, Erl_2);
			OE_compute_beta_j(Beta_j, rhol_0, rhol_1, rhol_2, ml_0, ml_1, ml_2, Eel_0, Eel_1, Eel_2, Eil_0, Eil_1, Eil_2, Erl_0, Erl_1, Erl_2);

			OE_procedure(rhol_1, rhol_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(ml_1, ml_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Eel_1, Eel_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Eil_1, Eil_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);
			OE_procedure(Erl_1, Erl_2, Beta_j, Sigma_0, Sigma_1, Sigma_2, tau, h_x);

			// update boundary
			update_boundary(rhol_0, rhol_1, rhol_2); update_boundary(ml_0, ml_1, ml_2);
			update_boundary(Eel_0, Eel_1, Eel_2); update_boundary(Eil_0, Eil_1, Eil_2); update_boundary(Erl_0, Erl_1, Erl_2);
		}

		
		
		delete[] Sigma_0; delete[] Sigma_1; delete[] Sigma_2; delete[] Beta_j;  delete[] alpha_c;
		delete[] coef_L; delete[] coef_R;
		delete[] S_0; delete[] S_1; delete[] S_2;

		delete[] rhofTP_0; delete[] rhofTP_1; delete[] rhofTP_2; delete[] rhosTP_0; delete[] rhosTP_1; delete[] rhosTP_2;
		delete[] rhotTP_0; delete[] rhotTP_1; delete[] rhotTP_2; delete[] rholTP_0; delete[] rholTP_1; delete[] rholTP_2;

		delete[] mfTP_0; delete[] mfTP_1; delete[] mfTP_2; delete[] msTP_0; delete[] msTP_1; delete[] msTP_2;
		delete[] mtTP_0; delete[] mtTP_1; delete[] mtTP_2; delete[] mlTP_0; delete[] mlTP_1; delete[] mlTP_2;

		delete[] EefTP_0; delete[] EefTP_1; delete[] EefTP_2; delete[] EesTP_0; delete[] EesTP_1; delete[] EesTP_2;
		delete[] EetTP_0; delete[] EetTP_1; delete[] EetTP_2; delete[] EelTP_0; delete[] EelTP_1; delete[] EelTP_2;

		delete[] EifTP_0; delete[] EifTP_1; delete[] EifTP_2; delete[] EisTP_0; delete[] EisTP_1; delete[] EisTP_2;
		delete[] EitTP_0; delete[] EitTP_1; delete[] EitTP_2; delete[] EilTP_0; delete[] EilTP_1; delete[] EilTP_2;

		delete[] ErfTP_0; delete[] ErfTP_1; delete[] ErfTP_2; delete[] ErsTP_0; delete[] ErsTP_1; delete[] ErsTP_2;
		delete[] ErtTP_0; delete[] ErtTP_1; delete[] ErtTP_2; delete[] ErlTP_0; delete[] ErlTP_1; delete[] ErlTP_2;

		double rho_whole = 0.0; double rhonew_whole = 0.0;
		for (int i = 1; i < N_x + 1; i++)
		{
			rho_whole = rho_whole + rho_0[i];
		    rhonew_whole = rhonew_whole + rhol_0[i];
		}
		if (abs(rho_whole - rhonew_whole) > rho_c)
		{
			rho_c = abs(rho_whole - rhonew_whole);
		}
		//cout << abs(rho_whole - rhonew_whole) << endl;

		double m_whole = 0.0; double mnew_whole = 0.0;
		for (int i = 1; i < N_x + 1; i++)
		{
			m_whole = m_whole + m_0[i];
			mnew_whole = mnew_whole + ml_0[i];
		}
		if (abs(m_whole - mnew_whole) > m_c)
		{
			m_c = abs(m_whole - mnew_whole);
		}
		//cout << abs(m_whole - mnew_whole) << endl;

		for (int i = 0; i < N_x + 2; i++)
		{
			rho_0[i] = rhol_0[i]; rho_1[i] = rhol_1[i]; rho_2[i] = rhol_2[i];
			m_0[i] = ml_0[i]; m_1[i] = ml_1[i]; m_2[i] = ml_2[i];
			Ee_0[i] = Eel_0[i]; Ee_1[i] = Eel_1[i]; Ee_2[i] = Eel_2[i];
			Ei_0[i] = Eil_0[i]; Ei_1[i] = Eil_1[i]; Ei_2[i] = Eil_2[i];
			Er_0[i] = Erl_0[i]; Er_1[i] = Erl_1[i]; Er_2[i] = Erl_2[i];
			E_new[i] = Ee_0[i] + Ei_0[i] + Er_0[i];
		}
		
		double E_whole = 0.0; double Enew_whole = 0.0;
		for (int i = 1; i < N_x + 1; i++)
		{
			E_whole = E_whole + E[i];
			Enew_whole = Enew_whole + E_new[i];
		}
		if (abs(E_whole - Enew_whole) > E_c)
		{
			E_c = abs(E_whole - Enew_whole);
		}
		//cout << abs(E_whole - Enew_whole) << endl;


		t = t + tau;

		delete[] rhof_0; delete[] rhof_1; delete[] rhof_2; delete[] mf_0; delete[] mf_1; delete[] mf_2;
		delete[] Eef_0; delete[] Eef_1; delete[] Eef_2; delete[] Eif_0; delete[] Eif_1; delete[] Eif_2;
		delete[] Erf_0; delete[] Erf_1; delete[] Erf_2;

		delete[] rhos_0; delete[] rhos_1; delete[] rhos_2; delete[] ms_0; delete[] ms_1; delete[] ms_2;
		delete[] Ees_0; delete[] Ees_1; delete[] Ees_2; delete[] Eis_0; delete[] Eis_1; delete[] Eis_2;
		delete[] Ers_0; delete[] Ers_1; delete[] Ers_2;

		delete[] rhot_0; delete[] rhot_1; delete[] rhot_2; delete[] mt_0; delete[] mt_1; delete[] mt_2;
		delete[] Eet_0; delete[] Eet_1; delete[] Eet_2; delete[] Eit_0; delete[] Eit_1; delete[] Eit_2;
		delete[] Ert_0; delete[] Ert_1; delete[] Ert_2;

		delete[] rhol_0; delete[] rhol_1; delete[] rhol_2; delete[] ml_0; delete[] ml_1; delete[] ml_2;
		delete[] Eel_0; delete[] Eel_1; delete[] Eel_2; delete[] Eil_0; delete[] Eil_1; delete[] Eil_2;
		delete[] Erl_0; delete[] Erl_1; delete[] Erl_2;
	}

	cout << Tau << endl;

	delete[] X_1; delete[] X_2; delete[] X_3; delete[] X_4; delete[] X_5; delete[] E; delete[] E_new;
	delete[] s_j; delete[] beta_j;

	if (t == t_fin)
	{
		cout << "The computing of numerical solution is now finished ! \n" << endl;
	}

	// whether conservative
	cout << "the case of rho: " << rho_c << endl;
	cout << "the case of m: " << m_c << endl;
	cout << "the case of E: " << E_c << endl;

	// compute error
	cout << "Start computing the error ..." << endl;
	//error of rho
	double* rho_ev = rho_err(rho_0, rho_1, rho_2, x, h_x, t);
	cout << "L1 error of rho under mesh " << N_x << " = " << rho_ev[0] << endl;
	cout << "L2 error of rho under mesh " << N_x << " = " << rho_ev[1] << endl;
	cout << "Li error of rho under mesh " << N_x << " = " << rho_ev[2] << "\n" << endl;
	delete[] rho_0; delete[] rho_1; delete[] rho_2; delete[] rho_ev;

	//error of m
	double* m_ev = m_err(m_0, m_1, m_2, x, h_x, t);
	cout << "L1 error of m under mesh " << N_x << " = " << m_ev[0] << endl;
	cout << "L2 error of m under mesh " << N_x << " = " << m_ev[1] << endl;
	cout << "Li error of m under mesh " << N_x << " = " << m_ev[2] << "\n" << endl;
	delete[] m_0; delete[] m_1; delete[] m_2; delete[] m_ev;

	//error of Ee
	double* Ee_ev = Ee_err(Ee_0, Ee_1, Ee_2, x, h_x, t);
	cout << "L1 error of Ee under mesh " << N_x << " = " << Ee_ev[0] << endl;
	cout << "L2 error of Ee under mesh " << N_x << " = " << Ee_ev[1] << endl;
	cout << "Li error of Ee under mesh " << N_x << " = " << Ee_ev[2] << "\n" << endl;
	delete[] Ee_0; delete[] Ee_1; delete[] Ee_2; delete[] Ee_ev;

	//error of Ei
	double* Ei_ev = Ei_err(Ei_0, Ei_1, Ei_2, x, h_x, t);
	cout << "L1 error of Ei under mesh " << N_x << " = " << Ei_ev[0] << endl;
	cout << "L2 error of Ei under mesh " << N_x << " = " << Ei_ev[1] << endl;
	cout << "Li error of Ei under mesh " << N_x << " = " << Ei_ev[2] << "\n" << endl;
	delete[] Ei_0; delete[] Ei_1; delete[] Ei_2; delete[] Ei_ev;

	//error of Er
	double* Er_ev = Er_err(Er_0, Er_1, Er_2, x, h_x, t);
	cout << "L1 error of Er under mesh " << N_x << " = " << Er_ev[0] << endl;
	cout << "L2 error of Er under mesh " << N_x << " = " << Er_ev[1] << endl;
	cout << "Li error of Er under mesh " << N_x << " = " << Er_ev[2] << "\n" << endl;
	delete[] Er_0; delete[] Er_1; delete[] Er_2; delete[] Er_ev;

	delete[] x;
	cout << "The computing of error is now finished ! \n" << endl;

	// timing ended
	end = std::time(NULL);
	cout << "Time has passed " << end - start << " seconds." << endl;
	return 0;
}
