#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG.h"
#include "time.h"

void piece_Linear_coef_matrix(Eigen::SparseMatrix<double>& mat, int M_rows, int M_cols, int rows, int cols, double alpha_0, double beta_0, double h, double m_index, double coff)
{
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(M_cols);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (i == j + 1)
			{
				//M1  M5  M9
				//主对角线下方对角线
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h * h)));
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * (6.0 - 3.0 * alpha_0) / (h * h)));
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * (5.0 * (alpha_0 - 3.0) - 15.0 + 60.0 * beta_0) / (h * h)));

				// M2  M6
				//主对角线下方对角线
				tripletList.push_back(T(i, j + cols, 0.0 + coff * (alpha_0 - 1.0) / (h * h)));
				tripletList.push_back(T(i + rows, j + 2 * cols, 0.0 + coff * (9.0 - 3.0 * (alpha_0 - 1.0) - 36.0 * beta_0) / (h * h)));

				// M4  M8
				//主对角线下方对角线
				tripletList.push_back(T(i + rows, j, 0.0 + coff * 3.0 * (1.0 - alpha_0) / (h * h)));
				tripletList.push_back(T(i + 2 * rows, j + cols, 0.0 + coff * (5.0 * alpha_0 - 20.0) / (h * h)));

				// M3
				//主对角线下方对角线
				tripletList.push_back(T(i, j + 2 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / (h * h)));

				// M7
				//主对角线下方对角线
				tripletList.push_back(T(i + 2 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h * h)));
			}
			if (i == j)
			{
				//M1  M5  M9
				//主对角线
				tripletList.push_back(T(i, j, m_index + coff * -2.0 * alpha_0 / (h * h)));
				tripletList.push_back(T(i + rows, j + cols, m_index + coff * -6.0 * alpha_0 / (h * h)));
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, m_index + coff * (-30.0 - 10.0 * (alpha_0 - 3.0) - 120.0 * beta_0) / (h * h)));

				// M3
				//主对角线
				tripletList.push_back(T(i, j + 2 * cols, 0.0 + coff * (6.0 - 2.0 * alpha_0 - 24.0 * beta_0) / (h * h)));

				// M7
				//主对角线
				tripletList.push_back(T(i + 2 * rows, j, 0.0 + coff * 10.0 * (3.0 - alpha_0) / (h * h)));
			}
			if (i == j - 1)
			{
				//M1  M5  M9
				//主对角线上方对角线
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h * h)));
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * (6.0 - 3.0 * alpha_0) / (h * h)));
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * (5.0 * (alpha_0 - 3.0) - 15.0 + 60.0 * beta_0) / (h * h)));

				// M2  M6
				//主对角线上方对角线
				tripletList.push_back(T(i, j + cols, 0.0 + coff * (1.0 - alpha_0) / (h * h)));
				tripletList.push_back(T(i + rows, j + 2 * cols, 0.0 + coff * (3.0 * (alpha_0 - 1.0) - 9.0 + 36.0 * beta_0) / (h * h)));

				// M4  M8
				//主对角线上方对角线
				tripletList.push_back(T(i + rows, j, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / (h * h)));
				tripletList.push_back(T(i + 2 * rows, j + cols, 0.0 + coff * (20.0 - 5.0 * alpha_0) / (h * h)));

				// M3
				//主对角线上方对角线
				tripletList.push_back(T(i, j + 2 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / (h * h)));

				// M7
				//主对角线上方对角线
				tripletList.push_back(T(i + 2 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h * h)));
			}
		}
	}
	//Mass matrix M1
	tripletList.push_back(T(0, cols - 1, 0.0 + coff * alpha_0 / (h * h))); tripletList.push_back(T(rows - 1, 0, 0.0 + coff * alpha_0 / (h * h)));
	//Mass matrix M5
	tripletList.push_back(T(rows, 2 * cols - 1, 0.0 + coff * (6.0 - 3.0 * alpha_0) / (h * h))); tripletList.push_back(T(2 * rows - 1, cols, 0.0 + coff * (6.0 - 3.0 * alpha_0) / (h * h)));
	//Mass matrix M9
	tripletList.push_back(T(2 * rows, 3 * cols - 1, 0.0 + coff * (5.0 * (alpha_0 - 3.0) - 15.0 + 60.0 * beta_0) / (h * h))); tripletList.push_back(T(3 * rows - 1, 2 * cols, 0.0 + coff * (5.0 * (alpha_0 - 3.0) - 15.0 + 60.0 * beta_0) / (h * h)));


	// Mass matrix M2
	tripletList.push_back(T(0, 2 * cols - 1, 0.0 + coff * (alpha_0 - 1.0) / (h * h))); tripletList.push_back(T(rows - 1, cols, 0.0 + coff * (1.0 - alpha_0) / (h * h)));
	// Mass matrix M6
	tripletList.push_back(T(rows, 3 * cols - 1, 0.0 + coff * (9.0 - 3.0 * (alpha_0 - 1.0) - 36.0 * beta_0) / (h * h))); tripletList.push_back(T(2 * rows - 1, 2 * cols, 0.0 + coff * (3.0 * (alpha_0 - 1.0) - 9.0 + 36.0 * beta_0) / (h * h)));


	// Mass matrix M4
	tripletList.push_back(T(rows, cols - 1, 0.0 + coff * 3.0 * (1.0 - alpha_0) / (h * h))); tripletList.push_back(T(2 * rows - 1, 0, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / (h * h)));
	// Mass matrix M8
	tripletList.push_back(T(2 * rows, 2 * cols - 1, 0.0 + coff * (5.0 * alpha_0 - 20.0) / (h * h))); tripletList.push_back(T(3 * rows - 1, cols, 0.0 + coff * (20.0 - 5.0 * alpha_0) / (h * h)));


	// Mass matrix M3
	tripletList.push_back(T(0, 3 * cols - 1, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / (h * h))); tripletList.push_back(T(rows - 1, 2 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / (h * h)));

	// Mass matrix M7
	tripletList.push_back(T(2 * rows, cols - 1, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h * h))); tripletList.push_back(T(3 * rows - 1, 0, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h * h)));


	mat.setFromTriplets(tripletList.begin(), tripletList.end());

}


void IMEX3_FS(Eigen::VectorXd& Uf_vec, double* Uf_0, double* Uf_1, double* Uf_2, Eigen::VectorXd& U_vec, Eigen::VectorXd& N_f, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, double* U_0, double* U_1, double* U_2, double tau, double a0, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(Mrows); Eigen::VectorXd UR1(Mrows); Eigen::VectorXd Uf_N(Mrows);

	for (int j = 1; j < mrows + 1; j++)
	{
		UTP_vec(j - 1) = UTP_0[j];
		UTP_vec(j - 1 + mrows) = UTP_1[j];
		UTP_vec(j - 1 + 2 * mrows) = UTP_2[j];

		U_vec(j - 1) = U_0[j];
		U_vec(j - 1 + mrows) = U_1[j];
		U_vec(j - 1 + 2 * mrows) = U_2[j];
	}
	UR1 = M * U_vec;

	for (int i = 0; i < Mrows; i++)
	{
		N_f(i) = UTP_vec(i) - a0 * UR1(i);
		Uf_N(i) = U_vec(i) + tau * gamma * N_f(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < Mrows; i++)
		{
			Uf_vec(i) = Uf_N(i);
		}
	}
	else
	{
		Uf_vec = solver.solve(Uf_N);
	}

	for (int j = 0; j < mrows; j++)
	{
		Uf_0[j + 1] = Uf_vec(j);
		Uf_1[j + 1] = Uf_vec(j + mrows);
		Uf_2[j + 1] = Uf_vec(j + 2 * mrows);
	}
}

void IMEX3_SS(Eigen::VectorXd& Us_vec, double* Us_0, double* Us_1, double* Us_2, Eigen::VectorXd& N_s, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd N_f, double tau, double a0, double a1, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(Mrows); Eigen::VectorXd UR2(Mrows); Eigen::VectorXd Us_N(Mrows);

	for (int j = 1; j < mrows + 1; j++)
	{
		UTP_vec(j - 1) = UTP_0[j];
		UTP_vec(j - 1 + mrows) = UTP_1[j];
		UTP_vec(j - 1 + 2 * mrows) = UTP_2[j];
	}

	UR2 = M * Uf_vec;

	for (int i = 0; i < Mrows; i++)
	{
		N_s(i) = UTP_vec(i) - a0 * UR2(i);
		Us_N(i) = U_vec(i) + tau * (1.0 - gamma) / 2.0 * a0 * UR2(i) + tau * ((1.0 + gamma) / 2.0 - a1) * N_f(i) + tau * a1 * N_s(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < Mrows; i++)
		{
			Us_vec(i) = Us_N(i);
		}
	}
	else
	{
		Us_vec = solver.solve(Us_N);
	}

	for (int j = 0; j < mrows; j++)
	{
		Us_0[j + 1] = Us_vec(j);
		Us_1[j + 1] = Us_vec(j + mrows);
		Us_2[j + 1] = Us_vec(j + 2 * mrows);
	}
}

void IMEX3_TS(Eigen::VectorXd& Ut_vec, double* Ut_0, double* Ut_1, double* Ut_2, Eigen::VectorXd& N_t, Eigen::SparseMatrix<double> M, int mrows, int Mrows, double* UTP_0, double* UTP_1, double* UTP_2, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd Us_vec, Eigen::VectorXd N_s, double tau, double a0, double b1, double b2, double a2, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(Mrows); Eigen::VectorXd UR2(Mrows); Eigen::VectorXd UR3(Mrows); Eigen::VectorXd Ut_N(Mrows);

	for (int j = 1; j < mrows + 1; j++)
	{
		UTP_vec(j - 1) = UTP_0[j];
		UTP_vec(j - 1 + mrows) = UTP_1[j];
		UTP_vec(j - 1 + 2 * mrows) = UTP_2[j];
	}

	UR2 = M * Uf_vec; UR3 = M * Us_vec;

	for (int i = 0; i < Mrows; i++)
	{
		N_t(i) = UTP_vec(i) - a0 * UR3(i);
		Ut_N(i) = U_vec(i) + tau * b1 * a0 * UR2(i) + tau * b2 * a0 * UR3(i) + tau * (1.0 - a2) * N_s(i) + tau * a2 * N_t(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < Mrows; i++)
		{
			Ut_vec(i) = Ut_N(i);
		}
	}
	else
	{
		Ut_vec = solver.solve(Ut_N);
	}

	for (int j = 0; j < mrows; j++)
	{
		Ut_0[j + 1] = Ut_vec(j);
		Ut_1[j + 1] = Ut_vec(j + mrows);
		Ut_2[j + 1] = Ut_vec(j + 2 * mrows);
	}
}

void compute_st(double* s_j, double* rho_0, double* m_0, double* Ee_0)
{
	double temp_1 = a * omega_ei * pow(C_ve, 4) + a * omega_ei * pow(C_ve, 3) * C_vi;
	double temp_2 = 4.0 * a * omega_ei * omega_er * pow(C_ve, 4) * C_vi;


	double* e_e = new double[N_x]; double* s1 = new double[N_x]; double* s2 = new double[N_x]; double* s = new double[N_x];
	double* alp_1 = new double[N_x]; double* alp_2 = new double[N_x];

	for (int i = 0; i < N_x; i++)
	{
		double temp_rhou2 = pow(m_0[i + 1], 2) / rho_0[i + 1];
		e_e[i] = (Ee_0[i + 1] - temp_rhou2 / 6.0) / rho_0[i + 1];
	}

	for (int i = 0; i < N_x; i++)
	{
		double temp_s1_1 = temp_1 + 4.0 * a * omega_er * C_vi * pow(e_e[i], 3) + omega_er * pow(C_ve, 4) * C_vi * rho_0[i + 1];
		double temp_s1_2 = temp_2 * (4.0 * a * pow(e_e[i], 3) + pow(C_ve, 4) * rho_0[i + 1] + pow(C_ve, 3) * C_vi * rho_0[i + 1]);
		s1[i] = sqrt(pow(temp_s1_1, 2) - temp_s1_2);
		s2[i] = -1.0 * temp_s1_1;
		s[i] = 2.0 * a * pow(C_ve, 4) * C_vi * rho_0[i + 1];
		alp_1[i] = (s2[i] - s1[i]) / s[i];
		alp_2[i] = (s2[i] + s1[i]) / s[i];
	}

	for (int i = 0; i < N_x; i++)
	{
		double temp_max = abs(alp_1[i]);
		if (abs(alp_2[i]) > temp_max)
		{
			temp_max = abs(alp_2[i]);
		}
		s_j[i] = temp_max;
	}

	delete[] e_e; delete[] s1; delete[] s2; delete[] s; delete[] alp_1; delete[] alp_2;
}

double compute_sj_max(double* s_j)
{
	double sj_max = 0.0;
	for (int i = 0; i < N_x; i++)
	{
		if (s_j[i] > sj_max)
		{
			sj_max = s_j[i];
		}
	}
	return sj_max;
}

double comput_rho_min(double* rho)
{
	double temp_min = 10.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		if (rho[i] < temp_min)
		{
			temp_min = rho[i];
		}
	}
	return temp_min;
}

double compute_dj(double rho_min)
{
	double temp_max = 0.0;
	if (kappa_e / (C_ve * rho_min) > kappa_i / (C_vi * rho_min))
	{
		temp_max = kappa_e / (C_ve * rho_min);
	}
	else
	{
		temp_max = kappa_i / (C_vi * rho_min);
	}
	if (kappa_r / a > temp_max)
	{
		temp_max = kappa_r / a;
	}
	return temp_max;
}

