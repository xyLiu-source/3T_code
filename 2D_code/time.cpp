#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG2D.h"
#include "time.h"

// time marching
void IMEX3_f(double** uf_0, double** uf_1, double** uf_2, double** uf_3, double** uf_4, double** uf_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double** ufTP_0, double** ufTP_1, double** ufTP_2, double** ufTP_3, double** ufTP_4, double** ufTP_5, double tau)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			uf_0[i][j] = u_0[i][j] + tau * gamma * ufTP_0[i][j];
			uf_1[i][j] = u_1[i][j] + tau * gamma * ufTP_1[i][j];
			uf_2[i][j] = u_2[i][j] + tau * gamma * ufTP_2[i][j];
			uf_3[i][j] = u_3[i][j] + tau * gamma * ufTP_3[i][j];
			uf_4[i][j] = u_4[i][j] + tau * gamma * ufTP_4[i][j];
			uf_5[i][j] = u_5[i][j] + tau * gamma * ufTP_5[i][j];
		}
	}
}

void IMEX3_s(double** us_0, double** us_1, double** us_2, double** us_3, double** us_4, double** us_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double** uf_0, double** uf_1, double** uf_2, double** uf_3, double** uf_4, double** uf_5, double** ufTP_0, double** ufTP_1, double** ufTP_2, double** ufTP_3, double** ufTP_4, double** ufTP_5, double** usTP_0, double** usTP_1, double** usTP_2, double** usTP_3, double** usTP_4, double** usTP_5, double tau, double a1)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			us_0[i][j] = u_0[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_0[i][j] + tau * a1 * usTP_0[i][j];
			us_1[i][j] = u_1[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_1[i][j] + tau * a1 * usTP_1[i][j];
			us_2[i][j] = u_2[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_2[i][j] + tau * a1 * usTP_2[i][j];
			us_3[i][j] = u_3[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_3[i][j] + tau * a1 * usTP_3[i][j];
			us_4[i][j] = u_4[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_4[i][j] + tau * a1 * usTP_4[i][j];
			us_5[i][j] = u_5[i][j] + tau * ((1.0 + gamma) / 2.0 - a1) * ufTP_5[i][j] + tau * a1 * usTP_5[i][j];
		}
	}
}

void IMEX3_t(double** ut_0, double** ut_1, double** ut_2, double** ut_3, double** ut_4, double** ut_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double** uf_0, double** uf_1, double** uf_2, double** uf_3, double** uf_4, double** uf_5, double** us_0, double** us_1, double** us_2, double** us_3, double** us_4, double** us_5, double** usTP_0, double** usTP_1, double** usTP_2, double** usTP_3, double** usTP_4, double** usTP_5, double** utTP_0, double** utTP_1, double** utTP_2, double** utTP_3, double** utTP_4, double** utTP_5, double tau, double a2)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			ut_0[i][j] = u_0[i][j] + tau * (1.0 - a2) * usTP_0[i][j] + tau * a2 * utTP_0[i][j];
			ut_1[i][j] = u_1[i][j] + tau * (1.0 - a2) * usTP_1[i][j] + tau * a2 * utTP_1[i][j];
			ut_2[i][j] = u_2[i][j] + tau * (1.0 - a2) * usTP_2[i][j] + tau * a2 * utTP_2[i][j];
			ut_3[i][j] = u_3[i][j] + tau * (1.0 - a2) * usTP_3[i][j] + tau * a2 * utTP_3[i][j];
			ut_4[i][j] = u_4[i][j] + tau * (1.0 - a2) * usTP_4[i][j] + tau * a2 * utTP_4[i][j];
			ut_5[i][j] = u_5[i][j] + tau * (1.0 - a2) * usTP_5[i][j] + tau * a2 * utTP_5[i][j];
		}
	}
}

void IMEX3_l(double** ul_0, double** ul_1, double** ul_2, double** ul_3, double** ul_4, double** ul_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double** uf_0, double** uf_1, double** uf_2, double** uf_3, double** uf_4, double** uf_5, double** us_0, double** us_1, double** us_2, double** us_3, double** us_4, double** us_5, double** ut_0, double** ut_1, double** ut_2, double** ut_3, double** ut_4, double** ut_5, double** usTP_0, double** usTP_1, double** usTP_2, double** usTP_3, double** usTP_4, double** usTP_5, double** utTP_0, double** utTP_1, double** utTP_2, double** utTP_3, double** utTP_4, double** utTP_5, double** ulTP_0, double** ulTP_1, double** ulTP_2, double** ulTP_3, double** ulTP_4, double** ulTP_5, double tau, double b1, double b2)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			ul_0[i][j] = u_0[i][j] + tau * b1 * usTP_0[i][j] + tau * b2 * utTP_0[i][j] + tau * gamma * ulTP_0[i][j];
			ul_1[i][j] = u_1[i][j] + tau * b1 * usTP_1[i][j] + tau * b2 * utTP_1[i][j] + tau * gamma * ulTP_1[i][j];
			ul_2[i][j] = u_2[i][j] + tau * b1 * usTP_2[i][j] + tau * b2 * utTP_2[i][j] + tau * gamma * ulTP_2[i][j];
			ul_3[i][j] = u_3[i][j] + tau * b1 * usTP_3[i][j] + tau * b2 * utTP_3[i][j] + tau * gamma * ulTP_3[i][j];
			ul_4[i][j] = u_4[i][j] + tau * b1 * usTP_4[i][j] + tau * b2 * utTP_4[i][j] + tau * gamma * ulTP_4[i][j];
			ul_5[i][j] = u_5[i][j] + tau * b1 * usTP_5[i][j] + tau * b2 * utTP_5[i][j] + tau * gamma * ulTP_5[i][j];
		}
	}
}


void IMEX3_FS(Eigen::VectorXd& Uf_vec, double** Uf_0, double** Uf_1, double** Uf_2, double** Uf_3, double** Uf_4, double** Uf_5, Eigen::VectorXd& U_vec, Eigen::VectorXd& N_f, Eigen::SparseMatrix<double> M, int m_rows, int M_rows, double** UTP_0, double** UTP_1, double** UTP_2, double** UTP_3, double** UTP_4, double** UTP_5, double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double tau, double a0, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(M_rows); Eigen::VectorXd UR1(M_rows); Eigen::VectorXd Uf_N(M_rows);

	double* UTP0_vec = new double[m_rows]; double* UTP1_vec = new double[m_rows]; double* UTP2_vec = new double[m_rows];
	double* UTP3_vec = new double[m_rows]; double* UTP4_vec = new double[m_rows]; double* UTP5_vec = new double[m_rows];

	double* U0_vec = new double[m_rows]; double* U1_vec = new double[m_rows]; double* U2_vec = new double[m_rows];
	double* U3_vec = new double[m_rows]; double* U4_vec = new double[m_rows]; double* U5_vec = new double[m_rows];

	Matrix_to_vector(UTP0_vec, UTP_0, N_y + 2, N_x + 2); Matrix_to_vector(UTP1_vec, UTP_1, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP2_vec, UTP_2, N_y + 2, N_x + 2); Matrix_to_vector(UTP3_vec, UTP_3, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP4_vec, UTP_4, N_y + 2, N_x + 2); Matrix_to_vector(UTP5_vec, UTP_5, N_y + 2, N_x + 2);

	Matrix_to_vector(U0_vec, U_0, N_y + 2, N_x + 2); Matrix_to_vector(U1_vec, U_1, N_y + 2, N_x + 2);
	Matrix_to_vector(U2_vec, U_2, N_y + 2, N_x + 2); Matrix_to_vector(U3_vec, U_3, N_y + 2, N_x + 2);
	Matrix_to_vector(U4_vec, U_4, N_y + 2, N_x + 2); Matrix_to_vector(U5_vec, U_5, N_y + 2, N_x + 2);

	for (int j = 0; j < m_rows; j++)
	{
		UTP_vec(j) = UTP0_vec[j];
		UTP_vec(j + m_rows) = UTP1_vec[j];
		UTP_vec(j + 2 * m_rows) = UTP2_vec[j];
		UTP_vec(j + 3 * m_rows) = UTP3_vec[j];
		UTP_vec(j + 4 * m_rows) = UTP4_vec[j];
		UTP_vec(j + 5 * m_rows) = UTP5_vec[j];

		U_vec(j) = U0_vec[j];
		U_vec(j + m_rows) = U1_vec[j];
		U_vec(j + 2 * m_rows) = U2_vec[j];
		U_vec(j + 3 * m_rows) = U3_vec[j];
		U_vec(j + 4 * m_rows) = U4_vec[j];
		U_vec(j + 5 * m_rows) = U5_vec[j];
	}

	UR1 = M * U_vec;

	for (int i = 0; i < M_rows; i++)
	{
		N_f(i) = UTP_vec(i) - a0 * UR1(i);
		Uf_N(i) = U_vec(i) + tau * gamma * N_f(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < M_rows; i++)
		{
			Uf_vec(i) = Uf_N(i);
		}
	}
	else
	{
		Uf_vec = solver.solve(Uf_N);
	}

	double* Uf0_vec = new double[m_rows]; double* Uf1_vec = new double[m_rows]; double* Uf2_vec = new double[m_rows];
	double* Uf3_vec = new double[m_rows]; double* Uf4_vec = new double[m_rows]; double* Uf5_vec = new double[m_rows];

	for (int j = 0; j < m_rows; j++)
	{
		Uf0_vec[j] = Uf_vec(j);
		Uf1_vec[j] = Uf_vec(j + m_rows);
		Uf2_vec[j] = Uf_vec(j + 2 * m_rows);
		Uf3_vec[j] = Uf_vec(j + 3 * m_rows);
		Uf4_vec[j] = Uf_vec(j + 4 * m_rows);
		Uf5_vec[j] = Uf_vec(j + 5 * m_rows);
	}

	vector_to_Matrix(Uf_0, Uf0_vec, N_y + 2, N_x + 2); vector_to_Matrix(Uf_1, Uf1_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Uf_2, Uf2_vec, N_y + 2, N_x + 2); vector_to_Matrix(Uf_3, Uf3_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Uf_4, Uf4_vec, N_y + 2, N_x + 2); vector_to_Matrix(Uf_5, Uf5_vec, N_y + 2, N_x + 2);

	delete[] UTP0_vec; delete[] UTP1_vec; delete[] UTP2_vec; delete[] UTP3_vec; delete[] UTP4_vec; delete[] UTP5_vec;
	delete[] U0_vec; delete[] U1_vec; delete[] U2_vec; delete[] U3_vec; delete[] U4_vec; delete[] U5_vec;
	delete[] Uf0_vec; delete[] Uf1_vec; delete[] Uf2_vec; delete[] Uf3_vec; delete[] Uf4_vec; delete[] Uf5_vec;
}

void IMEX3_SS(Eigen::VectorXd& Us_vec, double** Us_0, double** Us_1, double** Us_2, double** Us_3, double** Us_4, double** Us_5, Eigen::VectorXd& N_s, Eigen::SparseMatrix<double> M, int m_rows, int M_rows, double** UTP_0, double** UTP_1, double** UTP_2, double** UTP_3, double** UTP_4, double** UTP_5, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd N_f, double tau, double a0, double a1, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(M_rows); Eigen::VectorXd UR2(M_rows); Eigen::VectorXd Us_N(M_rows);

	double* UTP0_vec = new double[m_rows]; double* UTP1_vec = new double[m_rows]; double* UTP2_vec = new double[m_rows];
	double* UTP3_vec = new double[m_rows]; double* UTP4_vec = new double[m_rows]; double* UTP5_vec = new double[m_rows];

	Matrix_to_vector(UTP0_vec, UTP_0, N_y + 2, N_x + 2); Matrix_to_vector(UTP1_vec, UTP_1, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP2_vec, UTP_2, N_y + 2, N_x + 2); Matrix_to_vector(UTP3_vec, UTP_3, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP4_vec, UTP_4, N_y + 2, N_x + 2); Matrix_to_vector(UTP5_vec, UTP_5, N_y + 2, N_x + 2);

	for (int j = 0; j < m_rows; j++)
	{
		UTP_vec(j) = UTP0_vec[j];
		UTP_vec(j + m_rows) = UTP1_vec[j];
		UTP_vec(j + 2 * m_rows) = UTP2_vec[j];
		UTP_vec(j + 3 * m_rows) = UTP3_vec[j];
		UTP_vec(j + 4 * m_rows) = UTP4_vec[j];
		UTP_vec(j + 5 * m_rows) = UTP5_vec[j];
	}

	UR2 = M * Uf_vec;

	for (int i = 0; i < M_rows; i++)
	{
		N_s(i) = UTP_vec(i) - a0 * UR2(i);
		Us_N(i) = U_vec(i) + tau * (1.0 - gamma) / 2.0 * a0 * UR2(i) + tau * ((1.0 + gamma) / 2.0 - a1) * N_f(i) + tau * a1 * N_s(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < M_rows; i++)
		{
			Us_vec(i) = Us_N(i);
		}
	}
	else
	{
		Us_vec = solver.solve(Us_N);
	}

	double* Us0_vec = new double[m_rows]; double* Us1_vec = new double[m_rows]; double* Us2_vec = new double[m_rows];
	double* Us3_vec = new double[m_rows]; double* Us4_vec = new double[m_rows]; double* Us5_vec = new double[m_rows];

	for (int j = 0; j < m_rows; j++)
	{
		Us0_vec[j] = Us_vec(j);
		Us1_vec[j] = Us_vec(j + m_rows);
		Us2_vec[j] = Us_vec(j + 2 * m_rows);
		Us3_vec[j] = Us_vec(j + 3 * m_rows);
		Us4_vec[j] = Us_vec(j + 4 * m_rows);
		Us5_vec[j] = Us_vec(j + 5 * m_rows);
	}

	vector_to_Matrix(Us_0, Us0_vec, N_y + 2, N_x + 2); vector_to_Matrix(Us_1, Us1_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Us_2, Us2_vec, N_y + 2, N_x + 2); vector_to_Matrix(Us_3, Us3_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Us_4, Us4_vec, N_y + 2, N_x + 2); vector_to_Matrix(Us_5, Us5_vec, N_y + 2, N_x + 2);

	delete[] UTP0_vec; delete[] UTP1_vec; delete[] UTP2_vec; delete[] UTP3_vec; delete[] UTP4_vec; delete[] UTP5_vec;
	delete[] Us0_vec; delete[] Us1_vec; delete[] Us2_vec; delete[] Us3_vec; delete[] Us4_vec; delete[] Us5_vec;
}

void IMEX3_TS(Eigen::VectorXd& Ut_vec, double** Ut_0, double** Ut_1, double** Ut_2, double** Ut_3, double** Ut_4, double** Ut_5, Eigen::VectorXd& N_t, Eigen::SparseMatrix<double> M, int m_rows, int M_rows, double** UTP_0, double** UTP_1, double** UTP_2, double** UTP_3, double** UTP_4, double** UTP_5, Eigen::VectorXd U_vec, Eigen::VectorXd Uf_vec, Eigen::VectorXd Us_vec, Eigen::VectorXd N_s, double tau, double a0, double b1, double b2, double a2, Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
{
	Eigen::VectorXd UTP_vec(M_rows); Eigen::VectorXd UR2(M_rows); Eigen::VectorXd UR3(M_rows); Eigen::VectorXd Ut_N(M_rows);

	double* UTP0_vec = new double[m_rows]; double* UTP1_vec = new double[m_rows]; double* UTP2_vec = new double[m_rows];
	double* UTP3_vec = new double[m_rows]; double* UTP4_vec = new double[m_rows]; double* UTP5_vec = new double[m_rows];

	Matrix_to_vector(UTP0_vec, UTP_0, N_y + 2, N_x + 2); Matrix_to_vector(UTP1_vec, UTP_1, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP2_vec, UTP_2, N_y + 2, N_x + 2); Matrix_to_vector(UTP3_vec, UTP_3, N_y + 2, N_x + 2);
	Matrix_to_vector(UTP4_vec, UTP_4, N_y + 2, N_x + 2); Matrix_to_vector(UTP5_vec, UTP_5, N_y + 2, N_x + 2);

	for (int j = 0; j < m_rows; j++)
	{
		UTP_vec(j) = UTP0_vec[j];
		UTP_vec(j + m_rows) = UTP1_vec[j];
		UTP_vec(j + 2 * m_rows) = UTP2_vec[j];
		UTP_vec(j + 3 * m_rows) = UTP3_vec[j];
		UTP_vec(j + 4 * m_rows) = UTP4_vec[j];
		UTP_vec(j + 5 * m_rows) = UTP5_vec[j];
	}

	UR2 = M * Uf_vec; UR3 = M * Us_vec;

	for (int i = 0; i < M_rows; i++)
	{
		N_t(i) = UTP_vec(i) - a0 * UR3(i);
		Ut_N(i) = U_vec(i) + tau * b1 * a0 * UR2(i) + tau * b2 * a0 * UR3(i) + tau * (1.0 - a2) * N_s(i) + tau * a2 * N_t(i);
	}

	if (a0 == 0.0)
	{
		for (int i = 0; i < M_rows; i++)
		{
			Ut_vec(i) = Ut_N(i);
		}
	}
	else
	{
		Ut_vec = solver.solve(Ut_N);
	}

	double* Ut0_vec = new double[m_rows]; double* Ut1_vec = new double[m_rows]; double* Ut2_vec = new double[m_rows];
	double* Ut3_vec = new double[m_rows]; double* Ut4_vec = new double[m_rows]; double* Ut5_vec = new double[m_rows];

	for (int j = 0; j < m_rows; j++)
	{
		Ut0_vec[j] = Ut_vec(j);
		Ut1_vec[j] = Ut_vec(j + m_rows);
		Ut2_vec[j] = Ut_vec(j + 2 * m_rows);
		Ut3_vec[j] = Ut_vec(j + 3 * m_rows);
		Ut4_vec[j] = Ut_vec(j + 4 * m_rows);
		Ut5_vec[j] = Ut_vec(j + 5 * m_rows);
	}

	vector_to_Matrix(Ut_0, Ut0_vec, N_y + 2, N_x + 2); vector_to_Matrix(Ut_1, Ut1_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Ut_2, Ut2_vec, N_y + 2, N_x + 2); vector_to_Matrix(Ut_3, Ut3_vec, N_y + 2, N_x + 2);
	vector_to_Matrix(Ut_4, Ut4_vec, N_y + 2, N_x + 2); vector_to_Matrix(Ut_5, Ut5_vec, N_y + 2, N_x + 2);

	delete[] UTP0_vec; delete[] UTP1_vec; delete[] UTP2_vec; delete[] UTP3_vec; delete[] UTP4_vec; delete[] UTP5_vec;
	delete[] Ut0_vec; delete[] Ut1_vec; delete[] Ut2_vec; delete[] Ut3_vec; delete[] Ut4_vec; delete[] Ut5_vec;
}

void piece_Linear_coef_matrix(Eigen::SparseMatrix<double>& mat, int M_rows, int M_cols, int rows, int cols, double alpha_0, double beta_0, double h_x, double h_y, double m_index, double coff)
{
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(M_cols);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j == i - N_x) // 主对角线元素下方元素（关于y） 
			{
				// M_00
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h_y * h_y)));
				// M_02
				tripletList.push_back(T(i, j + 2 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_y, 2)));
				// M_05 
				tripletList.push_back(T(i, j + 5 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_y, 2)));

				// M_11
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * alpha_0 / pow(h_y, 2)));
				// M_13
				tripletList.push_back(T(i + rows, j + 3 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_y, 2)));

				// M_20
				tripletList.push_back(T(i + 2 * rows, j, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_y, 2)));
				// M_22
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));
				// M_25
				tripletList.push_back(T(i + 2 * rows, j + 5 * cols, 0.0 + coff * (9.0 + 3.0 * (1.0 - alpha_0) - 36.0 * beta_0) / pow(h_y, 2)));

				// M_31
				tripletList.push_back(T(i + 3 * rows, j + cols, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_y, 2)));
				// M_33
				tripletList.push_back(T(i + 3 * rows, j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));

				// M_44
				tripletList.push_back(T(i + 4 * rows, j + 4 * cols, 0.0 + coff * alpha_0 / (h_y * h_y)));

				// M_50
				tripletList.push_back(T(i + 5 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_y * h_y)));
				// M_52
				tripletList.push_back(T(i + 5 * rows, j + 2 * cols, 0.0 + coff * (-5.0 + 5.0 * (alpha_0 - 3.0)) / pow(h_y, 2)));
				// M_55
				tripletList.push_back(T(i + 5 * rows, j + 5 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_y * h_y)));

			}
			if (j == i - 1 && fmod(i, N_x) != 0) // 主对角线元素下方元素（关于x）
			{
				// M_00
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h_x * h_x)));
				// M_01
				tripletList.push_back(T(i, j + cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_x, 2)));
				// M_04
				tripletList.push_back(T(i, j + 4 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_x, 2)));

				// M_10
				tripletList.push_back(T(i + rows, j, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_x, 2)));
				// M_11
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));
				// M_14
				tripletList.push_back(T(i + rows, j + 4 * cols, 0.0 + coff * (9.0 + 3.0 * (1.0 - alpha_0) - 36.0 * beta_0) / pow(h_x, 2)));

				// M_22
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));
				// M_23
				tripletList.push_back(T(i + 2 * rows, j + 3 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_x, 2)));

				// M_32
				tripletList.push_back(T(i + 3 * rows, j + 2 * cols, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_x, 2)));
				// M_33
				tripletList.push_back(T(i + 3 * rows, j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));

				// M_40
				tripletList.push_back(T(i + 4 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_x * h_x)));
				// M_41
				tripletList.push_back(T(i + 4 * rows, j + cols, 0.0 + coff * (-5.0 + 5.0 * (alpha_0 - 3.0)) / pow(h_x, 2)));
				// M_44
				tripletList.push_back(T(i + 4 * rows, j + 4 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_x * h_x)));

				// M_55
				tripletList.push_back(T(i + 5 * rows, j + 5 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));

			}
			if (j == i) // 主对角线元素
			{
				// M_00
				tripletList.push_back(T(i, j, m_index + coff * (-2.0 * alpha_0 / (h_x * h_x) - 2.0 * alpha_0 / (h_y * h_y))));
				// M_04
				tripletList.push_back(T(i, j + 4 * cols, 0.0 + coff * (-2.0 * alpha_0 + 6.0 - 24.0 * beta_0) / pow(h_x, 2)));
				// M_05 
				tripletList.push_back(T(i, j + 5 * cols, 0.0 + coff * (-2.0 * alpha_0 + 6.0 - 24.0 * beta_0) / pow(h_y, 2)));

				// M_11
				tripletList.push_back(T(i + rows, j + cols, m_index + coff * ((-6.0 + 6.0 * (1.0 - alpha_0)) / pow(h_x, 2) - 2.0 * alpha_0 / pow(h_y, 2))));

				// M_22
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, m_index + coff * ((-6.0 + 6.0 * (1.0 - alpha_0)) / pow(h_y, 2) - 2.0 * alpha_0 / pow(h_x, 2))));

				// M_33
				tripletList.push_back(T(i + 3 * rows, j + 3 * cols, m_index + coff * ((-6.0 + 6.0 * (1.0 - alpha_0)) / pow(h_x, 2) + (-6.0 + 6.0 * (1.0 - alpha_0)) / pow(h_y, 2))));

				// M_40
				tripletList.push_back(T(i + 4 * rows, j, 0.0 + coff * 10.0 * (3.0 - alpha_0) / (h_x * h_x)));
				// M_44
				tripletList.push_back(T(i + 4 * rows, j + 4 * cols, m_index + coff * ((-30.0 + 10.0 * (3.0 - alpha_0) - 120.0 * beta_0) / (h_x * h_x) - 2.0 * alpha_0 / (h_y * h_y))));

				// M_50
				tripletList.push_back(T(i + 5 * rows, j, 0.0 + coff * 10.0 * (3.0 - alpha_0) / (h_y * h_y)));
				// M_55
				tripletList.push_back(T(i + 5 * rows, j + 5 * cols, m_index + coff * ((-30.0 + 10.0 * (3.0 - alpha_0) - 120.0 * beta_0) / (h_y * h_y) - 2.0 * alpha_0 / (h_x * h_x))));

			}
			if (j == i + 1 && fmod(j, N_x) != 0) // 主对角线元素上方元素（关于x）
			{
				// M_00
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h_x * h_x)));
				// M_01
				tripletList.push_back(T(i, j + cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_x, 2)));
				// M_04
				tripletList.push_back(T(i, j + 4 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_x, 2)));

				// M_10
				tripletList.push_back(T(i + rows, j, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_x, 2)));
				// M_11
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));
				// M_14
				tripletList.push_back(T(i + rows, j + 4 * cols, 0.0 + coff * (-9.0 + 3.0 * (alpha_0 - 1.0) + 36.0 * beta_0) / pow(h_x, 2)));

				// M_22
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));
				// M_23
				tripletList.push_back(T(i + 2 * rows, j + 3 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_x, 2)));

				// M_32
				tripletList.push_back(T(i + 3 * rows, j + 2 * cols, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_x, 2)));
				// M_33
				tripletList.push_back(T(i + 3 * rows, j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));

				// M_40
				tripletList.push_back(T(i + 4 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_x * h_x)));
				// M_41
				tripletList.push_back(T(i + 4 * rows, j + cols, 0.0 + coff * (5.0 + 5.0 * (3.0 - alpha_0)) / pow(h_x, 2)));
				// M_44
				tripletList.push_back(T(i + 4 * rows, j + 4 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_x * h_x)));

				// M_55
				tripletList.push_back(T(i + 5 * rows, j + 5 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));

			}
			if (j == i + N_x) // 主对角线元素上方元素（关于y） 
			{
				// M_00
				tripletList.push_back(T(i, j, 0.0 + coff * alpha_0 / (h_y * h_y)));
				// M_02
				tripletList.push_back(T(i, j + 2 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_y, 2)));
				// M_05 
				tripletList.push_back(T(i, j + 5 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_y, 2)));

				// M_11
				tripletList.push_back(T(i + rows, j + cols, 0.0 + coff * alpha_0 / pow(h_y, 2)));
				// M_13
				tripletList.push_back(T(i + rows, j + 3 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_y, 2)));

				// M_20
				tripletList.push_back(T(i + 2 * rows, j, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_y, 2)));
				// M_22
				tripletList.push_back(T(i + 2 * rows, j + 2 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));
				// M_25
				tripletList.push_back(T(i + 2 * rows, j + 5 * cols, 0.0 + coff * (-9.0 + 3.0 * (alpha_0 - 1.0) + 36.0 * beta_0) / pow(h_y, 2)));

				// M_31
				tripletList.push_back(T(i + 3 * rows, j + cols, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_y, 2)));
				// M_33
				tripletList.push_back(T(i + 3 * rows, j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));

				// M_44
				tripletList.push_back(T(i + 4 * rows, j + 4 * cols, 0.0 + coff * alpha_0 / (h_y * h_y)));

				// M_50
				tripletList.push_back(T(i + 5 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_y * h_y)));
				// M_52
				tripletList.push_back(T(i + 5 * rows, j + 2 * cols, 0.0 + coff * (5.0 + 5.0 * (3.0 - alpha_0)) / pow(h_y, 2)));
				// M_55
				tripletList.push_back(T(i + 5 * rows, j + 5 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_y * h_y)));

			}
		}
	}


	for (int i = 0; i < N_y; i++)
	{
		// 主对角线元素下方元素（关于x）
		// M_00
		tripletList.push_back(T(i * N_x, i * N_x + N_x - 1, 0.0 + coff * alpha_0 / (h_x * h_x)));
		// M_01
		tripletList.push_back(T(i * N_x, i * N_x + N_x - 1 + cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_x, 2)));
		// M_04
		tripletList.push_back(T(i * N_x, i * N_x + N_x - 1 + 4 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_x, 2)));


		// M_10
		tripletList.push_back(T(i * N_x + rows, i * N_x + N_x - 1, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_x, 2)));
		// M_11
		tripletList.push_back(T(i * N_x + rows, i * N_x + N_x - 1 + cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));
		// M_14
		tripletList.push_back(T(i * N_x + rows, i * N_x + N_x - 1 + 4 * cols, 0.0 + coff * (9.0 + 3.0 * (1.0 - alpha_0) - 36.0 * beta_0) / pow(h_x, 2)));

		// M_22
		tripletList.push_back(T(i * N_x + 2 * rows, i * N_x + N_x - 1 + 2 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));
		// M_23
		tripletList.push_back(T(i * N_x + 2 * rows, i * N_x + N_x - 1 + 3 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_x, 2)));

		// M_32
		tripletList.push_back(T(i * N_x + 3 * rows, i * N_x + N_x - 1 + 2 * cols, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_x, 2)));
		// M_33
		tripletList.push_back(T(i * N_x + 3 * rows, i * N_x + N_x - 1 + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));

		// M_40
		tripletList.push_back(T(i * N_x + 4 * rows, i * N_x + N_x - 1, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_x * h_x)));
		// M_41
		tripletList.push_back(T(i * N_x + 4 * rows, i * N_x + N_x - 1 + cols, 0.0 + coff * (-5.0 + 5.0 * (alpha_0 - 3.0)) / pow(h_x, 2)));
		// M_44
		tripletList.push_back(T(i * N_x + 4 * rows, i * N_x + N_x - 1 + 4 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_x * h_x)));

		// M_55
		tripletList.push_back(T(i * N_x + 5 * rows, i * N_x + N_x - 1 + 5 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));



		// 主对角线元素上方元素（关于x）
		// M_00
		tripletList.push_back(T(i * N_x + N_x - 1, i * N_x, 0.0 + coff * alpha_0 / (h_x * h_x)));
		// M_01
		tripletList.push_back(T(i * N_x + N_x - 1, i * N_x + cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_x, 2)));
		// M_04
		tripletList.push_back(T(i * N_x + N_x - 1, i * N_x + 4 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_x, 2)));

		// M_10
		tripletList.push_back(T(i * N_x + N_x - 1 + rows, i * N_x, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_x, 2)));
		// M_11
		tripletList.push_back(T(i * N_x + N_x - 1 + rows, i * N_x + cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));
		// M_14
		tripletList.push_back(T(i * N_x + N_x - 1 + rows, i * N_x + 4 * cols, 0.0 + coff * (-9.0 + 3.0 * (alpha_0 - 1.0) + 36.0 * beta_0) / pow(h_x, 2)));

		// M_22
		tripletList.push_back(T(i * N_x + N_x - 1 + 2 * rows, i * N_x + 2 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));
		// M_23
		tripletList.push_back(T(i * N_x + N_x - 1 + 2 * rows, i * N_x + 3 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_x, 2)));

		// M_32
		tripletList.push_back(T(i * N_x + N_x - 1 + 3 * rows, i * N_x + 2 * cols, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_x, 2)));
		// M_33
		tripletList.push_back(T(i * N_x + N_x - 1 + 3 * rows, i * N_x + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_x, 2)));

		// M_40
		tripletList.push_back(T(i * N_x + N_x - 1 + 4 * rows, i * N_x, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_x * h_x)));
		// M_41
		tripletList.push_back(T(i * N_x + N_x - 1 + 4 * rows, i * N_x + cols, 0.0 + coff * (5.0 + 5.0 * (3.0 - alpha_0)) / pow(h_x, 2)));
		// M_44
		tripletList.push_back(T(i * N_x + N_x - 1 + 4 * rows, i * N_x + 4 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_x * h_x)));

		// M_55
		tripletList.push_back(T(i * N_x + N_x - 1 + 5 * rows, i * N_x + 5 * cols, 0.0 + coff * alpha_0 / (h_x * h_x)));

	}

	for (int j = 0; j < N_x; j++)
	{
		// 主对角线元素下方元素（关于y）
		// M_00
		tripletList.push_back(T(j, N_x * (N_y - 1) + j, 0.0 + coff * alpha_0 / (h_y * h_y)));
		// M_02
		tripletList.push_back(T(j, N_x * (N_y - 1) + j + 2 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_y, 2)));
		// M_05 
		tripletList.push_back(T(j, N_x * (N_y - 1) + j + 5 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_y, 2)));

		// M_11
		tripletList.push_back(T(j + rows, N_x * (N_y - 1) + j + cols, 0.0 + coff * alpha_0 / pow(h_y, 2)));
		// M_13
		tripletList.push_back(T(j + rows, N_x * (N_y - 1) + j + 3 * cols, 0.0 + coff * (alpha_0 - 1.0) / pow(h_y, 2)));

		// M_20
		tripletList.push_back(T(j + 2 * rows, N_x * (N_y - 1) + j, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_y, 2)));
		// M_22
		tripletList.push_back(T(j + 2 * rows, N_x * (N_y - 1) + j + 2 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));
		// M_25
		tripletList.push_back(T(j + 2 * rows, N_x * (N_y - 1) + j + 5 * cols, 0.0 + coff * (9.0 + 3.0 * (1.0 - alpha_0) - 36.0 * beta_0) / pow(h_y, 2)));

		// M_31
		tripletList.push_back(T(j + 3 * rows, N_x * (N_y - 1) + j + cols, 0.0 + coff * 3.0 * (1.0 - alpha_0) / pow(h_y, 2)));
		// M_33
		tripletList.push_back(T(j + 3 * rows, N_x * (N_y - 1) + j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));

		// M_44
		tripletList.push_back(T(j + 4 * rows, N_x * (N_y - 1) + j + 4 * cols, 0.0 + coff * alpha_0 / (h_y * h_y)));

		// M_50
		tripletList.push_back(T(j + 5 * rows, N_x * (N_y - 1) + j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_y * h_y)));
		// M_52
		tripletList.push_back(T(j + 5 * rows, N_x * (N_y - 1) + j + 2 * cols, 0.0 + coff * (-5.0 + 5.0 * (alpha_0 - 3.0)) / pow(h_y, 2)));
		// M_55
		tripletList.push_back(T(j + 5 * rows, N_x * (N_y - 1) + j + 5 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_y * h_y)));



		// 主对角线元素上方元素（关于y）
		// M_00
		tripletList.push_back(T(N_x * (N_y - 1) + j, j, 0.0 + coff * alpha_0 / (h_y * h_y)));
		// M_02
		tripletList.push_back(T(N_x * (N_y - 1) + j, j + 2 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_y, 2)));
		// M_05 
		tripletList.push_back(T(N_x * (N_y - 1) + j, j + 5 * cols, 0.0 + coff * (alpha_0 - 3.0 + 12.0 * beta_0) / pow(h_y, 2)));

		// M_11
		tripletList.push_back(T(N_x * (N_y - 1) + j + rows, j + cols, 0.0 + coff * alpha_0 / pow(h_y, 2)));
		// M_13
		tripletList.push_back(T(N_x * (N_y - 1) + j + rows, j + 3 * cols, 0.0 + coff * (1.0 - alpha_0) / pow(h_y, 2)));

		// M_20
		tripletList.push_back(T(N_x * (N_y - 1) + j + 2 * rows, j, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_y, 2)));
		// M_22
		tripletList.push_back(T(N_x * (N_y - 1) + j + 2 * rows, j + 2 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));
		// M_25
		tripletList.push_back(T(N_x * (N_y - 1) + j + 2 * rows, j + 5 * cols, 0.0 + coff * (-9.0 + 3.0 * (alpha_0 - 1.0) + 36.0 * beta_0) / pow(h_y, 2)));

		// M_31
		tripletList.push_back(T(N_x * (N_y - 1) + j + 3 * rows, j + cols, 0.0 + coff * 3.0 * (alpha_0 - 1.0) / pow(h_y, 2)));
		// M_33
		tripletList.push_back(T(N_x * (N_y - 1) + j + 3 * rows, j + 3 * cols, 0.0 + coff * (3.0 + 3.0 * (1.0 - alpha_0)) / pow(h_y, 2)));

		// M_44
		tripletList.push_back(T(N_x * (N_y - 1) + j + 4 * rows, j + 4 * cols, 0.0 + coff * alpha_0 / (h_y * h_y)));

		// M_50
		tripletList.push_back(T(N_x * (N_y - 1) + j + 5 * rows, j, 0.0 + coff * 5.0 * (alpha_0 - 3.0) / (h_y * h_y)));
		// M_52
		tripletList.push_back(T(N_x * (N_y - 1) + j + 5 * rows, j + 2 * cols, 0.0 + coff * (5.0 + 5.0 * (3.0 - alpha_0)) / pow(h_y, 2)));
		// M_55
		tripletList.push_back(T(N_x * (N_y - 1) + j + 5 * rows, j + 5 * cols, 0.0 + coff * (-15.0 + 5.0 * (alpha_0 - 3.0) + 60.0 * beta_0) / (h_y * h_y)));

	}

	mat.setFromTriplets(tripletList.begin(), tripletList.end());

}


// time step
double comput_rho_min(double** rho)
{
	double temp_min = 10.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			if (rho[i][j] < temp_min)
			{
				temp_min = rho[i][j];
			}
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

void compute_st(double** s_j, double** rho_0, double** m1_0, double** m2_0, double** Ee_0)
{
	double temp_1 = a * omega_ei * pow(C_ve, 4) + a * omega_ei * pow(C_ve, 3) * C_vi;
	double temp_2 = 4.0 * a * omega_ei * omega_er * pow(C_ve, 4) * C_vi;

	double** e_e; createMatrix(N_y, N_x, e_e); double** s1; createMatrix(N_y, N_x, s1); double** s2; createMatrix(N_y, N_x, s2);
	double** s; createMatrix(N_y, N_x, s); double** alp_1; createMatrix(N_y, N_x, alp_1); double** alp_2; createMatrix(N_y, N_x, alp_2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_rhou2 = pow(m1_0[i + 1][j + 1], 2) / rho_0[i + 1][j + 1];
			double temp_rhov2 = pow(m2_0[i + 1][j + 1], 2) / rho_0[i + 1][j + 1];
			e_e[i][j] = (Ee_0[i + 1][j + 1] - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / rho_0[i + 1][j + 1];
		}
	}

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_s1_1 = temp_1 + 4.0 * a * omega_er * C_vi * pow(e_e[i][j], 3) + omega_er * pow(C_ve, 4) * C_vi * rho_0[i + 1][j + 1];
			double temp_s1_2 = temp_2 * (4.0 * a * pow(e_e[i][j], 3) + pow(C_ve, 4) * rho_0[i + 1][j + 1] + pow(C_ve, 3) * C_vi * rho_0[i + 1][j + 1]);
			s1[i][j] = sqrt(pow(temp_s1_1, 2) - temp_s1_2);
			s2[i][j] = -1.0 * temp_s1_1;
			s[i][j] = 2.0 * a * pow(C_ve, 4) * C_vi * rho_0[i + 1][j + 1];
			alp_1[i][j] = (s2[i][j] - s1[i][j]) / s[i][j];
			alp_2[i][j] = (s2[i][j] + s1[i][j]) / s[i][j];
		}
	}

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_max = abs(alp_1[i][j]);
			if (abs(alp_2[i][j]) > temp_max)
			{
				temp_max = abs(alp_2[i][j]);
			}
			s_j[i][j] = temp_max;
		}
	}

	deleteMatrix(N_y, e_e); deleteMatrix(N_y, s1); deleteMatrix(N_y, s2); deleteMatrix(N_y, s);
	deleteMatrix(N_y, alp_1); deleteMatrix(N_y, alp_2);
}

double compute_sj_max(double** s_j)
{
	double sj_max = 0.0;
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x; j++)
		{
			if (s_j[i][j] > sj_max)
			{
				sj_max = s_j[i][j];
			}
		}
	}

	return sj_max;
}

