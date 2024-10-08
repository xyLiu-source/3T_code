#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG2D.h"

// 2D Modal DG general framework
void get_Gaussian_points_x(double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* x, double h_x)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		X_1[i] = x[i] + h_x / 2.0 * GQxi_1; X_2[i] = x[i] + h_x / 2.0 * GQxi_2; X_3[i] = x[i] + h_x / 2.0 * GQxi_3;
		X_4[i] = x[i] + h_x / 2.0 * GQxi_4; X_5[i] = x[i] + h_x / 2.0 * GQxi_5;
	}
}

void get_Gaussian_points_y(double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double* y, double h_y)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		Y_1[i] = y[i] + h_y / 2.0 * GQxi_1; Y_2[i] = y[i] + h_y / 2.0 * GQxi_2; Y_3[i] = y[i] + h_y / 2.0 * GQxi_3;
		Y_4[i] = y[i] + h_y / 2.0 * GQxi_4; Y_5[i] = y[i] + h_y / 2.0 * GQxi_5;
	}
}

void get_Gauss_quadrature_2D(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** u_11, double** u_21, double** u_31, double** u_41, double** u_51, double** u_12, double** u_22, double** u_32, double** u_42, double** u_52, double** u_13, double** u_23, double** u_33, double** u_43, double** u_53, double** u_14, double** u_24, double** u_34, double** u_44, double** u_54, double** u_15, double** u_25, double** u_35, double** u_45, double** u_55, double h_x, double h_y)
{
	double w_11 = weight_1 * weight_1; double w_12 = weight_1 * weight_2; double w_13 = weight_1 * weight_3;
	double w_14 = weight_1 * weight_4; double w_15 = weight_1 * weight_5;

	double w_21 = weight_2 * weight_1; double w_22 = weight_2 * weight_2; double w_23 = weight_2 * weight_3;
	double w_24 = weight_2 * weight_4; double w_25 = weight_2 * weight_5;

	double w_31 = weight_3 * weight_1; double w_32 = weight_3 * weight_2; double w_33 = weight_3 * weight_3;
	double w_34 = weight_3 * weight_4; double w_35 = weight_3 * weight_5;

	double w_41 = weight_4 * weight_1; double w_42 = weight_4 * weight_2; double w_43 = weight_4 * weight_3;
	double w_44 = weight_4 * weight_4; double w_45 = weight_4 * weight_5;

	double w_51 = weight_5 * weight_1; double w_52 = weight_5 * weight_2; double w_53 = weight_5 * weight_3;
	double w_54 = weight_5 * weight_4; double w_55 = weight_5 * weight_5;

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] + w_12 * u_12[i][j] + w_13 * u_13[i][j] + w_14 * u_14[i][j] + w_15 * u_15[i][j] + w_21 * u_21[i][j] + w_22 * u_22[i][j] + w_23 * u_23[i][j] + w_24 * u_24[i][j] + w_25 * u_25[i][j] + w_31 * u_31[i][j] + w_32 * u_32[i][j] + w_33 * u_33[i][j] + w_34 * u_34[i][j] + w_35 * u_35[i][j] + w_41 * u_41[i][j] + w_42 * u_42[i][j] + w_43 * u_43[i][j] + w_44 * u_44[i][j] + w_45 * u_45[i][j] + w_51 * u_51[i][j] + w_52 * u_52[i][j] + w_53 * u_53[i][j] + w_54 * u_54[i][j] + w_55 * u_55[i][j]);

			U_1[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_1 + w_13 * u_13[i][j] * GQxi_1 + w_14 * u_14[i][j] * GQxi_1 + w_15 * u_15[i][j] * GQxi_1 + w_21 * u_21[i][j] * GQxi_2 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_2 + w_24 * u_24[i][j] * GQxi_2 + w_25 * u_25[i][j] * GQxi_2 + w_31 * u_31[i][j] * GQxi_3 + w_32 * u_32[i][j] * GQxi_3 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_3 + w_35 * u_35[i][j] * GQxi_3 + w_41 * u_41[i][j] * GQxi_4 + w_42 * u_42[i][j] * GQxi_4 + w_43 * u_43[i][j] * GQxi_4 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_4 + w_51 * u_51[i][j] * GQxi_5 + w_52 * u_52[i][j] * GQxi_5 + w_53 * u_53[i][j] * GQxi_5 + w_54 * u_54[i][j] * GQxi_5 + w_55 * u_55[i][j] * GQxi_5);

			U_2[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_2 + w_13 * u_13[i][j] * GQxi_3 + w_14 * u_14[i][j] * GQxi_4 + w_15 * u_15[i][j] * GQxi_5 + w_21 * u_21[i][j] * GQxi_1 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_3 + w_24 * u_24[i][j] * GQxi_4 + w_25 * u_25[i][j] * GQxi_5 + w_31 * u_31[i][j] * GQxi_1 + w_32 * u_32[i][j] * GQxi_2 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_4 + w_35 * u_35[i][j] * GQxi_5 + w_41 * u_41[i][j] * GQxi_1 + w_42 * u_42[i][j] * GQxi_2 + w_43 * u_43[i][j] * GQxi_3 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_5 + w_51 * u_51[i][j] * GQxi_1 + w_52 * u_52[i][j] * GQxi_2 + w_53 * u_53[i][j] * GQxi_3 + w_54 * u_54[i][j] * GQxi_4 + w_55 * u_55[i][j] * GQxi_5);

			U_3[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] * GQxi_1 * GQxi_1 + w_12 * u_12[i][j] * GQxi_1 * GQxi_2 + w_13 * u_13[i][j] * GQxi_1 * GQxi_3 + w_14 * u_14[i][j] * GQxi_1 * GQxi_4 + w_15 * u_15[i][j] * GQxi_1 * GQxi_5 + w_21 * u_21[i][j] * GQxi_2 * GQxi_1 + w_22 * u_22[i][j] * GQxi_2 * GQxi_2 + w_23 * u_23[i][j] * GQxi_2 * GQxi_3 + w_24 * u_24[i][j] * GQxi_2 * GQxi_4 + w_25 * u_25[i][j] * GQxi_2 * GQxi_5 + w_31 * u_31[i][j] * GQxi_3 * GQxi_1 + w_32 * u_32[i][j] * GQxi_3 * GQxi_2 + w_33 * u_33[i][j] * GQxi_3 * GQxi_3 + w_34 * u_34[i][j] * GQxi_3 * GQxi_4 + w_35 * u_35[i][j] * GQxi_3 * GQxi_5 + w_41 * u_41[i][j] * GQxi_4 * GQxi_1 + w_42 * u_42[i][j] * GQxi_4 * GQxi_2 + w_43 * u_43[i][j] * GQxi_4 * GQxi_3 + w_44 * u_44[i][j] * GQxi_4 * GQxi_4 + w_45 * u_45[i][j] * GQxi_4 * GQxi_5 + w_51 * u_51[i][j] * GQxi_5 * GQxi_1 + w_52 * u_52[i][j] * GQxi_5 * GQxi_2 + w_53 * u_53[i][j] * GQxi_5 * GQxi_3 + w_54 * u_54[i][j] * GQxi_5 * GQxi_4 + w_55 * u_55[i][j] * GQxi_5 * GQxi_5);

			U_4[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_12 * u_12[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_13 * u_13[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_14 * u_14[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_15 * u_15[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_21 * u_21[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_22 * u_22[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_23 * u_23[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_24 * u_24[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_25 * u_25[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_31 * u_31[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_32 * u_32[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_33 * u_33[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_34 * u_34[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_35 * u_35[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_41 * u_41[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_42 * u_42[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_43 * u_43[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_44 * u_44[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_45 * u_45[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_51 * u_51[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_52 * u_52[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_53 * u_53[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_54 * u_54[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_55 * u_55[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0);

			U_5[i][j] = h_x * h_y / 4.0 * (w_11 * u_11[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_12 * u_12[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_13 * u_13[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_14 * u_14[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_15 * u_15[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_21 * u_21[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_22 * u_22[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_23 * u_23[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_24 * u_24[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_25 * u_25[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_31 * u_31[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_32 * u_32[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_33 * u_33[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_34 * u_34[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_35 * u_35[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_41 * u_41[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_42 * u_42[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_43 * u_43[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_44 * u_44[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_45 * u_45[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 + w_51 * u_51[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + w_52 * u_52[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + w_53 * u_53[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + w_54 * u_54[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + w_55 * u_55[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0);
		}
	}
}

void get_L2_projection_2D(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double h_x, double h_y)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = u_0[i][j] / (h_x * h_y);
			U_1[i][j] = 3.0 * u_1[i][j] / (h_x * h_y); U_2[i][j] = 3.0 * u_2[i][j] / (h_x * h_y);
			U_3[i][j] = 9.0 * u_3[i][j] / (h_x * h_y);
			U_4[i][j] = 5.0 * u_4[i][j] / (h_x * h_y); U_5[i][j] = 5.0 * u_5[i][j] / (h_x * h_y);
		}
	}
}

void createMatrix(int rows, int cols, double**& arr)
{
	arr = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (int j = 0; j < cols; j++)
		{
			arr[i][j] = 0.0;
		}
	}
}

void deleteMatrix(int rows, double** arr)
{
	for (int i = 0; i < rows; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;
}

void Matrix_to_vector(double* v, double** M, int M_rows, int M_cols)
{
	for (int i = 1; i < M_rows - 1; i++)
	{
		for (int j = 1; j < M_cols - 1; j++)
		{
			v[(i - 1) * (M_cols - 2) + (j - 1)] = M[i][j];
		}
	}
}

void vector_to_Matrix(double** M, double* v, int M_rows, int M_cols)
{
	for (int i = 1; i < M_rows - 1; i++)
	{
		for (int j = 1; j < M_cols - 1; j++)
		{
			M[i][j] = v[(i - 1) * (M_cols - 2) + (j - 1)];
		}
	}
}

void update_boundary(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5)
{
	// periodic
	for (int i = 1; i < N_y + 1; i++)
	{
		// left boundary
		u_0[i][0] = u_0[i][N_x]; u_1[i][0] = u_1[i][N_x]; u_2[i][0] = u_2[i][N_x];
		u_3[i][0] = u_3[i][N_x]; u_4[i][0] = u_4[i][N_x]; u_5[i][0] = u_5[i][N_x];
		// right boundary
		u_0[i][N_x + 1] = u_0[i][1]; u_1[i][N_x + 1] = u_1[i][1]; u_2[i][N_x + 1] = u_2[i][1];
		u_3[i][N_x + 1] = u_3[i][1]; u_4[i][N_x + 1] = u_4[i][1]; u_5[i][N_x + 1] = u_5[i][1];
	}

	for (int j = 1; j < N_x + 1; j++)
	{
		// Top boundary
		u_0[0][j] = u_0[N_y][j]; u_1[0][j] = u_1[N_y][j]; u_2[0][j] = u_2[N_y][j];
		u_3[0][j] = u_3[N_y][j]; u_4[0][j] = u_4[N_y][j]; u_5[0][j] = u_5[N_y][j];
		// Bottom boundary
		u_0[N_y + 1][j] = u_0[1][j]; u_1[N_y + 1][j] = u_1[1][j]; u_2[N_y + 1][j] = u_2[1][j];
		u_3[N_y + 1][j] = u_3[1][j]; u_4[N_y + 1][j] = u_4[1][j]; u_5[N_y + 1][j] = u_5[1][j];
	}
}

double compute_MAX(double** Matrix)
{
	double temp_max = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			if (abs(Matrix[i][j]) > temp_max)
			{
				temp_max = abs(Matrix[i][j]);
			}
		}
	}
	return temp_max;
}


// accuracy solution
double rho_ex(double x, double y, double t)
{
	return 1.0 + 0.2 * sin(x + y - t);
}

double u_ex(double x, double y, double t)
{
	return 0.7;
}

double v_ex(double x, double y, double t)
{
	return 0.3;
}

double rhoe_e_ex(double x, double y, double t)
{
	return 3.0;
}

double rhoe_i_ex(double x, double y, double t)
{
	return 3.0;
}

double rhoe_r_ex(double x, double y, double t)
{
	return 2.0;
}

double rho_prime_ex(double x, double y, double t)
{
	return 0.2 * cos(x + y - t);
}

double rho_prime2_ex(double x, double y, double t)
{
	return -0.2 * sin(x + y - t);
}


// model solution
double m1_ex(double x, double y, double t)
{
	return rho_ex(x, y, t) * u_ex(x, y, t);
}

double m2_ex(double x, double y, double t)
{
	return rho_ex(x, y, t) * v_ex(x, y, t);
}

double Ee_ex(double x, double y, double t)
{
	return rhoe_e_ex(x, y, t) + rho_ex(x, y, t) * u_ex(x, y, t) * u_ex(x, y, t) / 6.0 + rho_ex(x, y, t) * v_ex(x, y, t) * v_ex(x, y, t) / 6.0;
}

double Ei_ex(double x, double y, double t)
{
	return rhoe_i_ex(x, y, t) + rho_ex(x, y, t) * u_ex(x, y, t) * u_ex(x, y, t) / 6.0 + rho_ex(x, y, t) * v_ex(x, y, t) * v_ex(x, y, t) / 6.0;
}

double Er_ex(double x, double y, double t)
{
	return rhoe_r_ex(x, y, t) + rho_ex(x, y, t) * u_ex(x, y, t) * u_ex(x, y, t) / 6.0 + rho_ex(x, y, t) * v_ex(x, y, t) * v_ex(x, y, t) / 6.0;
}

double S3_ex(double x, double y, double t)
{
	double Tr4 = rhoe_r_ex(x, y, t) / a;
	double Te = rhoe_e_ex(x, y, t) / (C_ve * rho_ex(x, y, t));
	double Ti = rhoe_i_ex(x, y, t) / (C_vi * rho_ex(x, y, t));
	double Source_i = omega_ei * (Te - Ti);
	double Source_r = omega_er * (pow(Te, 4) - Tr4);
	double Source_e = Source_i + Source_r;
	double D1 = 2.0 * pow(rho_prime_ex(x, y, t), 2) / pow(rho_ex(x, y, t), 3) - rho_prime2_ex(x, y, t) / pow(rho_ex(x, y, t), 2);
	double Diff_e_x = kappa_e / C_ve * D1 * rhoe_e_ex(x, y, t);
	double Diff_e_y = kappa_e / C_ve * D1 * rhoe_e_ex(x, y, t);
	double S3 = Source_e - Diff_e_x - Diff_e_y;
	return S3;
}

double S4_ex(double x, double y, double t)
{
	double Te = rhoe_e_ex(x, y, t) / (C_ve * rho_ex(x, y, t));
	double Ti = rhoe_i_ex(x, y, t) / (C_vi * rho_ex(x, y, t));
	double Source_i = omega_ei * (Te - Ti);
	double D1 = 2.0 * pow(rho_prime_ex(x, y, t), 2) / pow(rho_ex(x, y, t), 3) - rho_prime2_ex(x, y, t) / pow(rho_ex(x, y, t), 2);
	double Diff_i_x = kappa_i / C_vi * D1 * rhoe_i_ex(x, y, t);
	double Diff_i_y = kappa_i / C_vi * D1 * rhoe_i_ex(x, y, t);
	double S4 = -1.0 * Source_i - Diff_i_x - Diff_i_y;
	return S4;
}

double S5_ex(double x, double y, double t)
{
	double Tr4 = rhoe_r_ex(x, y, t) / a;
	double Te = rhoe_e_ex(x, y, t) / (C_ve * rho_ex(x, y, t));
	double Source_r = omega_er * (pow(Te, 4) - Tr4);
	double S5 = -1.0 * Source_r;
	return S5;
}


// Compute error and order
double* rho_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(rho_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(rho_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(rho_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* m1_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(m1_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(m1_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(m1_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* m2_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(m2_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(m2_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(m2_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Ee_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(Ee_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(Ee_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(Ee_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Ei_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(Ei_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(Ei_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(Ei_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Er_err_2D(double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double* x, double* y, double h_x, double h_y, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			int n = 9; // the number of points selected in each cell along x and y directions
			double err = 0.0;
			for (int p = 0; p < n; p++)
			{
				for (int q = 0; q < n; q++)
				{
					switch (Mode)
					{
					case 0:
						err = abs(Er_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j]));
						break;
					case 1:
						err = abs(Er_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0)));
						break;
					case 2:
						err = abs(Er_ex(x[j] + h_x / 2.0 * (2.0 / (n - 1.0) * p - 1.0), y[i] + h_y / 2.0 * (2.0 / (n - 1.0) * q - 1.0), t) - (u_0[i][j] + u_1[i][j] * (2.0 / (n - 1.0) * p - 1.0) + u_2[i][j] * (2.0 / (n - 1.0) * q - 1.0) + u_3[i][j] * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * q - 1.0) + u_4[i][j] * (3.0 * (2.0 / (n - 1.0) * p - 1.0) * (2.0 / (n - 1.0) * p - 1.0) - 1.0) / 2.0 + u_5[i][j] * (3.0 * (2.0 / (n - 1.0) * q - 1.0) * (2.0 / (n - 1.0) * q - 1.0) - 1.0) / 2.0));
						break;
					default:
						cout << "please choose the correct mode !" << endl;
						break;
					}
					L1_error = L1_error + err / pow(n, 2);
					L2_error = L2_error + pow(err, 2) / pow(n, 2);
					if (err > Linf_error)
					{
						Linf_error = err;
					}
				}
			}
		}
	}
	L1_error = L1_error / (N_x * N_y);
	L2_error = sqrt(L2_error / (N_x * N_y));
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}
