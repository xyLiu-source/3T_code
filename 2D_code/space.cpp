#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG2D.h"
#include "space.h"

// space function
//Source
void compute_Talpha(double** u, double xi, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / (C_valpha * temp_rho);
		}
	}
}

void compute_Tr4(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Er = Er_0[i][j] + Er_1[i][j] * xi + Er_2[i][j] * eta + Er_3[i][j] * xi * eta + Er_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Er_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = (temp_Er - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / a;
		}
	}
}

void F_Si(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5)
{
	double** Te; createMatrix(N_y + 2, N_x + 2, Te); double** Ti; createMatrix(N_y + 2, N_x + 2, Ti);
	compute_Talpha(Te, xi, eta, C_ve, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	compute_Talpha(Ti, xi, eta, C_vi, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u[i][j] = omega_ei * (Te[i][j] - Ti[i][j]);
		}
	}
	deleteMatrix(N_y + 2, Te); deleteMatrix(N_y + 2, Ti);
}

void F_Sr(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Te; createMatrix(N_y + 2, N_x + 2, Te); double** Tr4; createMatrix(N_y + 2, N_x + 2, Tr4);
	compute_Talpha(Te, xi, eta, C_ve, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	compute_Tr4(Tr4, xi, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u[i][j] = omega_er * (pow(Te[i][j], 4) - Tr4[i][j]);
		}
	}
	deleteMatrix(N_y + 2, Te); deleteMatrix(N_y + 2, Tr4);
}

void Source_i(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double h_x, double h_y)
{
	int rows = N_y + 2; int cols = N_x + 2;
	double** FSi_11; createMatrix(rows, cols, FSi_11); double** FSi_12; createMatrix(rows, cols, FSi_12); double** FSi_13; createMatrix(rows, cols, FSi_13); double** FSi_14; createMatrix(rows, cols, FSi_14); double** FSi_15; createMatrix(rows, cols, FSi_15);
	double** FSi_21; createMatrix(rows, cols, FSi_21); double** FSi_22; createMatrix(rows, cols, FSi_22); double** FSi_23; createMatrix(rows, cols, FSi_23); double** FSi_24; createMatrix(rows, cols, FSi_24); double** FSi_25; createMatrix(rows, cols, FSi_25);
	double** FSi_31; createMatrix(rows, cols, FSi_31); double** FSi_32; createMatrix(rows, cols, FSi_32); double** FSi_33; createMatrix(rows, cols, FSi_33); double** FSi_34; createMatrix(rows, cols, FSi_34); double** FSi_35; createMatrix(rows, cols, FSi_35);
	double** FSi_41; createMatrix(rows, cols, FSi_41); double** FSi_42; createMatrix(rows, cols, FSi_42); double** FSi_43; createMatrix(rows, cols, FSi_43); double** FSi_44; createMatrix(rows, cols, FSi_44); double** FSi_45; createMatrix(rows, cols, FSi_45);
	double** FSi_51; createMatrix(rows, cols, FSi_51); double** FSi_52; createMatrix(rows, cols, FSi_52); double** FSi_53; createMatrix(rows, cols, FSi_53); double** FSi_54; createMatrix(rows, cols, FSi_54); double** FSi_55; createMatrix(rows, cols, FSi_55);
	F_Si(FSi_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	F_Si(FSi_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	F_Si(FSi_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	F_Si(FSi_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	F_Si(FSi_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	F_Si(FSi_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	get_Gauss_quadrature_2D(U_0, U_1, U_2, U_3, U_4, U_5, FSi_11, FSi_21, FSi_31, FSi_41, FSi_51, FSi_12, FSi_22, FSi_32, FSi_42, FSi_52, FSi_13, FSi_23, FSi_33, FSi_43, FSi_53, FSi_14, FSi_24, FSi_34, FSi_44, FSi_54, FSi_15, FSi_25, FSi_35, FSi_45, FSi_55, h_x, h_y);
	deleteMatrix(rows, FSi_11); deleteMatrix(rows, FSi_12); deleteMatrix(rows, FSi_13); deleteMatrix(rows, FSi_14); deleteMatrix(rows, FSi_15);
	deleteMatrix(rows, FSi_21); deleteMatrix(rows, FSi_22); deleteMatrix(rows, FSi_23); deleteMatrix(rows, FSi_24); deleteMatrix(rows, FSi_25);
	deleteMatrix(rows, FSi_31); deleteMatrix(rows, FSi_32); deleteMatrix(rows, FSi_33); deleteMatrix(rows, FSi_34); deleteMatrix(rows, FSi_35);
	deleteMatrix(rows, FSi_41); deleteMatrix(rows, FSi_42); deleteMatrix(rows, FSi_43); deleteMatrix(rows, FSi_44); deleteMatrix(rows, FSi_45);
	deleteMatrix(rows, FSi_51); deleteMatrix(rows, FSi_52); deleteMatrix(rows, FSi_53); deleteMatrix(rows, FSi_54); deleteMatrix(rows, FSi_55);
}

void Source_r(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	int rows = N_y + 2; int cols = N_x + 2;
	double** FSr_11; createMatrix(rows, cols, FSr_11); double** FSr_12; createMatrix(rows, cols, FSr_12); double** FSr_13; createMatrix(rows, cols, FSr_13); double** FSr_14; createMatrix(rows, cols, FSr_14); double** FSr_15; createMatrix(rows, cols, FSr_15);
	double** FSr_21; createMatrix(rows, cols, FSr_21); double** FSr_22; createMatrix(rows, cols, FSr_22); double** FSr_23; createMatrix(rows, cols, FSr_23); double** FSr_24; createMatrix(rows, cols, FSr_24); double** FSr_25; createMatrix(rows, cols, FSr_25);
	double** FSr_31; createMatrix(rows, cols, FSr_31); double** FSr_32; createMatrix(rows, cols, FSr_32); double** FSr_33; createMatrix(rows, cols, FSr_33); double** FSr_34; createMatrix(rows, cols, FSr_34); double** FSr_35; createMatrix(rows, cols, FSr_35);
	double** FSr_41; createMatrix(rows, cols, FSr_41); double** FSr_42; createMatrix(rows, cols, FSr_42); double** FSr_43; createMatrix(rows, cols, FSr_43); double** FSr_44; createMatrix(rows, cols, FSr_44); double** FSr_45; createMatrix(rows, cols, FSr_45);
	double** FSr_51; createMatrix(rows, cols, FSr_51); double** FSr_52; createMatrix(rows, cols, FSr_52); double** FSr_53; createMatrix(rows, cols, FSr_53); double** FSr_54; createMatrix(rows, cols, FSr_54); double** FSr_55; createMatrix(rows, cols, FSr_55);

	F_Sr(FSr_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_Sr(FSr_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_Sr(FSr_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_Sr(FSr_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_Sr(FSr_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F_Sr(FSr_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	get_Gauss_quadrature_2D(U_0, U_1, U_2, U_3, U_4, U_5, FSr_11, FSr_21, FSr_31, FSr_41, FSr_51, FSr_12, FSr_22, FSr_32, FSr_42, FSr_52, FSr_13, FSr_23, FSr_33, FSr_43, FSr_53, FSr_14, FSr_24, FSr_34, FSr_44, FSr_54, FSr_15, FSr_25, FSr_35, FSr_45, FSr_55, h_x, h_y);
	deleteMatrix(rows, FSr_11); deleteMatrix(rows, FSr_12); deleteMatrix(rows, FSr_13); deleteMatrix(rows, FSr_14); deleteMatrix(rows, FSr_15);
	deleteMatrix(rows, FSr_21); deleteMatrix(rows, FSr_22); deleteMatrix(rows, FSr_23); deleteMatrix(rows, FSr_24); deleteMatrix(rows, FSr_25);
	deleteMatrix(rows, FSr_31); deleteMatrix(rows, FSr_32); deleteMatrix(rows, FSr_33); deleteMatrix(rows, FSr_34); deleteMatrix(rows, FSr_35);
	deleteMatrix(rows, FSr_41); deleteMatrix(rows, FSr_42); deleteMatrix(rows, FSr_43); deleteMatrix(rows, FSr_44); deleteMatrix(rows, FSr_45);
	deleteMatrix(rows, FSr_51); deleteMatrix(rows, FSr_52); deleteMatrix(rows, FSr_53); deleteMatrix(rows, FSr_54); deleteMatrix(rows, FSr_55);
}

void Source_e(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** Si_0, double** Si_1, double** Si_2, double** Si_3, double** Si_4, double** Si_5, double** Sr_0, double** Sr_1, double** Sr_2, double** Sr_3, double** Sr_4, double** Sr_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = -1.0 * (Si_0[i][j] + Sr_0[i][j]);
			U_1[i][j] = -1.0 * (Si_1[i][j] + Sr_1[i][j]);
			U_2[i][j] = -1.0 * (Si_2[i][j] + Sr_2[i][j]);
			U_3[i][j] = -1.0 * (Si_3[i][j] + Sr_3[i][j]);
			U_4[i][j] = -1.0 * (Si_4[i][j] + Sr_4[i][j]);
			U_5[i][j] = -1.0 * (Si_5[i][j] + Sr_5[i][j]);
		}
	}
}


//Convection
void F_d(double** u, double xi, double eta, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u[i][j] = u_0[i][j] + u_1[i][j] * xi + u_2[i][j] * eta + u_3[i][j] * xi * eta + u_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + u_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
		}
	}
}

void F1_m1_fun(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double g_t = gamma_e + gamma_i + gamma_r - 3.0;
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ee = Ee_0[i][j] + Ee_1[i][j] * xi + Ee_2[i][j] * eta + Ee_3[i][j] * xi * eta + Ee_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ee_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ei = Ei_0[i][j] + Ei_1[i][j] * xi + Ei_2[i][j] * eta + Ei_3[i][j] * xi * eta + Ei_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ei_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Er = Er_0[i][j] + Er_1[i][j] * xi + Er_2[i][j] * eta + Er_3[i][j] * xi * eta + Er_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Er_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = temp_rhou2 + (gamma_e - 1.0) * temp_Ee + (gamma_i - 1.0) * temp_Ei + (gamma_r - 1.0) * temp_Er - g_t * (temp_rhou2 + temp_rhov2) / 6.0;
		}
	}
}

void F_m(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = temp_m1 * temp_m2 / temp_rho;
		}
	}
}

void F2_m2_fun(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double g_t = gamma_e + gamma_i + gamma_r - 3.0;
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ee = Ee_0[i][j] + Ee_1[i][j] * xi + Ee_2[i][j] * eta + Ee_3[i][j] * xi * eta + Ee_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ee_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ei = Ei_0[i][j] + Ei_1[i][j] * xi + Ei_2[i][j] * eta + Ei_3[i][j] * xi * eta + Ei_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ei_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Er = Er_0[i][j] + Er_1[i][j] * xi + Er_2[i][j] * eta + Er_3[i][j] * xi * eta + Er_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Er_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = temp_rhov2 + (gamma_e - 1.0) * temp_Ee + (gamma_i - 1.0) * temp_Ei + (gamma_r - 1.0) * temp_Er - g_t * (temp_rhou2 + temp_rhov2) / 6.0;
		}
	}
}

void F1_alpha(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double gamma_alpha)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = (temp_Ealpha + (gamma_alpha - 1.0) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0)) * temp_m1 / temp_rho;
		}
	}
}

void F2_alpha(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double gamma_alpha)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = (temp_Ealpha + (gamma_alpha - 1.0) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0)) * temp_m2 / temp_rho;
		}
	}
}

void compute_jump_x(double** u, double** F_plus_x, double** F_minus_x)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x + 2; j++)
		{
			u[i][j] = F_plus_x[i + 1][j] - F_minus_x[i + 1][j - 1];
		}
		u[i][0] = 0.0; u[i][N_x + 2] = 0.0;
	}
}

void compute_ave_x(double** u, double** F_plus_x, double** F_minus_x)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x + 2; j++)
		{
			u[i][j] = 0.5 * (F_plus_x[i + 1][j] + F_minus_x[i + 1][j - 1]);
		}
		u[i][0] = 0.0; u[i][N_x + 2] = 0.0;
	}
}

void compute_jump_y(double** u, double** F_plus_y, double** F_minus_y)
{
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 1; i < N_y + 2; i++)
		{
			u[i][j] = F_plus_y[i][j + 1] - F_minus_y[i - 1][j + 1];
		}
		u[0][j] = 0.0; u[N_y + 2][j] = 0.0;
	}
}

void compute_ave_y(double** u, double** F_plus_y, double** F_minus_y)
{
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 1; i < N_y + 2; i++)
		{
			u[i][j] = 0.5 * (F_plus_y[i][j + 1] + F_minus_y[i - 1][j + 1]);
		}
		u[0][j] = 0.0; u[N_y + 2][j] = 0.0;
	}
}

void C_s(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ee = Ee_0[i][j] + Ee_1[i][j] * xi + Ee_2[i][j] * eta + Ee_3[i][j] * xi * eta + Ee_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ee_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ei = Ei_0[i][j] + Ei_1[i][j] * xi + Ei_2[i][j] * eta + Ei_3[i][j] * xi * eta + Ei_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ei_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Er = Er_0[i][j] + Er_1[i][j] * xi + Er_2[i][j] * eta + Er_3[i][j] * xi * eta + Er_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Er_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double e_e = (temp_Ee - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / temp_rho;
			double e_i = (temp_Ei - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / temp_rho;
			double e_r = (temp_Er - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / temp_rho;
			u[i][j] = sqrt(gamma_e * (gamma_e - 1.0) * e_e + gamma_i * (gamma_i - 1.0) * e_i + gamma_r * (gamma_r - 1.0) * e_r);
		}
	}
}

void compute_upxi(double** u, double xi, double eta, double** u_1, double** u_3, double** u_4)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u[i][j] = u_1[i][j] + u_3[i][j] * eta + 3.0 * u_4[i][j] * xi;
		}
	}
}

void compute_upeta(double** u, double xi, double eta, double** u_2, double** u_3, double** u_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u[i][j] = u_2[i][j] + u_3[i][j] * xi + 3.0 * u_5[i][j] * eta;
		}
	}
}

void compute_F1_alphaC(double** alpha_c, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** m1_plus; createMatrix(N_y + 2, N_x + 2, m1_plus); double** m1_minus; createMatrix(N_y + 2, N_x + 2, m1_minus);
	double** rho_plus; createMatrix(N_y + 2, N_x + 2, rho_plus); double** rho_minus; createMatrix(N_y + 2, N_x + 2, rho_minus);
	F_d(m1_plus, -1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(m1_minus, 1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(rho_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5); F_d(rho_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u_plus[i][j] = m1_plus[i][j] / rho_plus[i][j];
			u_minus[i][j] = m1_minus[i][j] / rho_minus[i][j];
		}
	}
	double** Cs_plus; createMatrix(N_y + 2, N_x + 2, Cs_plus); double** Cs_minus; createMatrix(N_y + 2, N_x + 2, Cs_minus);
	C_s(Cs_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	C_s(Cs_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	//alpha_c(N_y, N_x + 1)
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 2; j++)
		{
			double temp_max = abs(u_plus[i][j] - Cs_plus[i][j]);
			if (abs(u_plus[i][j]) > temp_max)
			{
				temp_max = abs(u_plus[i][j]);
			}
			if (abs(u_plus[i][j] + Cs_plus[i][j]) > temp_max)
			{
				temp_max = abs(u_plus[i][j] + Cs_plus[i][j]);
			}
			if (abs(u_minus[i][j - 1] - Cs_minus[i][j - 1]) > temp_max)
			{
				temp_max = abs(u_minus[i][j - 1] - Cs_minus[i][j - 1]);
			}
			if (abs(u_minus[i][j - 1]) > temp_max)
			{
				temp_max = abs(u_minus[i][j - 1]);
			}
			if (abs(u_minus[i][j - 1] + Cs_minus[i][j - 1]) > temp_max)
			{
				temp_max = abs(u_minus[i][j - 1] + Cs_minus[i][j - 1]);
			}
			alpha_c[i - 1][j - 1] = temp_max;
		}
	}

	deleteMatrix(N_y + 2, m1_plus); deleteMatrix(N_y + 2, m1_minus);
	deleteMatrix(N_y + 2, rho_plus); deleteMatrix(N_y + 2, rho_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 2, Cs_plus); deleteMatrix(N_y + 2, Cs_minus);
}

void compute_F2_alphaC(double** alpha_c, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** m2_plus; createMatrix(N_y + 2, N_x + 2, m2_plus); double** m2_minus; createMatrix(N_y + 2, N_x + 2, m2_minus);
	double** rho_plus; createMatrix(N_y + 2, N_x + 2, rho_plus); double** rho_minus; createMatrix(N_y + 2, N_x + 2, rho_minus);
	F_d(m2_plus, xi, -1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(m2_minus, xi, 1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(rho_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5); F_d(rho_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			u_plus[i][j] = m2_plus[i][j] / rho_plus[i][j];
			u_minus[i][j] = m2_minus[i][j] / rho_minus[i][j];
		}
	}
	double** Cs_plus; createMatrix(N_y + 2, N_x + 2, Cs_plus); double** Cs_minus; createMatrix(N_y + 2, N_x + 2, Cs_minus);
	C_s(Cs_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	C_s(Cs_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	//alpha_c(N_y + 1, N_x)
	for (int j = 1; j < N_x + 1; j++)
	{
		for (int i = 1; i < N_y + 2; i++)
		{
			double temp_max = abs(u_plus[i][j] - Cs_plus[i][j]);
			if (abs(u_plus[i][j]) > temp_max)
			{
				temp_max = abs(u_plus[i][j]);
			}
			if (abs(u_plus[i][j] + Cs_plus[i][j]) > temp_max)
			{
				temp_max = abs(u_plus[i][j] + Cs_plus[i][j]);
			}
			if (abs(u_minus[i - 1][j] - Cs_minus[i - 1][j]) > temp_max)
			{
				temp_max = abs(u_minus[i - 1][j] - Cs_minus[i - 1][j]);
			}
			if (abs(u_minus[i - 1][j]) > temp_max)
			{
				temp_max = abs(u_minus[i - 1][j]);
			}
			if (abs(u_minus[i - 1][j] + Cs_minus[i - 1][j]) > temp_max)
			{
				temp_max = abs(u_minus[i - 1][j] + Cs_minus[i - 1][j]);
			}
			alpha_c[i - 1][j - 1] = temp_max;
		}
	}

	deleteMatrix(N_y + 2, m2_plus); deleteMatrix(N_y + 2, m2_minus);
	deleteMatrix(N_y + 2, rho_plus); deleteMatrix(N_y + 2, rho_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 2, Cs_plus); deleteMatrix(N_y + 2, Cs_minus);
}

void F1_d_hat(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** alpha_c)
{
	double** F1_plus; createMatrix(N_y + 2, N_x + 2, F1_plus); double** F1_minus; createMatrix(N_y + 2, N_x + 2, F1_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F_d(F1_plus, -1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_minus, 1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);

	F_d(u_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(u_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);

	double** F_ave; createMatrix(N_y, N_x + 3, F_ave); double** u_jump; createMatrix(N_y, N_x + 3, u_jump);
	compute_ave_x(F_ave, F1_plus, F1_minus); compute_jump_x(u_jump, u_plus, u_minus);
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			u[i][j] = F_ave[i][j + 1] - 0.5 * alpha_c[i][j] * u_jump[i][j + 1];
		}
	}
	deleteMatrix(N_y + 2, F1_plus); deleteMatrix(N_y + 2, F1_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y, F_ave); deleteMatrix(N_y, u_jump);
}

void F1_m1_hat(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** alpha_c)
{
	double** F1_plus; createMatrix(N_y + 2, N_x + 2, F1_plus); double** F1_minus; createMatrix(N_y + 2, N_x + 2, F1_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F1_m1_fun(F1_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_d(u_plus, -1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(u_minus, 1.0, eta, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);

	double** F_ave; createMatrix(N_y, N_x + 3, F_ave); double** u_jump; createMatrix(N_y, N_x + 3, u_jump);
	compute_ave_x(F_ave, F1_plus, F1_minus); compute_jump_x(u_jump, u_plus, u_minus);
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			u[i][j] = F_ave[i][j + 1] - 0.5 * alpha_c[i][j] * u_jump[i][j + 1];
		}
	}
	deleteMatrix(N_y + 2, F1_plus); deleteMatrix(N_y + 2, F1_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y, F_ave); deleteMatrix(N_y, u_jump);
}

void F1_m2_hat(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** alpha_c)
{
	double** F1_plus; createMatrix(N_y + 2, N_x + 2, F1_plus); double** F1_minus; createMatrix(N_y + 2, N_x + 2, F1_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F_m(F1_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_d(u_plus, -1.0, eta, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(u_minus, 1.0, eta, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	double** F_ave; createMatrix(N_y, N_x + 3, F_ave); double** u_jump; createMatrix(N_y, N_x + 3, u_jump);
	compute_ave_x(F_ave, F1_plus, F1_minus); compute_jump_x(u_jump, u_plus, u_minus);
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			u[i][j] = F_ave[i][j + 1] - 0.5 * alpha_c[i][j] * u_jump[i][j + 1];
		}
	}
	deleteMatrix(N_y + 2, F1_plus); deleteMatrix(N_y + 2, F1_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y, F_ave); deleteMatrix(N_y, u_jump);
}

void F1_alpha_hat(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double** alpha_c, double gamma_alpha)
{
	double** F1_plus; createMatrix(N_y + 2, N_x + 2, F1_plus); double** F1_minus; createMatrix(N_y + 2, N_x + 2, F1_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F1_alpha(F1_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F_d(u_plus, -1.0, eta, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	F_d(u_minus, 1.0, eta, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	double** F_ave; createMatrix(N_y, N_x + 3, F_ave); double** u_jump; createMatrix(N_y, N_x + 3, u_jump);
	compute_ave_x(F_ave, F1_plus, F1_minus); compute_jump_x(u_jump, u_plus, u_minus);
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			u[i][j] = F_ave[i][j + 1] - 0.5 * alpha_c[i][j] * u_jump[i][j + 1];
		}
	}
	deleteMatrix(N_y + 2, F1_plus); deleteMatrix(N_y + 2, F1_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y, F_ave); deleteMatrix(N_y, u_jump);
}

void F2_d_hat(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** alpha_c)
{
	double** F2_plus; createMatrix(N_y + 2, N_x + 2, F2_plus); double** F2_minus; createMatrix(N_y + 2, N_x + 2, F2_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F_d(F2_plus, xi, -1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_minus, xi, 1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_d(u_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(u_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);

	double** F_ave; createMatrix(N_y + 3, N_x, F_ave); double** u_jump; createMatrix(N_y + 3, N_x, u_jump);
	compute_ave_y(F_ave, F2_plus, F2_minus); compute_jump_y(u_jump, u_plus, u_minus);
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 1; i++)
		{
			u[i][j] = F_ave[i + 1][j] - 0.5 * alpha_c[i][j] * u_jump[i + 1][j];
		}
	}
	deleteMatrix(N_y + 2, F2_plus); deleteMatrix(N_y + 2, F2_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 3, F_ave); deleteMatrix(N_y + 3, u_jump);
}

void F2_m1_hat(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** alpha_c)
{
	double** F2_plus; createMatrix(N_y + 2, N_x + 2, F2_plus); double** F2_minus; createMatrix(N_y + 2, N_x + 2, F2_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F_m(F2_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_d(u_plus, xi, -1.0, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(u_minus, xi, 1.0, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);

	double** F_ave; createMatrix(N_y + 3, N_x, F_ave); double** u_jump; createMatrix(N_y + 3, N_x, u_jump);
	compute_ave_y(F_ave, F2_plus, F2_minus); compute_jump_y(u_jump, u_plus, u_minus);
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 1; i++)
		{
			u[i][j] = F_ave[i + 1][j] - 0.5 * alpha_c[i][j] * u_jump[i + 1][j];
		}
	}
	deleteMatrix(N_y + 2, F2_plus); deleteMatrix(N_y + 2, F2_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 3, F_ave); deleteMatrix(N_y + 3, u_jump);
}

void F2_m2_hat(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** alpha_c)
{
	double** F2_plus; createMatrix(N_y + 2, N_x + 2, F2_plus); double** F2_minus; createMatrix(N_y + 2, N_x + 2, F2_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F2_m2_fun(F2_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F_d(u_plus, xi, -1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(u_minus, xi, 1.0, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	double** F_ave; createMatrix(N_y + 3, N_x, F_ave); double** u_jump; createMatrix(N_y + 3, N_x, u_jump);
	compute_ave_y(F_ave, F2_plus, F2_minus); compute_jump_y(u_jump, u_plus, u_minus);
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 1; i++)
		{
			u[i][j] = F_ave[i + 1][j] - 0.5 * alpha_c[i][j] * u_jump[i + 1][j];
		}
	}
	deleteMatrix(N_y + 2, F2_plus); deleteMatrix(N_y + 2, F2_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 3, F_ave); deleteMatrix(N_y + 3, u_jump);
}

void F2_alpha_hat(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double** alpha_c, double gamma_alpha)
{
	double** F2_plus; createMatrix(N_y + 2, N_x + 2, F2_plus); double** F2_minus; createMatrix(N_y + 2, N_x + 2, F2_minus);
	double** u_plus; createMatrix(N_y + 2, N_x + 2, u_plus); double** u_minus; createMatrix(N_y + 2, N_x + 2, u_minus);

	F2_alpha(F2_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F_d(u_plus, xi, -1.0, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	F_d(u_minus, xi, 1.0, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	double** F_ave; createMatrix(N_y + 3, N_x, F_ave); double** u_jump; createMatrix(N_y + 3, N_x, u_jump);
	compute_ave_y(F_ave, F2_plus, F2_minus); compute_jump_y(u_jump, u_plus, u_minus);
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 1; i++)
		{
			u[i][j] = F_ave[i + 1][j] - 0.5 * alpha_c[i][j] * u_jump[i + 1][j];
		}
	}
	deleteMatrix(N_y + 2, F2_plus); deleteMatrix(N_y + 2, F2_minus);
	deleteMatrix(N_y + 2, u_plus); deleteMatrix(N_y + 2, u_minus);
	deleteMatrix(N_y + 3, F_ave); deleteMatrix(N_y + 3, u_jump);
}

void getGaussquadrature1D_1(double** u, double** F_1, double** F_2, double** F_3, double** F_4, double** F_5)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			u[i][j] = weight_1 * F_1[i][j] + weight_2 * F_2[i][j] + weight_3 * F_3[i][j] + weight_4 * F_4[i][j] + weight_5 * F_5[i][j];
		}
	}
}

void getGaussquadrature1D_p(double** u, double** F_1, double** F_2, double** F_3, double** F_4, double** F_5)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			u[i][j] = weight_1 * F_1[i][j] * GQxi_1 + weight_2 * F_2[i][j] * GQxi_2 + weight_3 * F_3[i][j] * GQxi_3 + weight_4 * F_4[i][j] * GQxi_4 + weight_5 * F_5[i][j] * GQxi_5;
		}
	}
}

void getGaussquadrature1D_p2(double** u, double** F_1, double** F_2, double** F_3, double** F_4, double** F_5)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			u[i][j] = weight_1 * F_1[i][j] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + weight_2 * F_2[i][j] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + weight_3 * F_3[i][j] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + weight_4 * F_4[i][j] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + weight_5 * F_5[i][j] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0;
		}
	}
}

void get_Gauss_quadrature_prime_2D_F1(double** U_0, double** U_1, double** U_2, double** u_11, double** u_21, double** u_31, double** u_41, double** u_51, double** u_12, double** u_22, double** u_32, double** u_42, double** u_52, double** u_13, double** u_23, double** u_33, double** u_43, double** u_53, double** u_14, double** u_24, double** u_34, double** u_44, double** u_54, double** u_15, double** u_25, double** u_35, double** u_45, double** u_55)
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

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i - 1][j - 1] = w_11 * u_11[i][j] + w_12 * u_12[i][j] + w_13 * u_13[i][j] + w_14 * u_14[i][j] + w_15 * u_15[i][j] + w_21 * u_21[i][j] + w_22 * u_22[i][j] + w_23 * u_23[i][j] + w_24 * u_24[i][j] + w_25 * u_25[i][j] + w_31 * u_31[i][j] + w_32 * u_32[i][j] + w_33 * u_33[i][j] + w_34 * u_34[i][j] + w_35 * u_35[i][j] + w_41 * u_41[i][j] + w_42 * u_42[i][j] + w_43 * u_43[i][j] + w_44 * u_44[i][j] + w_45 * u_45[i][j] + w_51 * u_51[i][j] + w_52 * u_52[i][j] + w_53 * u_53[i][j] + w_54 * u_54[i][j] + w_55 * u_55[i][j];

			U_1[i - 1][j - 1] = w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_1 + w_13 * u_13[i][j] * GQxi_1 + w_14 * u_14[i][j] * GQxi_1 + w_15 * u_15[i][j] * GQxi_1 + w_21 * u_21[i][j] * GQxi_2 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_2 + w_24 * u_24[i][j] * GQxi_2 + w_25 * u_25[i][j] * GQxi_2 + w_31 * u_31[i][j] * GQxi_3 + w_32 * u_32[i][j] * GQxi_3 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_3 + w_35 * u_35[i][j] * GQxi_3 + w_41 * u_41[i][j] * GQxi_4 + w_42 * u_42[i][j] * GQxi_4 + w_43 * u_43[i][j] * GQxi_4 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_4 + w_51 * u_51[i][j] * GQxi_5 + w_52 * u_52[i][j] * GQxi_5 + w_53 * u_53[i][j] * GQxi_5 + w_54 * u_54[i][j] * GQxi_5 + w_55 * u_55[i][j] * GQxi_5;

			U_2[i - 1][j - 1] = w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_2 + w_13 * u_13[i][j] * GQxi_3 + w_14 * u_14[i][j] * GQxi_4 + w_15 * u_15[i][j] * GQxi_5 + w_21 * u_21[i][j] * GQxi_1 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_3 + w_24 * u_24[i][j] * GQxi_4 + w_25 * u_25[i][j] * GQxi_5 + w_31 * u_31[i][j] * GQxi_1 + w_32 * u_32[i][j] * GQxi_2 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_4 + w_35 * u_35[i][j] * GQxi_5 + w_41 * u_41[i][j] * GQxi_1 + w_42 * u_42[i][j] * GQxi_2 + w_43 * u_43[i][j] * GQxi_3 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_5 + w_51 * u_51[i][j] * GQxi_1 + w_52 * u_52[i][j] * GQxi_2 + w_53 * u_53[i][j] * GQxi_3 + w_54 * u_54[i][j] * GQxi_4 + w_55 * u_55[i][j] * GQxi_5;
		}
	}
}

void get_Gauss_quadrature_prime_2D_F2(double** U_0, double** U_1, double** U_2, double** u_11, double** u_21, double** u_31, double** u_41, double** u_51, double** u_12, double** u_22, double** u_32, double** u_42, double** u_52, double** u_13, double** u_23, double** u_33, double** u_43, double** u_53, double** u_14, double** u_24, double** u_34, double** u_44, double** u_54, double** u_15, double** u_25, double** u_35, double** u_45, double** u_55)
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

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i - 1][j - 1] = w_11 * u_11[i][j] + w_12 * u_12[i][j] + w_13 * u_13[i][j] + w_14 * u_14[i][j] + w_15 * u_15[i][j] + w_21 * u_21[i][j] + w_22 * u_22[i][j] + w_23 * u_23[i][j] + w_24 * u_24[i][j] + w_25 * u_25[i][j] + w_31 * u_31[i][j] + w_32 * u_32[i][j] + w_33 * u_33[i][j] + w_34 * u_34[i][j] + w_35 * u_35[i][j] + w_41 * u_41[i][j] + w_42 * u_42[i][j] + w_43 * u_43[i][j] + w_44 * u_44[i][j] + w_45 * u_45[i][j] + w_51 * u_51[i][j] + w_52 * u_52[i][j] + w_53 * u_53[i][j] + w_54 * u_54[i][j] + w_55 * u_55[i][j];

			U_1[i - 1][j - 1] = w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_1 + w_13 * u_13[i][j] * GQxi_1 + w_14 * u_14[i][j] * GQxi_1 + w_15 * u_15[i][j] * GQxi_1 + w_21 * u_21[i][j] * GQxi_2 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_2 + w_24 * u_24[i][j] * GQxi_2 + w_25 * u_25[i][j] * GQxi_2 + w_31 * u_31[i][j] * GQxi_3 + w_32 * u_32[i][j] * GQxi_3 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_3 + w_35 * u_35[i][j] * GQxi_3 + w_41 * u_41[i][j] * GQxi_4 + w_42 * u_42[i][j] * GQxi_4 + w_43 * u_43[i][j] * GQxi_4 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_4 + w_51 * u_51[i][j] * GQxi_5 + w_52 * u_52[i][j] * GQxi_5 + w_53 * u_53[i][j] * GQxi_5 + w_54 * u_54[i][j] * GQxi_5 + w_55 * u_55[i][j] * GQxi_5;

			U_2[i - 1][j - 1] = w_11 * u_11[i][j] * GQxi_1 + w_12 * u_12[i][j] * GQxi_2 + w_13 * u_13[i][j] * GQxi_3 + w_14 * u_14[i][j] * GQxi_4 + w_15 * u_15[i][j] * GQxi_5 + w_21 * u_21[i][j] * GQxi_1 + w_22 * u_22[i][j] * GQxi_2 + w_23 * u_23[i][j] * GQxi_3 + w_24 * u_24[i][j] * GQxi_4 + w_25 * u_25[i][j] * GQxi_5 + w_31 * u_31[i][j] * GQxi_1 + w_32 * u_32[i][j] * GQxi_2 + w_33 * u_33[i][j] * GQxi_3 + w_34 * u_34[i][j] * GQxi_4 + w_35 * u_35[i][j] * GQxi_5 + w_41 * u_41[i][j] * GQxi_1 + w_42 * u_42[i][j] * GQxi_2 + w_43 * u_43[i][j] * GQxi_3 + w_44 * u_44[i][j] * GQxi_4 + w_45 * u_45[i][j] * GQxi_5 + w_51 * u_51[i][j] * GQxi_1 + w_52 * u_52[i][j] * GQxi_2 + w_53 * u_53[i][j] * GQxi_3 + w_54 * u_54[i][j] * GQxi_4 + w_55 * u_55[i][j] * GQxi_5;
		}
	}
}

void convection_d(double** con_0, double** con_1, double** con_2, double** con_3, double** con_4, double** con_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y)
{
	double** F1_hat1; createMatrix(N_y, N_x + 1, F1_hat1); double** F1_hat2; createMatrix(N_y, N_x + 1, F1_hat2);
	double** F1_hat3; createMatrix(N_y, N_x + 1, F1_hat3); double** F1_hat4; createMatrix(N_y, N_x + 1, F1_hat4);
	double** F1_hat5; createMatrix(N_y, N_x + 1, F1_hat5);

	double** F2_hat1; createMatrix(N_y + 1, N_x, F2_hat1); double** F2_hat2; createMatrix(N_y + 1, N_x, F2_hat2);
	double** F2_hat3; createMatrix(N_y + 1, N_x, F2_hat3); double** F2_hat4; createMatrix(N_y + 1, N_x, F2_hat4);
	double** F2_hat5; createMatrix(N_y + 1, N_x, F2_hat5);

	double** F1_11; createMatrix(N_y + 2, N_x + 2, F1_11); double** F1_12; createMatrix(N_y + 2, N_x + 2, F1_12); double** F1_13; createMatrix(N_y + 2, N_x + 2, F1_13);
	double** F1_14; createMatrix(N_y + 2, N_x + 2, F1_14); double** F1_15; createMatrix(N_y + 2, N_x + 2, F1_15);
	double** F1_21; createMatrix(N_y + 2, N_x + 2, F1_21); double** F1_22; createMatrix(N_y + 2, N_x + 2, F1_22); double** F1_23; createMatrix(N_y + 2, N_x + 2, F1_23);
	double** F1_24; createMatrix(N_y + 2, N_x + 2, F1_24); double** F1_25; createMatrix(N_y + 2, N_x + 2, F1_25);
	double** F1_31; createMatrix(N_y + 2, N_x + 2, F1_31); double** F1_32; createMatrix(N_y + 2, N_x + 2, F1_32); double** F1_33; createMatrix(N_y + 2, N_x + 2, F1_33);
	double** F1_34; createMatrix(N_y + 2, N_x + 2, F1_34); double** F1_35; createMatrix(N_y + 2, N_x + 2, F1_35);
	double** F1_41; createMatrix(N_y + 2, N_x + 2, F1_41); double** F1_42; createMatrix(N_y + 2, N_x + 2, F1_42); double** F1_43; createMatrix(N_y + 2, N_x + 2, F1_43);
	double** F1_44; createMatrix(N_y + 2, N_x + 2, F1_44); double** F1_45; createMatrix(N_y + 2, N_x + 2, F1_45);
	double** F1_51; createMatrix(N_y + 2, N_x + 2, F1_51); double** F1_52; createMatrix(N_y + 2, N_x + 2, F1_52); double** F1_53; createMatrix(N_y + 2, N_x + 2, F1_53);
	double** F1_54; createMatrix(N_y + 2, N_x + 2, F1_54); double** F1_55; createMatrix(N_y + 2, N_x + 2, F1_55);

	double** F2_11; createMatrix(N_y + 2, N_x + 2, F2_11); double** F2_12; createMatrix(N_y + 2, N_x + 2, F2_12); double** F2_13; createMatrix(N_y + 2, N_x + 2, F2_13);
	double** F2_14; createMatrix(N_y + 2, N_x + 2, F2_14); double** F2_15; createMatrix(N_y + 2, N_x + 2, F2_15);
	double** F2_21; createMatrix(N_y + 2, N_x + 2, F2_21); double** F2_22; createMatrix(N_y + 2, N_x + 2, F2_22); double** F2_23; createMatrix(N_y + 2, N_x + 2, F2_23);
	double** F2_24; createMatrix(N_y + 2, N_x + 2, F2_24); double** F2_25; createMatrix(N_y + 2, N_x + 2, F2_25);
	double** F2_31; createMatrix(N_y + 2, N_x + 2, F2_31); double** F2_32; createMatrix(N_y + 2, N_x + 2, F2_32); double** F2_33; createMatrix(N_y + 2, N_x + 2, F2_33);
	double** F2_34; createMatrix(N_y + 2, N_x + 2, F2_34); double** F2_35; createMatrix(N_y + 2, N_x + 2, F2_35);
	double** F2_41; createMatrix(N_y + 2, N_x + 2, F2_41); double** F2_42; createMatrix(N_y + 2, N_x + 2, F2_42); double** F2_43; createMatrix(N_y + 2, N_x + 2, F2_43);
	double** F2_44; createMatrix(N_y + 2, N_x + 2, F2_44); double** F2_45; createMatrix(N_y + 2, N_x + 2, F2_45);
	double** F2_51; createMatrix(N_y + 2, N_x + 2, F2_51); double** F2_52; createMatrix(N_y + 2, N_x + 2, F2_52); double** F2_53; createMatrix(N_y + 2, N_x + 2, F2_53);
	double** F2_54; createMatrix(N_y + 2, N_x + 2, F2_54); double** F2_55; createMatrix(N_y + 2, N_x + 2, F2_55);

	double** F1_m1; createMatrix(N_y, N_x, F1_m1); double** F1_m2; createMatrix(N_y, N_x, F1_m2); double** F1_m3; createMatrix(N_y, N_x, F1_m3);
	double** F1_m4; createMatrix(N_y, N_x, F1_m4); double** F1_m5; createMatrix(N_y, N_x, F1_m5);
	double** F1_p1; createMatrix(N_y, N_x, F1_p1); double** F1_p2; createMatrix(N_y, N_x, F1_p2); double** F1_p3; createMatrix(N_y, N_x, F1_p3);
	double** F1_p4; createMatrix(N_y, N_x, F1_p4); double** F1_p5; createMatrix(N_y, N_x, F1_p5);

	double** F2_m1; createMatrix(N_y, N_x, F2_m1); double** F2_m2; createMatrix(N_y, N_x, F2_m2); double** F2_m3; createMatrix(N_y, N_x, F2_m3);
	double** F2_m4; createMatrix(N_y, N_x, F2_m4); double** F2_m5; createMatrix(N_y, N_x, F2_m5);
	double** F2_p1; createMatrix(N_y, N_x, F2_p1); double** F2_p2; createMatrix(N_y, N_x, F2_p2); double** F2_p3; createMatrix(N_y, N_x, F2_p3);
	double** F2_p4; createMatrix(N_y, N_x, F2_p4); double** F2_p5; createMatrix(N_y, N_x, F2_p5);

	double** F1m_GQ_0; createMatrix(N_y, N_x, F1m_GQ_0); double** F1m_GQ_1; createMatrix(N_y, N_x, F1m_GQ_1);
	double** F1m_GQ_2; createMatrix(N_y, N_x, F1m_GQ_2);
	double** F1p_GQ_0; createMatrix(N_y, N_x, F1p_GQ_0); double** F1p_GQ_1; createMatrix(N_y, N_x, F1p_GQ_1);

	double** F2m_GQ_0; createMatrix(N_y, N_x, F2m_GQ_0); double** F2m_GQ_1; createMatrix(N_y, N_x, F2m_GQ_1);
	double** F2m_GQ_2; createMatrix(N_y, N_x, F2m_GQ_2);
	double** F2p_GQ_0; createMatrix(N_y, N_x, F2p_GQ_0); double** F2p_GQ_1; createMatrix(N_y, N_x, F2p_GQ_1);

	double** F1_GQ_0; createMatrix(N_y, N_x, F1_GQ_0); double** F1_GQ_1; createMatrix(N_y, N_x, F1_GQ_1);
	double** F1_GQ_2; createMatrix(N_y, N_x, F1_GQ_2);

	double** F2_GQ_0; createMatrix(N_y, N_x, F2_GQ_0); double** F2_GQ_1; createMatrix(N_y, N_x, F2_GQ_1);
	double** F2_GQ_2; createMatrix(N_y, N_x, F2_GQ_2);

	F1_d_hat(F1_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, F1_alphac);
	F1_d_hat(F1_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, F1_alphac);
	F1_d_hat(F1_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, F1_alphac);
	F1_d_hat(F1_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, F1_alphac);
	F1_d_hat(F1_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, F1_alphac);

	F2_d_hat(F2_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_d_hat(F2_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_d_hat(F2_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_d_hat(F2_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_d_hat(F2_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);


	F_d(F1_11, GQxi_1, GQxi_1, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_12, GQxi_1, GQxi_2, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_13, GQxi_1, GQxi_3, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_14, GQxi_1, GQxi_4, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_15, GQxi_1, GQxi_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_21, GQxi_2, GQxi_1, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_22, GQxi_2, GQxi_2, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_23, GQxi_2, GQxi_3, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_24, GQxi_2, GQxi_4, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_25, GQxi_2, GQxi_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_31, GQxi_3, GQxi_1, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_32, GQxi_3, GQxi_2, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_33, GQxi_3, GQxi_3, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_34, GQxi_3, GQxi_4, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_35, GQxi_3, GQxi_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_41, GQxi_4, GQxi_1, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_42, GQxi_4, GQxi_2, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_43, GQxi_4, GQxi_3, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_44, GQxi_4, GQxi_4, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_45, GQxi_4, GQxi_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_51, GQxi_5, GQxi_1, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_52, GQxi_5, GQxi_2, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_53, GQxi_5, GQxi_3, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5); F_d(F1_54, GQxi_5, GQxi_4, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	F_d(F1_55, GQxi_5, GQxi_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);


	F_d(F2_11, GQxi_1, GQxi_1, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_12, GQxi_1, GQxi_2, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_13, GQxi_1, GQxi_3, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_14, GQxi_1, GQxi_4, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_15, GQxi_1, GQxi_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_21, GQxi_2, GQxi_1, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_22, GQxi_2, GQxi_2, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_23, GQxi_2, GQxi_3, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_24, GQxi_2, GQxi_4, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_25, GQxi_2, GQxi_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_31, GQxi_3, GQxi_1, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_32, GQxi_3, GQxi_2, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_33, GQxi_3, GQxi_3, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_34, GQxi_3, GQxi_4, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_35, GQxi_3, GQxi_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_41, GQxi_4, GQxi_1, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_42, GQxi_4, GQxi_2, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_43, GQxi_4, GQxi_3, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_44, GQxi_4, GQxi_4, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_45, GQxi_4, GQxi_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_51, GQxi_5, GQxi_1, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_52, GQxi_5, GQxi_2, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_53, GQxi_5, GQxi_3, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5); F_d(F2_54, GQxi_5, GQxi_4, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_d(F2_55, GQxi_5, GQxi_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			F1_m1[i][j] = F1_hat1[i][j + 1] - F1_hat1[i][j]; F1_m2[i][j] = F1_hat2[i][j + 1] - F1_hat2[i][j];
			F1_m3[i][j] = F1_hat3[i][j + 1] - F1_hat3[i][j]; F1_m4[i][j] = F1_hat4[i][j + 1] - F1_hat4[i][j];
			F1_m5[i][j] = F1_hat5[i][j + 1] - F1_hat5[i][j];

			F1_p1[i][j] = F1_hat1[i][j + 1] + F1_hat1[i][j]; F1_p2[i][j] = F1_hat2[i][j + 1] + F1_hat2[i][j];
			F1_p3[i][j] = F1_hat3[i][j + 1] + F1_hat3[i][j]; F1_p4[i][j] = F1_hat4[i][j + 1] + F1_hat4[i][j];
			F1_p5[i][j] = F1_hat5[i][j + 1] + F1_hat5[i][j];
		}
	}

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y; i++)
		{
			F2_m1[i][j] = F2_hat1[i + 1][j] - F2_hat1[i][j]; F2_m2[i][j] = F2_hat2[i + 1][j] - F2_hat2[i][j];
			F2_m3[i][j] = F2_hat3[i + 1][j] - F2_hat3[i][j]; F2_m4[i][j] = F2_hat4[i + 1][j] - F2_hat4[i][j];
			F2_m5[i][j] = F2_hat5[i + 1][j] - F2_hat5[i][j];

			F2_p1[i][j] = F2_hat1[i + 1][j] + F2_hat1[i][j]; F2_p2[i][j] = F2_hat2[i + 1][j] + F2_hat2[i][j];
			F2_p3[i][j] = F2_hat3[i + 1][j] + F2_hat3[i][j]; F2_p4[i][j] = F2_hat4[i + 1][j] + F2_hat4[i][j];
			F2_p5[i][j] = F2_hat5[i + 1][j] + F2_hat5[i][j];
		}
	}

	getGaussquadrature1D_1(F1m_GQ_0, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5); getGaussquadrature1D_p(F1m_GQ_1, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_p2(F1m_GQ_2, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_1(F1p_GQ_0, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5); getGaussquadrature1D_p(F1p_GQ_1, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5);

	getGaussquadrature1D_1(F2m_GQ_0, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5); getGaussquadrature1D_p(F2m_GQ_1, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_p2(F2m_GQ_2, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_1(F2p_GQ_0, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5); getGaussquadrature1D_p(F2p_GQ_1, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5);


	get_Gauss_quadrature_prime_2D_F1(F1_GQ_0, F1_GQ_1, F1_GQ_2, F1_11, F1_21, F1_31, F1_41, F1_51, F1_12, F1_22, F1_32, F1_42, F1_52, F1_13, F1_23, F1_33, F1_43, F1_53, F1_14, F1_24, F1_34, F1_44, F1_54, F1_15, F1_25, F1_35, F1_45, F1_55);
	get_Gauss_quadrature_prime_2D_F2(F2_GQ_0, F2_GQ_1, F2_GQ_2, F2_11, F2_21, F2_31, F2_41, F2_51, F2_12, F2_22, F2_32, F2_42, F2_52, F2_13, F2_23, F2_33, F2_43, F2_53, F2_14, F2_24, F2_34, F2_44, F2_54, F2_15, F2_25, F2_35, F2_45, F2_55);


	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			con_0[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_0[i][j];
			con_1[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_1[i][j] + h_y / 2.0 * F1_GQ_0[i][j];
			con_2[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_0[i][j] + h_x / 2.0 * F2_GQ_0[i][j];
			con_3[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_1[i][j] + h_y / 2.0 * F1_GQ_2[i][j] + h_x / 2.0 * F2_GQ_1[i][j];
			con_4[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_2[i][j] + 3.0 * h_y / 2.0 * F1_GQ_1[i][j];
			con_5[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_2[i][j] - h_x / 2.0 * F2m_GQ_0[i][j] + 3.0 * h_x / 2.0 * F2_GQ_2[i][j];
		}
	}

	deleteMatrix(N_y, F1_hat1); deleteMatrix(N_y, F1_hat2); deleteMatrix(N_y, F1_hat3);
	deleteMatrix(N_y, F1_hat4); deleteMatrix(N_y, F1_hat5);

	deleteMatrix(N_y + 1, F2_hat1); deleteMatrix(N_y + 1, F2_hat2); deleteMatrix(N_y + 1, F2_hat3);
	deleteMatrix(N_y + 1, F2_hat4); deleteMatrix(N_y + 1, F2_hat5);

	deleteMatrix(N_y + 2, F1_11); deleteMatrix(N_y + 2, F1_12); deleteMatrix(N_y + 2, F1_13); deleteMatrix(N_y + 2, F1_14); deleteMatrix(N_y + 2, F1_15);
	deleteMatrix(N_y + 2, F1_21); deleteMatrix(N_y + 2, F1_22); deleteMatrix(N_y + 2, F1_23); deleteMatrix(N_y + 2, F1_24); deleteMatrix(N_y + 2, F1_25);
	deleteMatrix(N_y + 2, F1_31); deleteMatrix(N_y + 2, F1_32); deleteMatrix(N_y + 2, F1_33); deleteMatrix(N_y + 2, F1_34); deleteMatrix(N_y + 2, F1_35);
	deleteMatrix(N_y + 2, F1_41); deleteMatrix(N_y + 2, F1_42); deleteMatrix(N_y + 2, F1_43); deleteMatrix(N_y + 2, F1_44); deleteMatrix(N_y + 2, F1_45);
	deleteMatrix(N_y + 2, F1_51); deleteMatrix(N_y + 2, F1_52); deleteMatrix(N_y + 2, F1_53); deleteMatrix(N_y + 2, F1_54); deleteMatrix(N_y + 2, F1_55);

	deleteMatrix(N_y + 2, F2_11); deleteMatrix(N_y + 2, F2_12); deleteMatrix(N_y + 2, F2_13); deleteMatrix(N_y + 2, F2_14); deleteMatrix(N_y + 2, F2_15);
	deleteMatrix(N_y + 2, F2_21); deleteMatrix(N_y + 2, F2_22); deleteMatrix(N_y + 2, F2_23); deleteMatrix(N_y + 2, F2_24); deleteMatrix(N_y + 2, F2_25);
	deleteMatrix(N_y + 2, F2_31); deleteMatrix(N_y + 2, F2_32); deleteMatrix(N_y + 2, F2_33); deleteMatrix(N_y + 2, F2_34); deleteMatrix(N_y + 2, F2_35);
	deleteMatrix(N_y + 2, F2_41); deleteMatrix(N_y + 2, F2_42); deleteMatrix(N_y + 2, F2_43); deleteMatrix(N_y + 2, F2_44); deleteMatrix(N_y + 2, F2_45);
	deleteMatrix(N_y + 2, F2_51); deleteMatrix(N_y + 2, F2_52); deleteMatrix(N_y + 2, F2_53); deleteMatrix(N_y + 2, F2_54); deleteMatrix(N_y + 2, F2_55);

	deleteMatrix(N_y, F1_m1); deleteMatrix(N_y, F1_m2); deleteMatrix(N_y, F1_m3); deleteMatrix(N_y, F1_m4); deleteMatrix(N_y, F1_m5);
	deleteMatrix(N_y, F1_p1); deleteMatrix(N_y, F1_p2); deleteMatrix(N_y, F1_p3); deleteMatrix(N_y, F1_p4); deleteMatrix(N_y, F1_p5);

	deleteMatrix(N_y, F2_m1); deleteMatrix(N_y, F2_m2); deleteMatrix(N_y, F2_m3); deleteMatrix(N_y, F2_m4); deleteMatrix(N_y, F2_m5);
	deleteMatrix(N_y, F2_p1); deleteMatrix(N_y, F2_p2); deleteMatrix(N_y, F2_p3); deleteMatrix(N_y, F2_p4); deleteMatrix(N_y, F2_p5);

	deleteMatrix(N_y, F1m_GQ_0); deleteMatrix(N_y, F1m_GQ_1); deleteMatrix(N_y, F1m_GQ_2);
	deleteMatrix(N_y, F1p_GQ_0); deleteMatrix(N_y, F1p_GQ_1);

	deleteMatrix(N_y, F2m_GQ_0); deleteMatrix(N_y, F2m_GQ_1); deleteMatrix(N_y, F2m_GQ_2);
	deleteMatrix(N_y, F2p_GQ_0); deleteMatrix(N_y, F2p_GQ_1);

	deleteMatrix(N_y, F1_GQ_0); deleteMatrix(N_y, F1_GQ_1); deleteMatrix(N_y, F1_GQ_2);
	deleteMatrix(N_y, F2_GQ_0); deleteMatrix(N_y, F2_GQ_1); deleteMatrix(N_y, F2_GQ_2);
}

void convection_m1(double** con_0, double** con_1, double** con_2, double** con_3, double** con_4, double** con_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y)
{
	double** F1_hat1; createMatrix(N_y, N_x + 1, F1_hat1); double** F1_hat2; createMatrix(N_y, N_x + 1, F1_hat2);
	double** F1_hat3; createMatrix(N_y, N_x + 1, F1_hat3); double** F1_hat4; createMatrix(N_y, N_x + 1, F1_hat4);
	double** F1_hat5; createMatrix(N_y, N_x + 1, F1_hat5);

	double** F2_hat1; createMatrix(N_y + 1, N_x, F2_hat1); double** F2_hat2; createMatrix(N_y + 1, N_x, F2_hat2);
	double** F2_hat3; createMatrix(N_y + 1, N_x, F2_hat3); double** F2_hat4; createMatrix(N_y + 1, N_x, F2_hat4);
	double** F2_hat5; createMatrix(N_y + 1, N_x, F2_hat5);

	double** F1_11; createMatrix(N_y + 2, N_x + 2, F1_11); double** F1_12; createMatrix(N_y + 2, N_x + 2, F1_12); double** F1_13; createMatrix(N_y + 2, N_x + 2, F1_13);
	double** F1_14; createMatrix(N_y + 2, N_x + 2, F1_14); double** F1_15; createMatrix(N_y + 2, N_x + 2, F1_15);
	double** F1_21; createMatrix(N_y + 2, N_x + 2, F1_21); double** F1_22; createMatrix(N_y + 2, N_x + 2, F1_22); double** F1_23; createMatrix(N_y + 2, N_x + 2, F1_23);
	double** F1_24; createMatrix(N_y + 2, N_x + 2, F1_24); double** F1_25; createMatrix(N_y + 2, N_x + 2, F1_25);
	double** F1_31; createMatrix(N_y + 2, N_x + 2, F1_31); double** F1_32; createMatrix(N_y + 2, N_x + 2, F1_32); double** F1_33; createMatrix(N_y + 2, N_x + 2, F1_33);
	double** F1_34; createMatrix(N_y + 2, N_x + 2, F1_34); double** F1_35; createMatrix(N_y + 2, N_x + 2, F1_35);
	double** F1_41; createMatrix(N_y + 2, N_x + 2, F1_41); double** F1_42; createMatrix(N_y + 2, N_x + 2, F1_42); double** F1_43; createMatrix(N_y + 2, N_x + 2, F1_43);
	double** F1_44; createMatrix(N_y + 2, N_x + 2, F1_44); double** F1_45; createMatrix(N_y + 2, N_x + 2, F1_45);
	double** F1_51; createMatrix(N_y + 2, N_x + 2, F1_51); double** F1_52; createMatrix(N_y + 2, N_x + 2, F1_52); double** F1_53; createMatrix(N_y + 2, N_x + 2, F1_53);
	double** F1_54; createMatrix(N_y + 2, N_x + 2, F1_54); double** F1_55; createMatrix(N_y + 2, N_x + 2, F1_55);

	double** F2_11; createMatrix(N_y + 2, N_x + 2, F2_11); double** F2_12; createMatrix(N_y + 2, N_x + 2, F2_12); double** F2_13; createMatrix(N_y + 2, N_x + 2, F2_13);
	double** F2_14; createMatrix(N_y + 2, N_x + 2, F2_14); double** F2_15; createMatrix(N_y + 2, N_x + 2, F2_15);
	double** F2_21; createMatrix(N_y + 2, N_x + 2, F2_21); double** F2_22; createMatrix(N_y + 2, N_x + 2, F2_22); double** F2_23; createMatrix(N_y + 2, N_x + 2, F2_23);
	double** F2_24; createMatrix(N_y + 2, N_x + 2, F2_24); double** F2_25; createMatrix(N_y + 2, N_x + 2, F2_25);
	double** F2_31; createMatrix(N_y + 2, N_x + 2, F2_31); double** F2_32; createMatrix(N_y + 2, N_x + 2, F2_32); double** F2_33; createMatrix(N_y + 2, N_x + 2, F2_33);
	double** F2_34; createMatrix(N_y + 2, N_x + 2, F2_34); double** F2_35; createMatrix(N_y + 2, N_x + 2, F2_35);
	double** F2_41; createMatrix(N_y + 2, N_x + 2, F2_41); double** F2_42; createMatrix(N_y + 2, N_x + 2, F2_42); double** F2_43; createMatrix(N_y + 2, N_x + 2, F2_43);
	double** F2_44; createMatrix(N_y + 2, N_x + 2, F2_44); double** F2_45; createMatrix(N_y + 2, N_x + 2, F2_45);
	double** F2_51; createMatrix(N_y + 2, N_x + 2, F2_51); double** F2_52; createMatrix(N_y + 2, N_x + 2, F2_52); double** F2_53; createMatrix(N_y + 2, N_x + 2, F2_53);
	double** F2_54; createMatrix(N_y + 2, N_x + 2, F2_54); double** F2_55; createMatrix(N_y + 2, N_x + 2, F2_55);

	double** F1_m1; createMatrix(N_y, N_x, F1_m1); double** F1_m2; createMatrix(N_y, N_x, F1_m2); double** F1_m3; createMatrix(N_y, N_x, F1_m3);
	double** F1_m4; createMatrix(N_y, N_x, F1_m4); double** F1_m5; createMatrix(N_y, N_x, F1_m5);
	double** F1_p1; createMatrix(N_y, N_x, F1_p1); double** F1_p2; createMatrix(N_y, N_x, F1_p2); double** F1_p3; createMatrix(N_y, N_x, F1_p3);
	double** F1_p4; createMatrix(N_y, N_x, F1_p4); double** F1_p5; createMatrix(N_y, N_x, F1_p5);

	double** F2_m1; createMatrix(N_y, N_x, F2_m1); double** F2_m2; createMatrix(N_y, N_x, F2_m2); double** F2_m3; createMatrix(N_y, N_x, F2_m3);
	double** F2_m4; createMatrix(N_y, N_x, F2_m4); double** F2_m5; createMatrix(N_y, N_x, F2_m5);
	double** F2_p1; createMatrix(N_y, N_x, F2_p1); double** F2_p2; createMatrix(N_y, N_x, F2_p2); double** F2_p3; createMatrix(N_y, N_x, F2_p3);
	double** F2_p4; createMatrix(N_y, N_x, F2_p4); double** F2_p5; createMatrix(N_y, N_x, F2_p5);

	double** F1m_GQ_0; createMatrix(N_y, N_x, F1m_GQ_0); double** F1m_GQ_1; createMatrix(N_y, N_x, F1m_GQ_1);
	double** F1m_GQ_2; createMatrix(N_y, N_x, F1m_GQ_2);
	double** F1p_GQ_0; createMatrix(N_y, N_x, F1p_GQ_0); double** F1p_GQ_1; createMatrix(N_y, N_x, F1p_GQ_1);

	double** F2m_GQ_0; createMatrix(N_y, N_x, F2m_GQ_0); double** F2m_GQ_1; createMatrix(N_y, N_x, F2m_GQ_1);
	double** F2m_GQ_2; createMatrix(N_y, N_x, F2m_GQ_2);
	double** F2p_GQ_0; createMatrix(N_y, N_x, F2p_GQ_0); double** F2p_GQ_1; createMatrix(N_y, N_x, F2p_GQ_1);

	double** F1_GQ_0; createMatrix(N_y, N_x, F1_GQ_0); double** F1_GQ_1; createMatrix(N_y, N_x, F1_GQ_1);
	double** F1_GQ_2; createMatrix(N_y, N_x, F1_GQ_2);

	double** F2_GQ_0; createMatrix(N_y, N_x, F2_GQ_0); double** F2_GQ_1; createMatrix(N_y, N_x, F2_GQ_1);
	double** F2_GQ_2; createMatrix(N_y, N_x, F2_GQ_2);

	F1_m1_hat(F1_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac);
	F1_m1_hat(F1_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac);
	F1_m1_hat(F1_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac);
	F1_m1_hat(F1_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac);
	F1_m1_hat(F1_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac);

	F2_m1_hat(F2_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_m1_hat(F2_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_m1_hat(F2_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_m1_hat(F2_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);
	F2_m1_hat(F2_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F2_alphac);



	F1_m1_fun(F1_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F1_m1_fun(F1_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F1_m1_fun(F1_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F1_m1_fun(F1_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F1_m1_fun(F1_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F1_m1_fun(F1_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);


	F_m(F2_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F2_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F2_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F2_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F2_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F2_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			F1_m1[i][j] = F1_hat1[i][j + 1] - F1_hat1[i][j]; F1_m2[i][j] = F1_hat2[i][j + 1] - F1_hat2[i][j];
			F1_m3[i][j] = F1_hat3[i][j + 1] - F1_hat3[i][j]; F1_m4[i][j] = F1_hat4[i][j + 1] - F1_hat4[i][j];
			F1_m5[i][j] = F1_hat5[i][j + 1] - F1_hat5[i][j];

			F1_p1[i][j] = F1_hat1[i][j + 1] + F1_hat1[i][j]; F1_p2[i][j] = F1_hat2[i][j + 1] + F1_hat2[i][j];
			F1_p3[i][j] = F1_hat3[i][j + 1] + F1_hat3[i][j]; F1_p4[i][j] = F1_hat4[i][j + 1] + F1_hat4[i][j];
			F1_p5[i][j] = F1_hat5[i][j + 1] + F1_hat5[i][j];
		}
	}

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y; i++)
		{
			F2_m1[i][j] = F2_hat1[i + 1][j] - F2_hat1[i][j]; F2_m2[i][j] = F2_hat2[i + 1][j] - F2_hat2[i][j];
			F2_m3[i][j] = F2_hat3[i + 1][j] - F2_hat3[i][j]; F2_m4[i][j] = F2_hat4[i + 1][j] - F2_hat4[i][j];
			F2_m5[i][j] = F2_hat5[i + 1][j] - F2_hat5[i][j];

			F2_p1[i][j] = F2_hat1[i + 1][j] + F2_hat1[i][j]; F2_p2[i][j] = F2_hat2[i + 1][j] + F2_hat2[i][j];
			F2_p3[i][j] = F2_hat3[i + 1][j] + F2_hat3[i][j]; F2_p4[i][j] = F2_hat4[i + 1][j] + F2_hat4[i][j];
			F2_p5[i][j] = F2_hat5[i + 1][j] + F2_hat5[i][j];
		}
	}

	getGaussquadrature1D_1(F1m_GQ_0, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5); getGaussquadrature1D_p(F1m_GQ_1, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_p2(F1m_GQ_2, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_1(F1p_GQ_0, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5); getGaussquadrature1D_p(F1p_GQ_1, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5);

	getGaussquadrature1D_1(F2m_GQ_0, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5); getGaussquadrature1D_p(F2m_GQ_1, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_p2(F2m_GQ_2, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_1(F2p_GQ_0, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5); getGaussquadrature1D_p(F2p_GQ_1, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5);


	get_Gauss_quadrature_prime_2D_F1(F1_GQ_0, F1_GQ_1, F1_GQ_2, F1_11, F1_21, F1_31, F1_41, F1_51, F1_12, F1_22, F1_32, F1_42, F1_52, F1_13, F1_23, F1_33, F1_43, F1_53, F1_14, F1_24, F1_34, F1_44, F1_54, F1_15, F1_25, F1_35, F1_45, F1_55);
	get_Gauss_quadrature_prime_2D_F2(F2_GQ_0, F2_GQ_1, F2_GQ_2, F2_11, F2_21, F2_31, F2_41, F2_51, F2_12, F2_22, F2_32, F2_42, F2_52, F2_13, F2_23, F2_33, F2_43, F2_53, F2_14, F2_24, F2_34, F2_44, F2_54, F2_15, F2_25, F2_35, F2_45, F2_55);


	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			con_0[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_0[i][j];
			con_1[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_1[i][j] + h_y / 2.0 * F1_GQ_0[i][j];
			con_2[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_0[i][j] + h_x / 2.0 * F2_GQ_0[i][j];
			con_3[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_1[i][j] + h_y / 2.0 * F1_GQ_2[i][j] + h_x / 2.0 * F2_GQ_1[i][j];
			con_4[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_2[i][j] + 3.0 * h_y / 2.0 * F1_GQ_1[i][j];
			con_5[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_2[i][j] - h_x / 2.0 * F2m_GQ_0[i][j] + 3.0 * h_x / 2.0 * F2_GQ_2[i][j];
		}
	}

	deleteMatrix(N_y, F1_hat1); deleteMatrix(N_y, F1_hat2); deleteMatrix(N_y, F1_hat3);
	deleteMatrix(N_y, F1_hat4); deleteMatrix(N_y, F1_hat5);

	deleteMatrix(N_y + 1, F2_hat1); deleteMatrix(N_y + 1, F2_hat2); deleteMatrix(N_y + 1, F2_hat3);
	deleteMatrix(N_y + 1, F2_hat4); deleteMatrix(N_y + 1, F2_hat5);

	deleteMatrix(N_y + 2, F1_11); deleteMatrix(N_y + 2, F1_12); deleteMatrix(N_y + 2, F1_13); deleteMatrix(N_y + 2, F1_14); deleteMatrix(N_y + 2, F1_15);
	deleteMatrix(N_y + 2, F1_21); deleteMatrix(N_y + 2, F1_22); deleteMatrix(N_y + 2, F1_23); deleteMatrix(N_y + 2, F1_24); deleteMatrix(N_y + 2, F1_25);
	deleteMatrix(N_y + 2, F1_31); deleteMatrix(N_y + 2, F1_32); deleteMatrix(N_y + 2, F1_33); deleteMatrix(N_y + 2, F1_34); deleteMatrix(N_y + 2, F1_35);
	deleteMatrix(N_y + 2, F1_41); deleteMatrix(N_y + 2, F1_42); deleteMatrix(N_y + 2, F1_43); deleteMatrix(N_y + 2, F1_44); deleteMatrix(N_y + 2, F1_45);
	deleteMatrix(N_y + 2, F1_51); deleteMatrix(N_y + 2, F1_52); deleteMatrix(N_y + 2, F1_53); deleteMatrix(N_y + 2, F1_54); deleteMatrix(N_y + 2, F1_55);

	deleteMatrix(N_y + 2, F2_11); deleteMatrix(N_y + 2, F2_12); deleteMatrix(N_y + 2, F2_13); deleteMatrix(N_y + 2, F2_14); deleteMatrix(N_y + 2, F2_15);
	deleteMatrix(N_y + 2, F2_21); deleteMatrix(N_y + 2, F2_22); deleteMatrix(N_y + 2, F2_23); deleteMatrix(N_y + 2, F2_24); deleteMatrix(N_y + 2, F2_25);
	deleteMatrix(N_y + 2, F2_31); deleteMatrix(N_y + 2, F2_32); deleteMatrix(N_y + 2, F2_33); deleteMatrix(N_y + 2, F2_34); deleteMatrix(N_y + 2, F2_35);
	deleteMatrix(N_y + 2, F2_41); deleteMatrix(N_y + 2, F2_42); deleteMatrix(N_y + 2, F2_43); deleteMatrix(N_y + 2, F2_44); deleteMatrix(N_y + 2, F2_45);
	deleteMatrix(N_y + 2, F2_51); deleteMatrix(N_y + 2, F2_52); deleteMatrix(N_y + 2, F2_53); deleteMatrix(N_y + 2, F2_54); deleteMatrix(N_y + 2, F2_55);

	deleteMatrix(N_y, F1_m1); deleteMatrix(N_y, F1_m2); deleteMatrix(N_y, F1_m3); deleteMatrix(N_y, F1_m4); deleteMatrix(N_y, F1_m5);
	deleteMatrix(N_y, F1_p1); deleteMatrix(N_y, F1_p2); deleteMatrix(N_y, F1_p3); deleteMatrix(N_y, F1_p4); deleteMatrix(N_y, F1_p5);

	deleteMatrix(N_y, F2_m1); deleteMatrix(N_y, F2_m2); deleteMatrix(N_y, F2_m3); deleteMatrix(N_y, F2_m4); deleteMatrix(N_y, F2_m5);
	deleteMatrix(N_y, F2_p1); deleteMatrix(N_y, F2_p2); deleteMatrix(N_y, F2_p3); deleteMatrix(N_y, F2_p4); deleteMatrix(N_y, F2_p5);

	deleteMatrix(N_y, F1m_GQ_0); deleteMatrix(N_y, F1m_GQ_1); deleteMatrix(N_y, F1m_GQ_2);
	deleteMatrix(N_y, F1p_GQ_0); deleteMatrix(N_y, F1p_GQ_1);

	deleteMatrix(N_y, F2m_GQ_0); deleteMatrix(N_y, F2m_GQ_1); deleteMatrix(N_y, F2m_GQ_2);
	deleteMatrix(N_y, F2p_GQ_0); deleteMatrix(N_y, F2p_GQ_1);

	deleteMatrix(N_y, F1_GQ_0); deleteMatrix(N_y, F1_GQ_1); deleteMatrix(N_y, F1_GQ_2);
	deleteMatrix(N_y, F2_GQ_0); deleteMatrix(N_y, F2_GQ_1); deleteMatrix(N_y, F2_GQ_2);
}

void convection_m2(double** con_0, double** con_1, double** con_2, double** con_3, double** con_4, double** con_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y)
{
	double** F1_hat1; createMatrix(N_y, N_x + 1, F1_hat1); double** F1_hat2; createMatrix(N_y, N_x + 1, F1_hat2);
	double** F1_hat3; createMatrix(N_y, N_x + 1, F1_hat3); double** F1_hat4; createMatrix(N_y, N_x + 1, F1_hat4);
	double** F1_hat5; createMatrix(N_y, N_x + 1, F1_hat5);

	double** F2_hat1; createMatrix(N_y + 1, N_x, F2_hat1); double** F2_hat2; createMatrix(N_y + 1, N_x, F2_hat2);
	double** F2_hat3; createMatrix(N_y + 1, N_x, F2_hat3); double** F2_hat4; createMatrix(N_y + 1, N_x, F2_hat4);
	double** F2_hat5; createMatrix(N_y + 1, N_x, F2_hat5);

	double** F1_11; createMatrix(N_y + 2, N_x + 2, F1_11); double** F1_12; createMatrix(N_y + 2, N_x + 2, F1_12); double** F1_13; createMatrix(N_y + 2, N_x + 2, F1_13);
	double** F1_14; createMatrix(N_y + 2, N_x + 2, F1_14); double** F1_15; createMatrix(N_y + 2, N_x + 2, F1_15);
	double** F1_21; createMatrix(N_y + 2, N_x + 2, F1_21); double** F1_22; createMatrix(N_y + 2, N_x + 2, F1_22); double** F1_23; createMatrix(N_y + 2, N_x + 2, F1_23);
	double** F1_24; createMatrix(N_y + 2, N_x + 2, F1_24); double** F1_25; createMatrix(N_y + 2, N_x + 2, F1_25);
	double** F1_31; createMatrix(N_y + 2, N_x + 2, F1_31); double** F1_32; createMatrix(N_y + 2, N_x + 2, F1_32); double** F1_33; createMatrix(N_y + 2, N_x + 2, F1_33);
	double** F1_34; createMatrix(N_y + 2, N_x + 2, F1_34); double** F1_35; createMatrix(N_y + 2, N_x + 2, F1_35);
	double** F1_41; createMatrix(N_y + 2, N_x + 2, F1_41); double** F1_42; createMatrix(N_y + 2, N_x + 2, F1_42); double** F1_43; createMatrix(N_y + 2, N_x + 2, F1_43);
	double** F1_44; createMatrix(N_y + 2, N_x + 2, F1_44); double** F1_45; createMatrix(N_y + 2, N_x + 2, F1_45);
	double** F1_51; createMatrix(N_y + 2, N_x + 2, F1_51); double** F1_52; createMatrix(N_y + 2, N_x + 2, F1_52); double** F1_53; createMatrix(N_y + 2, N_x + 2, F1_53);
	double** F1_54; createMatrix(N_y + 2, N_x + 2, F1_54); double** F1_55; createMatrix(N_y + 2, N_x + 2, F1_55);

	double** F2_11; createMatrix(N_y + 2, N_x + 2, F2_11); double** F2_12; createMatrix(N_y + 2, N_x + 2, F2_12); double** F2_13; createMatrix(N_y + 2, N_x + 2, F2_13);
	double** F2_14; createMatrix(N_y + 2, N_x + 2, F2_14); double** F2_15; createMatrix(N_y + 2, N_x + 2, F2_15);
	double** F2_21; createMatrix(N_y + 2, N_x + 2, F2_21); double** F2_22; createMatrix(N_y + 2, N_x + 2, F2_22); double** F2_23; createMatrix(N_y + 2, N_x + 2, F2_23);
	double** F2_24; createMatrix(N_y + 2, N_x + 2, F2_24); double** F2_25; createMatrix(N_y + 2, N_x + 2, F2_25);
	double** F2_31; createMatrix(N_y + 2, N_x + 2, F2_31); double** F2_32; createMatrix(N_y + 2, N_x + 2, F2_32); double** F2_33; createMatrix(N_y + 2, N_x + 2, F2_33);
	double** F2_34; createMatrix(N_y + 2, N_x + 2, F2_34); double** F2_35; createMatrix(N_y + 2, N_x + 2, F2_35);
	double** F2_41; createMatrix(N_y + 2, N_x + 2, F2_41); double** F2_42; createMatrix(N_y + 2, N_x + 2, F2_42); double** F2_43; createMatrix(N_y + 2, N_x + 2, F2_43);
	double** F2_44; createMatrix(N_y + 2, N_x + 2, F2_44); double** F2_45; createMatrix(N_y + 2, N_x + 2, F2_45);
	double** F2_51; createMatrix(N_y + 2, N_x + 2, F2_51); double** F2_52; createMatrix(N_y + 2, N_x + 2, F2_52); double** F2_53; createMatrix(N_y + 2, N_x + 2, F2_53);
	double** F2_54; createMatrix(N_y + 2, N_x + 2, F2_54); double** F2_55; createMatrix(N_y + 2, N_x + 2, F2_55);

	double** F1_m1; createMatrix(N_y, N_x, F1_m1); double** F1_m2; createMatrix(N_y, N_x, F1_m2); double** F1_m3; createMatrix(N_y, N_x, F1_m3);
	double** F1_m4; createMatrix(N_y, N_x, F1_m4); double** F1_m5; createMatrix(N_y, N_x, F1_m5);
	double** F1_p1; createMatrix(N_y, N_x, F1_p1); double** F1_p2; createMatrix(N_y, N_x, F1_p2); double** F1_p3; createMatrix(N_y, N_x, F1_p3);
	double** F1_p4; createMatrix(N_y, N_x, F1_p4); double** F1_p5; createMatrix(N_y, N_x, F1_p5);

	double** F2_m1; createMatrix(N_y, N_x, F2_m1); double** F2_m2; createMatrix(N_y, N_x, F2_m2); double** F2_m3; createMatrix(N_y, N_x, F2_m3);
	double** F2_m4; createMatrix(N_y, N_x, F2_m4); double** F2_m5; createMatrix(N_y, N_x, F2_m5);
	double** F2_p1; createMatrix(N_y, N_x, F2_p1); double** F2_p2; createMatrix(N_y, N_x, F2_p2); double** F2_p3; createMatrix(N_y, N_x, F2_p3);
	double** F2_p4; createMatrix(N_y, N_x, F2_p4); double** F2_p5; createMatrix(N_y, N_x, F2_p5);

	double** F1m_GQ_0; createMatrix(N_y, N_x, F1m_GQ_0); double** F1m_GQ_1; createMatrix(N_y, N_x, F1m_GQ_1);
	double** F1m_GQ_2; createMatrix(N_y, N_x, F1m_GQ_2);
	double** F1p_GQ_0; createMatrix(N_y, N_x, F1p_GQ_0); double** F1p_GQ_1; createMatrix(N_y, N_x, F1p_GQ_1);

	double** F2m_GQ_0; createMatrix(N_y, N_x, F2m_GQ_0); double** F2m_GQ_1; createMatrix(N_y, N_x, F2m_GQ_1);
	double** F2m_GQ_2; createMatrix(N_y, N_x, F2m_GQ_2);
	double** F2p_GQ_0; createMatrix(N_y, N_x, F2p_GQ_0); double** F2p_GQ_1; createMatrix(N_y, N_x, F2p_GQ_1);

	double** F1_GQ_0; createMatrix(N_y, N_x, F1_GQ_0); double** F1_GQ_1; createMatrix(N_y, N_x, F1_GQ_1);
	double** F1_GQ_2; createMatrix(N_y, N_x, F1_GQ_2);

	double** F2_GQ_0; createMatrix(N_y, N_x, F2_GQ_0); double** F2_GQ_1; createMatrix(N_y, N_x, F2_GQ_1);
	double** F2_GQ_2; createMatrix(N_y, N_x, F2_GQ_2);

	F1_m2_hat(F1_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac);
	F1_m2_hat(F1_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac);
	F1_m2_hat(F1_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac);
	F1_m2_hat(F1_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac);
	F1_m2_hat(F1_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac);

	F2_m2_hat(F2_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F2_alphac);
	F2_m2_hat(F2_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F2_alphac);
	F2_m2_hat(F2_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F2_alphac);
	F2_m2_hat(F2_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F2_alphac);
	F2_m2_hat(F2_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F2_alphac);


	F2_m2_fun(F2_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F2_m2_fun(F2_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F2_m2_fun(F2_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F2_m2_fun(F2_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	F2_m2_fun(F2_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	F2_m2_fun(F2_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);


	F_m(F1_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F1_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F1_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F1_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	F_m(F1_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	F_m(F1_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			F1_m1[i][j] = F1_hat1[i][j + 1] - F1_hat1[i][j]; F1_m2[i][j] = F1_hat2[i][j + 1] - F1_hat2[i][j];
			F1_m3[i][j] = F1_hat3[i][j + 1] - F1_hat3[i][j]; F1_m4[i][j] = F1_hat4[i][j + 1] - F1_hat4[i][j];
			F1_m5[i][j] = F1_hat5[i][j + 1] - F1_hat5[i][j];

			F1_p1[i][j] = F1_hat1[i][j + 1] + F1_hat1[i][j]; F1_p2[i][j] = F1_hat2[i][j + 1] + F1_hat2[i][j];
			F1_p3[i][j] = F1_hat3[i][j + 1] + F1_hat3[i][j]; F1_p4[i][j] = F1_hat4[i][j + 1] + F1_hat4[i][j];
			F1_p5[i][j] = F1_hat5[i][j + 1] + F1_hat5[i][j];
		}
	}

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y; i++)
		{
			F2_m1[i][j] = F2_hat1[i + 1][j] - F2_hat1[i][j]; F2_m2[i][j] = F2_hat2[i + 1][j] - F2_hat2[i][j];
			F2_m3[i][j] = F2_hat3[i + 1][j] - F2_hat3[i][j]; F2_m4[i][j] = F2_hat4[i + 1][j] - F2_hat4[i][j];
			F2_m5[i][j] = F2_hat5[i + 1][j] - F2_hat5[i][j];

			F2_p1[i][j] = F2_hat1[i + 1][j] + F2_hat1[i][j]; F2_p2[i][j] = F2_hat2[i + 1][j] + F2_hat2[i][j];
			F2_p3[i][j] = F2_hat3[i + 1][j] + F2_hat3[i][j]; F2_p4[i][j] = F2_hat4[i + 1][j] + F2_hat4[i][j];
			F2_p5[i][j] = F2_hat5[i + 1][j] + F2_hat5[i][j];
		}
	}

	getGaussquadrature1D_1(F1m_GQ_0, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5); getGaussquadrature1D_p(F1m_GQ_1, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_p2(F1m_GQ_2, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_1(F1p_GQ_0, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5); getGaussquadrature1D_p(F1p_GQ_1, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5);

	getGaussquadrature1D_1(F2m_GQ_0, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5); getGaussquadrature1D_p(F2m_GQ_1, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_p2(F2m_GQ_2, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_1(F2p_GQ_0, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5); getGaussquadrature1D_p(F2p_GQ_1, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5);


	get_Gauss_quadrature_prime_2D_F1(F1_GQ_0, F1_GQ_1, F1_GQ_2, F1_11, F1_21, F1_31, F1_41, F1_51, F1_12, F1_22, F1_32, F1_42, F1_52, F1_13, F1_23, F1_33, F1_43, F1_53, F1_14, F1_24, F1_34, F1_44, F1_54, F1_15, F1_25, F1_35, F1_45, F1_55);
	get_Gauss_quadrature_prime_2D_F2(F2_GQ_0, F2_GQ_1, F2_GQ_2, F2_11, F2_21, F2_31, F2_41, F2_51, F2_12, F2_22, F2_32, F2_42, F2_52, F2_13, F2_23, F2_33, F2_43, F2_53, F2_14, F2_24, F2_34, F2_44, F2_54, F2_15, F2_25, F2_35, F2_45, F2_55);


	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			con_0[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_0[i][j];
			con_1[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_1[i][j] + h_y / 2.0 * F1_GQ_0[i][j];
			con_2[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_0[i][j] + h_x / 2.0 * F2_GQ_0[i][j];
			con_3[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_1[i][j] + h_y / 2.0 * F1_GQ_2[i][j] + h_x / 2.0 * F2_GQ_1[i][j];
			con_4[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_2[i][j] + 3.0 * h_y / 2.0 * F1_GQ_1[i][j];
			con_5[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_2[i][j] - h_x / 2.0 * F2m_GQ_0[i][j] + 3.0 * h_x / 2.0 * F2_GQ_2[i][j];
		}
	}

	deleteMatrix(N_y, F1_hat1); deleteMatrix(N_y, F1_hat2); deleteMatrix(N_y, F1_hat3);
	deleteMatrix(N_y, F1_hat4); deleteMatrix(N_y, F1_hat5);

	deleteMatrix(N_y + 1, F2_hat1); deleteMatrix(N_y + 1, F2_hat2); deleteMatrix(N_y + 1, F2_hat3);
	deleteMatrix(N_y + 1, F2_hat4); deleteMatrix(N_y + 1, F2_hat5);

	deleteMatrix(N_y + 2, F1_11); deleteMatrix(N_y + 2, F1_12); deleteMatrix(N_y + 2, F1_13); deleteMatrix(N_y + 2, F1_14); deleteMatrix(N_y + 2, F1_15);
	deleteMatrix(N_y + 2, F1_21); deleteMatrix(N_y + 2, F1_22); deleteMatrix(N_y + 2, F1_23); deleteMatrix(N_y + 2, F1_24); deleteMatrix(N_y + 2, F1_25);
	deleteMatrix(N_y + 2, F1_31); deleteMatrix(N_y + 2, F1_32); deleteMatrix(N_y + 2, F1_33); deleteMatrix(N_y + 2, F1_34); deleteMatrix(N_y + 2, F1_35);
	deleteMatrix(N_y + 2, F1_41); deleteMatrix(N_y + 2, F1_42); deleteMatrix(N_y + 2, F1_43); deleteMatrix(N_y + 2, F1_44); deleteMatrix(N_y + 2, F1_45);
	deleteMatrix(N_y + 2, F1_51); deleteMatrix(N_y + 2, F1_52); deleteMatrix(N_y + 2, F1_53); deleteMatrix(N_y + 2, F1_54); deleteMatrix(N_y + 2, F1_55);

	deleteMatrix(N_y + 2, F2_11); deleteMatrix(N_y + 2, F2_12); deleteMatrix(N_y + 2, F2_13); deleteMatrix(N_y + 2, F2_14); deleteMatrix(N_y + 2, F2_15);
	deleteMatrix(N_y + 2, F2_21); deleteMatrix(N_y + 2, F2_22); deleteMatrix(N_y + 2, F2_23); deleteMatrix(N_y + 2, F2_24); deleteMatrix(N_y + 2, F2_25);
	deleteMatrix(N_y + 2, F2_31); deleteMatrix(N_y + 2, F2_32); deleteMatrix(N_y + 2, F2_33); deleteMatrix(N_y + 2, F2_34); deleteMatrix(N_y + 2, F2_35);
	deleteMatrix(N_y + 2, F2_41); deleteMatrix(N_y + 2, F2_42); deleteMatrix(N_y + 2, F2_43); deleteMatrix(N_y + 2, F2_44); deleteMatrix(N_y + 2, F2_45);
	deleteMatrix(N_y + 2, F2_51); deleteMatrix(N_y + 2, F2_52); deleteMatrix(N_y + 2, F2_53); deleteMatrix(N_y + 2, F2_54); deleteMatrix(N_y + 2, F2_55);

	deleteMatrix(N_y, F1_m1); deleteMatrix(N_y, F1_m2); deleteMatrix(N_y, F1_m3); deleteMatrix(N_y, F1_m4); deleteMatrix(N_y, F1_m5);
	deleteMatrix(N_y, F1_p1); deleteMatrix(N_y, F1_p2); deleteMatrix(N_y, F1_p3); deleteMatrix(N_y, F1_p4); deleteMatrix(N_y, F1_p5);

	deleteMatrix(N_y, F2_m1); deleteMatrix(N_y, F2_m2); deleteMatrix(N_y, F2_m3); deleteMatrix(N_y, F2_m4); deleteMatrix(N_y, F2_m5);
	deleteMatrix(N_y, F2_p1); deleteMatrix(N_y, F2_p2); deleteMatrix(N_y, F2_p3); deleteMatrix(N_y, F2_p4); deleteMatrix(N_y, F2_p5);

	deleteMatrix(N_y, F1m_GQ_0); deleteMatrix(N_y, F1m_GQ_1); deleteMatrix(N_y, F1m_GQ_2);
	deleteMatrix(N_y, F1p_GQ_0); deleteMatrix(N_y, F1p_GQ_1);

	deleteMatrix(N_y, F2m_GQ_0); deleteMatrix(N_y, F2m_GQ_1); deleteMatrix(N_y, F2m_GQ_2);
	deleteMatrix(N_y, F2p_GQ_0); deleteMatrix(N_y, F2p_GQ_1);

	deleteMatrix(N_y, F1_GQ_0); deleteMatrix(N_y, F1_GQ_1); deleteMatrix(N_y, F1_GQ_2);
	deleteMatrix(N_y, F2_GQ_0); deleteMatrix(N_y, F2_GQ_1); deleteMatrix(N_y, F2_GQ_2);
}

void convection_Ealpha(double** con_0, double** con_1, double** con_2, double** con_3, double** con_4, double** con_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double gamma_alpha, double** F1_alphac, double** F2_alphac, double h_x, double h_y)
{
	double** F1_hat1; createMatrix(N_y, N_x + 1, F1_hat1); double** F1_hat2; createMatrix(N_y, N_x + 1, F1_hat2);
	double** F1_hat3; createMatrix(N_y, N_x + 1, F1_hat3); double** F1_hat4; createMatrix(N_y, N_x + 1, F1_hat4);
	double** F1_hat5; createMatrix(N_y, N_x + 1, F1_hat5);

	double** F2_hat1; createMatrix(N_y + 1, N_x, F2_hat1); double** F2_hat2; createMatrix(N_y + 1, N_x, F2_hat2);
	double** F2_hat3; createMatrix(N_y + 1, N_x, F2_hat3); double** F2_hat4; createMatrix(N_y + 1, N_x, F2_hat4);
	double** F2_hat5; createMatrix(N_y + 1, N_x, F2_hat5);

	double** F1_11; createMatrix(N_y + 2, N_x + 2, F1_11); double** F1_12; createMatrix(N_y + 2, N_x + 2, F1_12); double** F1_13; createMatrix(N_y + 2, N_x + 2, F1_13);
	double** F1_14; createMatrix(N_y + 2, N_x + 2, F1_14); double** F1_15; createMatrix(N_y + 2, N_x + 2, F1_15);
	double** F1_21; createMatrix(N_y + 2, N_x + 2, F1_21); double** F1_22; createMatrix(N_y + 2, N_x + 2, F1_22); double** F1_23; createMatrix(N_y + 2, N_x + 2, F1_23);
	double** F1_24; createMatrix(N_y + 2, N_x + 2, F1_24); double** F1_25; createMatrix(N_y + 2, N_x + 2, F1_25);
	double** F1_31; createMatrix(N_y + 2, N_x + 2, F1_31); double** F1_32; createMatrix(N_y + 2, N_x + 2, F1_32); double** F1_33; createMatrix(N_y + 2, N_x + 2, F1_33);
	double** F1_34; createMatrix(N_y + 2, N_x + 2, F1_34); double** F1_35; createMatrix(N_y + 2, N_x + 2, F1_35);
	double** F1_41; createMatrix(N_y + 2, N_x + 2, F1_41); double** F1_42; createMatrix(N_y + 2, N_x + 2, F1_42); double** F1_43; createMatrix(N_y + 2, N_x + 2, F1_43);
	double** F1_44; createMatrix(N_y + 2, N_x + 2, F1_44); double** F1_45; createMatrix(N_y + 2, N_x + 2, F1_45);
	double** F1_51; createMatrix(N_y + 2, N_x + 2, F1_51); double** F1_52; createMatrix(N_y + 2, N_x + 2, F1_52); double** F1_53; createMatrix(N_y + 2, N_x + 2, F1_53);
	double** F1_54; createMatrix(N_y + 2, N_x + 2, F1_54); double** F1_55; createMatrix(N_y + 2, N_x + 2, F1_55);

	double** F2_11; createMatrix(N_y + 2, N_x + 2, F2_11); double** F2_12; createMatrix(N_y + 2, N_x + 2, F2_12); double** F2_13; createMatrix(N_y + 2, N_x + 2, F2_13);
	double** F2_14; createMatrix(N_y + 2, N_x + 2, F2_14); double** F2_15; createMatrix(N_y + 2, N_x + 2, F2_15);
	double** F2_21; createMatrix(N_y + 2, N_x + 2, F2_21); double** F2_22; createMatrix(N_y + 2, N_x + 2, F2_22); double** F2_23; createMatrix(N_y + 2, N_x + 2, F2_23);
	double** F2_24; createMatrix(N_y + 2, N_x + 2, F2_24); double** F2_25; createMatrix(N_y + 2, N_x + 2, F2_25);
	double** F2_31; createMatrix(N_y + 2, N_x + 2, F2_31); double** F2_32; createMatrix(N_y + 2, N_x + 2, F2_32); double** F2_33; createMatrix(N_y + 2, N_x + 2, F2_33);
	double** F2_34; createMatrix(N_y + 2, N_x + 2, F2_34); double** F2_35; createMatrix(N_y + 2, N_x + 2, F2_35);
	double** F2_41; createMatrix(N_y + 2, N_x + 2, F2_41); double** F2_42; createMatrix(N_y + 2, N_x + 2, F2_42); double** F2_43; createMatrix(N_y + 2, N_x + 2, F2_43);
	double** F2_44; createMatrix(N_y + 2, N_x + 2, F2_44); double** F2_45; createMatrix(N_y + 2, N_x + 2, F2_45);
	double** F2_51; createMatrix(N_y + 2, N_x + 2, F2_51); double** F2_52; createMatrix(N_y + 2, N_x + 2, F2_52); double** F2_53; createMatrix(N_y + 2, N_x + 2, F2_53);
	double** F2_54; createMatrix(N_y + 2, N_x + 2, F2_54); double** F2_55; createMatrix(N_y + 2, N_x + 2, F2_55);

	double** F1_m1; createMatrix(N_y, N_x, F1_m1); double** F1_m2; createMatrix(N_y, N_x, F1_m2); double** F1_m3; createMatrix(N_y, N_x, F1_m3);
	double** F1_m4; createMatrix(N_y, N_x, F1_m4); double** F1_m5; createMatrix(N_y, N_x, F1_m5);
	double** F1_p1; createMatrix(N_y, N_x, F1_p1); double** F1_p2; createMatrix(N_y, N_x, F1_p2); double** F1_p3; createMatrix(N_y, N_x, F1_p3);
	double** F1_p4; createMatrix(N_y, N_x, F1_p4); double** F1_p5; createMatrix(N_y, N_x, F1_p5);

	double** F2_m1; createMatrix(N_y, N_x, F2_m1); double** F2_m2; createMatrix(N_y, N_x, F2_m2); double** F2_m3; createMatrix(N_y, N_x, F2_m3);
	double** F2_m4; createMatrix(N_y, N_x, F2_m4); double** F2_m5; createMatrix(N_y, N_x, F2_m5);
	double** F2_p1; createMatrix(N_y, N_x, F2_p1); double** F2_p2; createMatrix(N_y, N_x, F2_p2); double** F2_p3; createMatrix(N_y, N_x, F2_p3);
	double** F2_p4; createMatrix(N_y, N_x, F2_p4); double** F2_p5; createMatrix(N_y, N_x, F2_p5);

	double** F1m_GQ_0; createMatrix(N_y, N_x, F1m_GQ_0); double** F1m_GQ_1; createMatrix(N_y, N_x, F1m_GQ_1);
	double** F1m_GQ_2; createMatrix(N_y, N_x, F1m_GQ_2);
	double** F1p_GQ_0; createMatrix(N_y, N_x, F1p_GQ_0); double** F1p_GQ_1; createMatrix(N_y, N_x, F1p_GQ_1);

	double** F2m_GQ_0; createMatrix(N_y, N_x, F2m_GQ_0); double** F2m_GQ_1; createMatrix(N_y, N_x, F2m_GQ_1);
	double** F2m_GQ_2; createMatrix(N_y, N_x, F2m_GQ_2);
	double** F2p_GQ_0; createMatrix(N_y, N_x, F2p_GQ_0); double** F2p_GQ_1; createMatrix(N_y, N_x, F2p_GQ_1);

	double** F1_GQ_0; createMatrix(N_y, N_x, F1_GQ_0); double** F1_GQ_1; createMatrix(N_y, N_x, F1_GQ_1);
	double** F1_GQ_2; createMatrix(N_y, N_x, F1_GQ_2);

	double** F2_GQ_0; createMatrix(N_y, N_x, F2_GQ_0); double** F2_GQ_1; createMatrix(N_y, N_x, F2_GQ_1);
	double** F2_GQ_2; createMatrix(N_y, N_x, F2_GQ_2);

	F1_alpha_hat(F1_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F1_alphac, gamma_alpha);
	F1_alpha_hat(F1_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F1_alphac, gamma_alpha);
	F1_alpha_hat(F1_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F1_alphac, gamma_alpha);
	F1_alpha_hat(F1_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F1_alphac, gamma_alpha);
	F1_alpha_hat(F1_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F1_alphac, gamma_alpha);

	F2_alpha_hat(F2_hat1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F2_alphac, gamma_alpha);
	F2_alpha_hat(F2_hat2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F2_alphac, gamma_alpha);
	F2_alpha_hat(F2_hat3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F2_alphac, gamma_alpha);
	F2_alpha_hat(F2_hat4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F2_alphac, gamma_alpha);
	F2_alpha_hat(F2_hat5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, F2_alphac, gamma_alpha);

	F1_alpha(F1_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F1_alpha(F1_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F1_alpha(F1_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F1_alpha(F1_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F1_alpha(F1_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F1_alpha(F1_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);


	F2_alpha(F2_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F2_alpha(F2_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F2_alpha(F2_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F2_alpha(F2_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	F2_alpha(F2_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);
	F2_alpha(F2_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5, gamma_alpha);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			F1_m1[i][j] = F1_hat1[i][j + 1] - F1_hat1[i][j]; F1_m2[i][j] = F1_hat2[i][j + 1] - F1_hat2[i][j];
			F1_m3[i][j] = F1_hat3[i][j + 1] - F1_hat3[i][j]; F1_m4[i][j] = F1_hat4[i][j + 1] - F1_hat4[i][j];
			F1_m5[i][j] = F1_hat5[i][j + 1] - F1_hat5[i][j];

			F1_p1[i][j] = F1_hat1[i][j + 1] + F1_hat1[i][j]; F1_p2[i][j] = F1_hat2[i][j + 1] + F1_hat2[i][j];
			F1_p3[i][j] = F1_hat3[i][j + 1] + F1_hat3[i][j]; F1_p4[i][j] = F1_hat4[i][j + 1] + F1_hat4[i][j];
			F1_p5[i][j] = F1_hat5[i][j + 1] + F1_hat5[i][j];
		}
	}

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y; i++)
		{
			F2_m1[i][j] = F2_hat1[i + 1][j] - F2_hat1[i][j]; F2_m2[i][j] = F2_hat2[i + 1][j] - F2_hat2[i][j];
			F2_m3[i][j] = F2_hat3[i + 1][j] - F2_hat3[i][j]; F2_m4[i][j] = F2_hat4[i + 1][j] - F2_hat4[i][j];
			F2_m5[i][j] = F2_hat5[i + 1][j] - F2_hat5[i][j];

			F2_p1[i][j] = F2_hat1[i + 1][j] + F2_hat1[i][j]; F2_p2[i][j] = F2_hat2[i + 1][j] + F2_hat2[i][j];
			F2_p3[i][j] = F2_hat3[i + 1][j] + F2_hat3[i][j]; F2_p4[i][j] = F2_hat4[i + 1][j] + F2_hat4[i][j];
			F2_p5[i][j] = F2_hat5[i + 1][j] + F2_hat5[i][j];
		}
	}

	getGaussquadrature1D_1(F1m_GQ_0, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5); getGaussquadrature1D_p(F1m_GQ_1, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_p2(F1m_GQ_2, F1_m1, F1_m2, F1_m3, F1_m4, F1_m5);
	getGaussquadrature1D_1(F1p_GQ_0, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5); getGaussquadrature1D_p(F1p_GQ_1, F1_p1, F1_p2, F1_p3, F1_p4, F1_p5);

	getGaussquadrature1D_1(F2m_GQ_0, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5); getGaussquadrature1D_p(F2m_GQ_1, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_p2(F2m_GQ_2, F2_m1, F2_m2, F2_m3, F2_m4, F2_m5);
	getGaussquadrature1D_1(F2p_GQ_0, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5); getGaussquadrature1D_p(F2p_GQ_1, F2_p1, F2_p2, F2_p3, F2_p4, F2_p5);


	get_Gauss_quadrature_prime_2D_F1(F1_GQ_0, F1_GQ_1, F1_GQ_2, F1_11, F1_21, F1_31, F1_41, F1_51, F1_12, F1_22, F1_32, F1_42, F1_52, F1_13, F1_23, F1_33, F1_43, F1_53, F1_14, F1_24, F1_34, F1_44, F1_54, F1_15, F1_25, F1_35, F1_45, F1_55);
	get_Gauss_quadrature_prime_2D_F2(F2_GQ_0, F2_GQ_1, F2_GQ_2, F2_11, F2_21, F2_31, F2_41, F2_51, F2_12, F2_22, F2_32, F2_42, F2_52, F2_13, F2_23, F2_33, F2_43, F2_53, F2_14, F2_24, F2_34, F2_44, F2_54, F2_15, F2_25, F2_35, F2_45, F2_55);


	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			con_0[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_0[i][j];
			con_1[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_1[i][j] + h_y / 2.0 * F1_GQ_0[i][j];
			con_2[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_0[i][j] + h_x / 2.0 * F2_GQ_0[i][j];
			con_3[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1p_GQ_1[i][j] - h_x / 2.0 * F2p_GQ_1[i][j] + h_y / 2.0 * F1_GQ_2[i][j] + h_x / 2.0 * F2_GQ_1[i][j];
			con_4[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_0[i][j] - h_x / 2.0 * F2m_GQ_2[i][j] + 3.0 * h_y / 2.0 * F1_GQ_1[i][j];
			con_5[i + 1][j + 1] = -1.0 * h_y / 2.0 * F1m_GQ_2[i][j] - h_x / 2.0 * F2m_GQ_0[i][j] + 3.0 * h_x / 2.0 * F2_GQ_2[i][j];
		}
	}

	deleteMatrix(N_y, F1_hat1); deleteMatrix(N_y, F1_hat2); deleteMatrix(N_y, F1_hat3);
	deleteMatrix(N_y, F1_hat4); deleteMatrix(N_y, F1_hat5);

	deleteMatrix(N_y + 1, F2_hat1); deleteMatrix(N_y + 1, F2_hat2); deleteMatrix(N_y + 1, F2_hat3);
	deleteMatrix(N_y + 1, F2_hat4); deleteMatrix(N_y + 1, F2_hat5);

	deleteMatrix(N_y + 2, F1_11); deleteMatrix(N_y + 2, F1_12); deleteMatrix(N_y + 2, F1_13); deleteMatrix(N_y + 2, F1_14); deleteMatrix(N_y + 2, F1_15);
	deleteMatrix(N_y + 2, F1_21); deleteMatrix(N_y + 2, F1_22); deleteMatrix(N_y + 2, F1_23); deleteMatrix(N_y + 2, F1_24); deleteMatrix(N_y + 2, F1_25);
	deleteMatrix(N_y + 2, F1_31); deleteMatrix(N_y + 2, F1_32); deleteMatrix(N_y + 2, F1_33); deleteMatrix(N_y + 2, F1_34); deleteMatrix(N_y + 2, F1_35);
	deleteMatrix(N_y + 2, F1_41); deleteMatrix(N_y + 2, F1_42); deleteMatrix(N_y + 2, F1_43); deleteMatrix(N_y + 2, F1_44); deleteMatrix(N_y + 2, F1_45);
	deleteMatrix(N_y + 2, F1_51); deleteMatrix(N_y + 2, F1_52); deleteMatrix(N_y + 2, F1_53); deleteMatrix(N_y + 2, F1_54); deleteMatrix(N_y + 2, F1_55);

	deleteMatrix(N_y + 2, F2_11); deleteMatrix(N_y + 2, F2_12); deleteMatrix(N_y + 2, F2_13); deleteMatrix(N_y + 2, F2_14); deleteMatrix(N_y + 2, F2_15);
	deleteMatrix(N_y + 2, F2_21); deleteMatrix(N_y + 2, F2_22); deleteMatrix(N_y + 2, F2_23); deleteMatrix(N_y + 2, F2_24); deleteMatrix(N_y + 2, F2_25);
	deleteMatrix(N_y + 2, F2_31); deleteMatrix(N_y + 2, F2_32); deleteMatrix(N_y + 2, F2_33); deleteMatrix(N_y + 2, F2_34); deleteMatrix(N_y + 2, F2_35);
	deleteMatrix(N_y + 2, F2_41); deleteMatrix(N_y + 2, F2_42); deleteMatrix(N_y + 2, F2_43); deleteMatrix(N_y + 2, F2_44); deleteMatrix(N_y + 2, F2_45);
	deleteMatrix(N_y + 2, F2_51); deleteMatrix(N_y + 2, F2_52); deleteMatrix(N_y + 2, F2_53); deleteMatrix(N_y + 2, F2_54); deleteMatrix(N_y + 2, F2_55);

	deleteMatrix(N_y, F1_m1); deleteMatrix(N_y, F1_m2); deleteMatrix(N_y, F1_m3); deleteMatrix(N_y, F1_m4); deleteMatrix(N_y, F1_m5);
	deleteMatrix(N_y, F1_p1); deleteMatrix(N_y, F1_p2); deleteMatrix(N_y, F1_p3); deleteMatrix(N_y, F1_p4); deleteMatrix(N_y, F1_p5);

	deleteMatrix(N_y, F2_m1); deleteMatrix(N_y, F2_m2); deleteMatrix(N_y, F2_m3); deleteMatrix(N_y, F2_m4); deleteMatrix(N_y, F2_m5);
	deleteMatrix(N_y, F2_p1); deleteMatrix(N_y, F2_p2); deleteMatrix(N_y, F2_p3); deleteMatrix(N_y, F2_p4); deleteMatrix(N_y, F2_p5);

	deleteMatrix(N_y, F1m_GQ_0); deleteMatrix(N_y, F1m_GQ_1); deleteMatrix(N_y, F1m_GQ_2);
	deleteMatrix(N_y, F1p_GQ_0); deleteMatrix(N_y, F1p_GQ_1);

	deleteMatrix(N_y, F2m_GQ_0); deleteMatrix(N_y, F2m_GQ_1); deleteMatrix(N_y, F2m_GQ_2);
	deleteMatrix(N_y, F2p_GQ_0); deleteMatrix(N_y, F2p_GQ_1);

	deleteMatrix(N_y, F1_GQ_0); deleteMatrix(N_y, F1_GQ_1); deleteMatrix(N_y, F1_GQ_2);
	deleteMatrix(N_y, F2_GQ_0); deleteMatrix(N_y, F2_GQ_1); deleteMatrix(N_y, F2_GQ_2);
}

void compute_F1_alphaC_max(double** alpha_c, double** alphac1, double** alphac2, double** alphac3, double** alphac4, double** alphac5)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			double temp_max = alphac1[i][j];
			if (alphac2[i][j] > temp_max)
			{
				temp_max = alphac2[i][j];
			}
			if (alphac3[i][j] > temp_max)
			{
				temp_max = alphac3[i][j];
			}
			if (alphac4[i][j] > temp_max)
			{
				temp_max = alphac4[i][j];
			}
			if (alphac5[i][j] > temp_max)
			{
				temp_max = alphac5[i][j];
			}
			alpha_c[i][j] = temp_max;
		}
	}
}

void compute_F2_alphaC_max(double** alpha_c, double** alphac1, double** alphac2, double** alphac3, double** alphac4, double** alphac5)
{
	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 1; i++)
		{
			double temp_max = alphac1[i][j];
			if (alphac2[i][j] > temp_max)
			{
				temp_max = alphac2[i][j];
			}
			if (alphac3[i][j] > temp_max)
			{
				temp_max = alphac3[i][j];
			}
			if (alphac4[i][j] > temp_max)
			{
				temp_max = alphac4[i][j];
			}
			if (alphac5[i][j] > temp_max)
			{
				temp_max = alphac5[i][j];
			}
			alpha_c[i][j] = temp_max;
		}
	}
}

double compute_alphaCx_max(double** alphaC_x)
{
	double temp_max = 0.0;
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			if (alphaC_x[i][j] > temp_max)
			{
				temp_max = alphaC_x[i][j];
			}
		}
	}
	return temp_max;
}

double compute_alphaCx_min(double** alphaC_x)
{
	double temp_min = 100.0;
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 1; j++)
		{
			if (alphaC_x[i][j] < temp_min)
			{
				temp_min = alphaC_x[i][j];
			}
		}
	}
	return temp_min;
}

double compute_alphaCy_max(double** alphaC_y)
{
	double temp_max = 0.0;
	for (int i = 0; i < N_y + 1; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			if (alphaC_y[i][j] > temp_max)
			{
				temp_max = alphaC_y[i][j];
			}
		}
	}
	return temp_max;
}

double compute_alphaCy_min(double** alphaC_y)
{
	double temp_min = 100.0;
	for (int i = 0; i < N_y + 1; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			if (alphaC_y[i][j] < temp_min)
			{
				temp_min = alphaC_y[i][j];
			}
		}
	}
	return temp_min;
}


// Nonconservetive
void p_alpha(double** u, double xi, double eta, double gamma_alpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			u[i][j] = (gamma_alpha - 1.0) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0);
		}
	}
}

void p_alphaxi(double** u, double xi, double eta, double gamma_alpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ealphaxi = Ealpha_1[i][j] + Ealpha_3[i][j] * eta + 3.0 * Ealpha_4[i][j] * xi;
			double temp_rhoxi = rho_1[i][j] + rho_3[i][j] * eta + 3.0 * rho_4[i][j] * xi;
			double temp_m1xi = m1_1[i][j] + m1_3[i][j] * eta + 3.0 * m1_4[i][j] * xi;
			double temp_m2xi = m2_1[i][j] + m2_3[i][j] * eta + 3.0 * m2_4[i][j] * xi;
			u[i][j] = (gamma_alpha - 1.0) * (temp_Ealphaxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
		}
	}
}

void p_alphaeta(double** u, double xi, double eta, double gamma_alpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ealphaxi = Ealpha_2[i][j] + Ealpha_3[i][j] * xi + 3.0 * Ealpha_5[i][j] * eta;
			double temp_rhoxi = rho_2[i][j] + rho_3[i][j] * xi + 3.0 * rho_5[i][j] * eta;
			double temp_m1xi = m1_2[i][j] + m1_3[i][j] * xi + 3.0 * m1_5[i][j] * eta;
			double temp_m2xi = m2_2[i][j] + m2_3[i][j] * xi + 3.0 * m2_5[i][j] * eta;
			u[i][j] = (gamma_alpha - 1.0) * (temp_Ealphaxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
		}
	}
}

void F_Ne(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** pe_xi; createMatrix(N_y + 2, N_x + 2, pe_xi); double** pi_xi; createMatrix(N_y + 2, N_x + 2, pi_xi);
	double** pr_xi; createMatrix(N_y + 2, N_x + 2, pr_xi); double** pe_eta; createMatrix(N_y + 2, N_x + 2, pe_eta);
	double** pi_eta; createMatrix(N_y + 2, N_x + 2, pi_eta); double** pr_eta; createMatrix(N_y + 2, N_x + 2, pr_eta);

	p_alphaxi(pe_xi, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaxi(pi_xi, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaxi(pr_xi, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	p_alphaeta(pe_eta, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaeta(pi_eta, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaeta(pr_eta, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_u = temp_m1 / temp_rho;
			double temp_v = temp_m2 / temp_rho;

			u[i][j] = 2.0 / (3.0 * h_x) * temp_u * (2 * pe_xi[i][j] - pi_xi[i][j] - pr_xi[i][j]) + 2.0 / (3.0 * h_y) * temp_v * (2 * pe_eta[i][j] - pi_eta[i][j] - pr_eta[i][j]);
		}
	}

	deleteMatrix(N_y + 2, pe_xi); deleteMatrix(N_y + 2, pi_xi); deleteMatrix(N_y + 2, pr_xi);
	deleteMatrix(N_y + 2, pe_eta); deleteMatrix(N_y + 2, pi_eta); deleteMatrix(N_y + 2, pr_eta);
}

void F_Ni(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** pe_xi; createMatrix(N_y + 2, N_x + 2, pe_xi); double** pi_xi; createMatrix(N_y + 2, N_x + 2, pi_xi);
	double** pr_xi; createMatrix(N_y + 2, N_x + 2, pr_xi); double** pe_eta; createMatrix(N_y + 2, N_x + 2, pe_eta);
	double** pi_eta; createMatrix(N_y + 2, N_x + 2, pi_eta); double** pr_eta; createMatrix(N_y + 2, N_x + 2, pr_eta);

	p_alphaxi(pe_xi, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaxi(pi_xi, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaxi(pr_xi, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	p_alphaeta(pe_eta, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaeta(pi_eta, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaeta(pr_eta, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_u = temp_m1 / temp_rho;
			double temp_v = temp_m2 / temp_rho;

			u[i][j] = 2.0 / (3.0 * h_x) * temp_u * (2 * pi_xi[i][j] - pe_xi[i][j] - pr_xi[i][j]) + 2.0 / (3.0 * h_y) * temp_v * (2 * pi_eta[i][j] - pe_eta[i][j] - pr_eta[i][j]);
		}
	}

	deleteMatrix(N_y + 2, pe_xi); deleteMatrix(N_y + 2, pi_xi); deleteMatrix(N_y + 2, pr_xi);
	deleteMatrix(N_y + 2, pe_eta); deleteMatrix(N_y + 2, pi_eta); deleteMatrix(N_y + 2, pr_eta);
}

void F_Nr(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** pe_xi; createMatrix(N_y + 2, N_x + 2, pe_xi); double** pi_xi; createMatrix(N_y + 2, N_x + 2, pi_xi);
	double** pr_xi; createMatrix(N_y + 2, N_x + 2, pr_xi); double** pe_eta; createMatrix(N_y + 2, N_x + 2, pe_eta);
	double** pi_eta; createMatrix(N_y + 2, N_x + 2, pi_eta); double** pr_eta; createMatrix(N_y + 2, N_x + 2, pr_eta);

	p_alphaxi(pe_xi, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaxi(pi_xi, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaxi(pr_xi, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	p_alphaeta(pe_eta, xi, eta, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	p_alphaeta(pi_eta, xi, eta, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	p_alphaeta(pr_eta, xi, eta, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_u = temp_m1 / temp_rho;
			double temp_v = temp_m2 / temp_rho;

			u[i][j] = 2.0 / (3.0 * h_x) * temp_u * (2 * pr_xi[i][j] - pe_xi[i][j] - pi_xi[i][j]) + 2.0 / (3.0 * h_y) * temp_v * (2 * pr_eta[i][j] - pe_eta[i][j] - pi_eta[i][j]);
		}
	}

	deleteMatrix(N_y + 2, pe_xi); deleteMatrix(N_y + 2, pi_xi); deleteMatrix(N_y + 2, pr_xi);
	deleteMatrix(N_y + 2, pe_eta); deleteMatrix(N_y + 2, pi_eta); deleteMatrix(N_y + 2, pr_eta);
}

void getGaussquadrature(double* U_0, double* U_1, double* U_2, double* u_x1, double* u_x2, double* u_x3, double* u_x4, double* u_x5, int N_index, double h)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = h / 2.0 * (weight_1 * u_x1[i] + weight_2 * u_x2[i] + weight_3 * u_x3[i] + weight_4 * u_x4[i] + weight_5 * u_x5[i]);
		U_1[i] = h / 2.0 * (weight_1 * u_x1[i] * GQxi_1 + weight_2 * u_x2[i] * GQxi_2 + weight_3 * u_x3[i] * GQxi_3 + weight_4 * u_x4[i] * GQxi_4 + weight_5 * u_x5[i] * GQxi_5);
		U_2[i] = h / 2.0 * (weight_1 * u_x1[i] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + weight_2 * u_x2[i] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + weight_3 * u_x3[i] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + weight_4 * u_x4[i] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + weight_5 * u_x5[i] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0);
	}
}

void getL2projection(double* U_0, double* U_1, double* U_2, double* u_0, double* u_1, double* u_2, int N_index, double h)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = u_0[i] / h;
		U_1[i] = 3.0 * u_1[i] / h;
		U_2[i] = 5.0 * u_2[i] / h;
	}
}

void compute_uhalf(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m_0, double** m_1, double** m_2, double** m_3, double** m_4, double** m_5, double h_x, int index_mode)
{
	double** m1; createMatrix(N_y + 2, N_x + 2, m1);
	double** m2; createMatrix(N_y + 2, N_x + 2, m2);
	double** m3; createMatrix(N_y + 2, N_x + 2, m3);
	double** m4; createMatrix(N_y + 2, N_x + 2, m4);
	double** m5; createMatrix(N_y + 2, N_x + 2, m5);

	double** rho1; createMatrix(N_y + 2, N_x + 2, rho1);
	double** rho2; createMatrix(N_y + 2, N_x + 2, rho2);
	double** rho3; createMatrix(N_y + 2, N_x + 2, rho3);
	double** rho4; createMatrix(N_y + 2, N_x + 2, rho4);
	double** rho5; createMatrix(N_y + 2, N_x + 2, rho5);

	F_d(m1, GQxi_1, eta, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m2, GQxi_2, eta, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m3, GQxi_3, eta, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m4, GQxi_4, eta, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m5, GQxi_5, eta, m_0, m_1, m_2, m_3, m_4, m_5);

	F_d(rho1, GQxi_1, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho2, GQxi_2, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho3, GQxi_3, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho4, GQxi_4, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho5, GQxi_5, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);

	// u(N_y, N_x + 1)
	for (int i = 0; i < N_y; i++)
	{
		double* u_1 = new double[N_x + 2]; double* u_2 = new double[N_x + 2]; double* u_3 = new double[N_x + 2]; double* u_4 = new double[N_x + 2]; double* u_5 = new double[N_x + 2];
		double* u_GQ_0 = new double[N_x + 2]; double* u_GQ_1 = new double[N_x + 2]; double* u_GQ_2 = new double[N_x + 2];
		double* omega_0 = new double[N_x + 2]; double* omega_1 = new double[N_x + 2]; double* omega_2 = new double[N_x + 2];

		for (int j = 0; j < N_x + 2; j++)
		{
			u_1[j] = m1[i + 1][j] / rho1[i + 1][j];
			u_2[j] = m2[i + 1][j] / rho2[i + 1][j];
			u_3[j] = m3[i + 1][j] / rho3[i + 1][j];
			u_4[j] = m4[i + 1][j] / rho4[i + 1][j];
			u_5[j] = m5[i + 1][j] / rho5[i + 1][j];
		}

		getGaussquadrature(u_GQ_0, u_GQ_1, u_GQ_2, u_1, u_2, u_3, u_4, u_5, N_x + 2, h_x);
		getL2projection(omega_0, omega_1, omega_2, u_GQ_0, u_GQ_1, u_GQ_2, N_x + 2, h_x);
		if (index_mode == 0)
		{
			for (int j = 0; j < N_x + 1; j++)
			{
				u[i][j] = (omega_0[j] + omega_0[j + 1]) / 2.0;
			}
		}

		if (index_mode == 1)
		{
			for (int j = 0; j < N_x + 1; j++)
			{
				u[i][j] = omega_0[j] + omega_1[j] + (omega_1[j + 1] - omega_1[j]) / 6.0 - (omega_0[j] + omega_1[j] - omega_0[j + 1] + omega_1[j + 1]) / 2.0;
			}
		}

		if (index_mode == 2)
		{
			double* w_3 = new double[N_x + 1]; double* w_4 = new double[N_x + 1]; double* w_5 = new double[N_x + 1];
			for (int j = 0; j < N_x + 1; j++)
			{
				w_3[j] = 7.0 * (omega_1[j + 1] - 3 * omega_2[j + 1] - omega_1[j] - 3 * omega_2[j]) / 64.0 + 5.0 * (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j]) / 6.0 - (omega_2[j + 1] - omega_2[j]) / 15.0;
				w_4[j] = -3.0 * (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j]) / 8.0 + 3.0 * (omega_2[j + 1] - omega_2[j]) / 40.0 - (omega_1[j + 1] - 3 * omega_2[j + 1] - omega_1[j] - 3 * omega_2[j]) / 64.0;
				w_5[j] = (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j] - omega_2[j + 1] / 5.0 + omega_2[j] / 5.0) / 24.0;
			}

			for (int j = 0; j < N_x + 1; j++)
			{
				u[i][j] = omega_0[j] + omega_1[j] + omega_2[j] + w_3[j] + w_4[j] + w_5[j];
			}
			delete[] w_3; delete[] w_4; delete[] w_5;
		}

		delete[] u_1; delete[] u_2; delete[] u_3; delete[] u_4; delete[] u_5;
		delete[] u_GQ_0; delete[] u_GQ_1; delete[] u_GQ_2;
		delete[] omega_0; delete[] omega_1; delete[] omega_2;
	}
	deleteMatrix(N_y + 2, m1); deleteMatrix(N_y + 2, m2); deleteMatrix(N_y + 2, m3); deleteMatrix(N_y + 2, m4); deleteMatrix(N_y + 2, m5);
	deleteMatrix(N_y + 2, rho1); deleteMatrix(N_y + 2, rho2); deleteMatrix(N_y + 2, rho3); deleteMatrix(N_y + 2, rho4); deleteMatrix(N_y + 2, rho5);
}

void compute_vhalf(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m_0, double** m_1, double** m_2, double** m_3, double** m_4, double** m_5, double h_y, int index_mode)
{
	double** m1; createMatrix(N_y + 2, N_x + 2, m1);
	double** m2; createMatrix(N_y + 2, N_x + 2, m2);
	double** m3; createMatrix(N_y + 2, N_x + 2, m3);
	double** m4; createMatrix(N_y + 2, N_x + 2, m4);
	double** m5; createMatrix(N_y + 2, N_x + 2, m5);

	double** rho1; createMatrix(N_y + 2, N_x + 2, rho1);
	double** rho2; createMatrix(N_y + 2, N_x + 2, rho2);
	double** rho3; createMatrix(N_y + 2, N_x + 2, rho3);
	double** rho4; createMatrix(N_y + 2, N_x + 2, rho4);
	double** rho5; createMatrix(N_y + 2, N_x + 2, rho5);

	F_d(m1, xi, GQxi_1, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m2, xi, GQxi_2, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m3, xi, GQxi_3, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m4, xi, GQxi_4, m_0, m_1, m_2, m_3, m_4, m_5);
	F_d(m5, xi, GQxi_5, m_0, m_1, m_2, m_3, m_4, m_5);

	F_d(rho1, xi, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho2, xi, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho3, xi, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho4, xi, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);
	F_d(rho5, xi, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);

	// u(N_y + 1, N_x)
	for (int i = 0; i < N_x; i++)
	{
		double* u_1 = new double[N_y + 2]; double* u_2 = new double[N_y + 2]; double* u_3 = new double[N_y + 2]; double* u_4 = new double[N_y + 2]; double* u_5 = new double[N_y + 2];
		double* u_GQ_0 = new double[N_y + 2]; double* u_GQ_1 = new double[N_y + 2]; double* u_GQ_2 = new double[N_y + 2];
		double* omega_0 = new double[N_y + 2]; double* omega_1 = new double[N_y + 2]; double* omega_2 = new double[N_y + 2];

		for (int j = 0; j < N_y + 2; j++)
		{
			u_1[j] = m1[j][i + 1] / rho1[j][i + 1];
			u_2[j] = m2[j][i + 1] / rho2[j][i + 1];
			u_3[j] = m3[j][i + 1] / rho3[j][i + 1];
			u_4[j] = m4[j][i + 1] / rho4[j][i + 1];
			u_5[j] = m5[j][i + 1] / rho5[j][i + 1];
		}

		getGaussquadrature(u_GQ_0, u_GQ_1, u_GQ_2, u_1, u_2, u_3, u_4, u_5, N_y + 2, h_y);
		getL2projection(omega_0, omega_1, omega_2, u_GQ_0, u_GQ_1, u_GQ_2, N_y + 2, h_y);
		if (index_mode == 0)
		{
			for (int j = 0; j < N_y + 1; j++)
			{
				u[j][i] = (omega_0[j] + omega_0[j + 1]) / 2.0;
			}
		}

		if (index_mode == 1)
		{
			for (int j = 0; j < N_y + 1; j++)
			{
				u[j][i] = omega_0[j] + omega_1[j] + (omega_1[j + 1] - omega_1[j]) / 6.0 - (omega_0[j] + omega_1[j] - omega_0[j + 1] + omega_1[j + 1]) / 2.0;
			}
		}

		if (index_mode == 2)
		{
			double* w_3 = new double[N_y + 1]; double* w_4 = new double[N_y + 1]; double* w_5 = new double[N_y + 1];
			for (int j = 0; j < N_y + 1; j++)
			{
				w_3[j] = 7.0 * (omega_1[j + 1] - 3 * omega_2[j + 1] - omega_1[j] - 3 * omega_2[j]) / 64.0 + 5.0 * (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j]) / 6.0 - (omega_2[j + 1] - omega_2[j]) / 15.0;
				w_4[j] = -3.0 * (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j]) / 8.0 + 3.0 * (omega_2[j + 1] - omega_2[j]) / 40.0 - (omega_1[j + 1] - 3 * omega_2[j + 1] - omega_1[j] - 3 * omega_2[j]) / 64.0;
				w_5[j] = (omega_0[j + 1] - omega_1[j + 1] + omega_2[j + 1] - omega_0[j] - omega_1[j] - omega_2[j] - omega_2[j + 1] / 5.0 + omega_2[j] / 5.0) / 24.0;
			}

			for (int j = 0; j < N_y + 1; j++)
			{
				u[j][i] = omega_0[j] + omega_1[j] + omega_2[j] + w_3[j] + w_4[j] + w_5[j];
			}
			delete[] w_3; delete[] w_4; delete[] w_5;
		}

		delete[] u_1; delete[] u_2; delete[] u_3; delete[] u_4; delete[] u_5;
		delete[] u_GQ_0; delete[] u_GQ_1; delete[] u_GQ_2;
		delete[] omega_0; delete[] omega_1; delete[] omega_2;
	}
	deleteMatrix(N_y + 2, m1); deleteMatrix(N_y + 2, m2); deleteMatrix(N_y + 2, m3); deleteMatrix(N_y + 2, m4); deleteMatrix(N_y + 2, m5);
	deleteMatrix(N_y + 2, rho1); deleteMatrix(N_y + 2, rho2); deleteMatrix(N_y + 2, rho3); deleteMatrix(N_y + 2, rho4); deleteMatrix(N_y + 2, rho5);
}

void compute_coefx(double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double h_x)
{
	double** u_half1; createMatrix(N_y, N_x + 1, u_half1);
	double** u_half2; createMatrix(N_y, N_x + 1, u_half2);
	double** u_half3; createMatrix(N_y, N_x + 1, u_half3);
	double** u_half4; createMatrix(N_y, N_x + 1, u_half4);
	double** u_half5; createMatrix(N_y, N_x + 1, u_half5);

	compute_uhalf(u_half1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, Mode);
	compute_uhalf(u_half2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, Mode);
	compute_uhalf(u_half3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, Mode);
	compute_uhalf(u_half4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, Mode);
	compute_uhalf(u_half5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, Mode);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_L1 = -1.0 * u_half1[i][j] / 3.0;
			if (temp_L1 < 0)
			{
				temp_L1 = 0.0;

			}
			coef_L1[i][j] = temp_L1;

			double temp_L2 = -1.0 * u_half2[i][j] / 3.0;
			if (temp_L2 < 0)
			{
				temp_L2 = 0.0;

			}
			coef_L2[i][j] = temp_L2;

			double temp_L3 = -1.0 * u_half3[i][j] / 3.0;
			if (temp_L3 < 0)
			{
				temp_L3 = 0.0;

			}
			coef_L3[i][j] = temp_L3;

			double temp_L4 = -1.0 * u_half4[i][j] / 3.0;
			if (temp_L4 < 0)
			{
				temp_L4 = 0.0;

			}
			coef_L4[i][j] = temp_L4;

			double temp_L5 = -1.0 * u_half5[i][j] / 3.0;
			if (temp_L5 < 0)
			{
				temp_L5 = 0.0;

			}
			coef_L5[i][j] = temp_L5;

			double temp_R1 = -1.0 * u_half1[i][j + 1] / 3.0;
			if (temp_R1 > 0)
			{
				temp_R1 = 0.0;
			}
			coef_R1[i][j] = temp_R1;

			double temp_R2 = -1.0 * u_half2[i][j + 1] / 3.0;
			if (temp_R2 > 0)
			{
				temp_R2 = 0.0;
			}
			coef_R2[i][j] = temp_R2;

			double temp_R3 = -1.0 * u_half3[i][j + 1] / 3.0;
			if (temp_R3 > 0)
			{
				temp_R3 = 0.0;
			}
			coef_R3[i][j] = temp_R3;

			double temp_R4 = -1.0 * u_half4[i][j + 1] / 3.0;
			if (temp_R4 > 0)
			{
				temp_R4 = 0.0;
			}
			coef_R4[i][j] = temp_R4;

			double temp_R5 = -1.0 * u_half5[i][j + 1] / 3.0;
			if (temp_R5 > 0)
			{
				temp_R5 = 0.0;
			}
			coef_R5[i][j] = temp_R5;
		}
	}
	deleteMatrix(N_y, u_half1); deleteMatrix(N_y, u_half2); deleteMatrix(N_y, u_half3); deleteMatrix(N_y, u_half4); deleteMatrix(N_y, u_half5);
}

void compute_coefy(double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double h_y)
{
	double** u_half1; createMatrix(N_y + 1, N_x, u_half1);
	double** u_half2; createMatrix(N_y + 1, N_x, u_half2);
	double** u_half3; createMatrix(N_y + 1, N_x, u_half3);
	double** u_half4; createMatrix(N_y + 1, N_x, u_half4);
	double** u_half5; createMatrix(N_y + 1, N_x, u_half5);

	compute_vhalf(u_half1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y, Mode);
	compute_vhalf(u_half2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y, Mode);
	compute_vhalf(u_half3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y, Mode);
	compute_vhalf(u_half4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y, Mode);
	compute_vhalf(u_half5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y, Mode);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y; i++)
		{
			double temp_L1 = -1.0 * u_half1[i][j] / 3.0;
			if (temp_L1 < 0)
			{
				temp_L1 = 0.0;

			}
			coef_B1[i][j] = temp_L1;

			double temp_L2 = -1.0 * u_half2[i][j] / 3.0;
			if (temp_L2 < 0)
			{
				temp_L2 = 0.0;

			}
			coef_B2[i][j] = temp_L2;

			double temp_L3 = -1.0 * u_half3[i][j] / 3.0;
			if (temp_L3 < 0)
			{
				temp_L3 = 0.0;

			}
			coef_B3[i][j] = temp_L3;

			double temp_L4 = -1.0 * u_half4[i][j] / 3.0;
			if (temp_L4 < 0)
			{
				temp_L4 = 0.0;

			}
			coef_B4[i][j] = temp_L4;

			double temp_L5 = -1.0 * u_half5[i][j] / 3.0;
			if (temp_L5 < 0)
			{
				temp_L5 = 0.0;

			}
			coef_B5[i][j] = temp_L5;

			double temp_R1 = -1.0 * u_half1[i + 1][j] / 3.0;
			if (temp_R1 > 0)
			{
				temp_R1 = 0.0;
			}
			coef_T1[i][j] = temp_R1;

			double temp_R2 = -1.0 * u_half2[i + 1][j] / 3.0;
			if (temp_R2 > 0)
			{
				temp_R2 = 0.0;
			}
			coef_T2[i][j] = temp_R2;

			double temp_R3 = -1.0 * u_half3[i + 1][j] / 3.0;
			if (temp_R3 > 0)
			{
				temp_R3 = 0.0;
			}
			coef_T3[i][j] = temp_R3;

			double temp_R4 = -1.0 * u_half4[i + 1][j] / 3.0;
			if (temp_R4 > 0)
			{
				temp_R4 = 0.0;
			}
			coef_T4[i][j] = temp_R4;

			double temp_R5 = -1.0 * u_half5[i + 1][j] / 3.0;
			if (temp_R5 > 0)
			{
				temp_R5 = 0.0;
			}
			coef_T5[i][j] = temp_R5;
		}
	}
	deleteMatrix(N_y + 1, u_half1); deleteMatrix(N_y + 1, u_half2); deleteMatrix(N_y + 1, u_half3); deleteMatrix(N_y + 1, u_half4); deleteMatrix(N_y + 1, u_half5);
}

void palpha_jumpx(double** u, double eta, double gamma_alpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** palpha_plus; createMatrix(N_y + 2, N_x + 2, palpha_plus);
	double** palpha_minus; createMatrix(N_y + 2, N_x + 2, palpha_minus);

	p_alpha(palpha_plus, -1.0, eta, gamma_alpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	p_alpha(palpha_minus, 1.0, eta, gamma_alpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_x(u, palpha_plus, palpha_minus);

	deleteMatrix(N_y + 2, palpha_plus); deleteMatrix(N_y + 2, palpha_minus);
}

void palpha_jumpy(double** u, double xi, double gamma_alpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** palpha_plus; createMatrix(N_y + 2, N_x + 2, palpha_plus);
	double** palpha_minus; createMatrix(N_y + 2, N_x + 2, palpha_minus);

	p_alpha(palpha_plus, xi, -1.0, gamma_alpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	p_alpha(palpha_minus, xi, 1.0, gamma_alpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_y(u, palpha_plus, palpha_minus);

	deleteMatrix(N_y + 2, palpha_plus); deleteMatrix(N_y + 2, palpha_minus);
}

void Nonconservetive_Ne(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** F_Ne11; createMatrix(N_y + 2, N_x + 2, F_Ne11); double** F_Ne12; createMatrix(N_y + 2, N_x + 2, F_Ne12);
	double** F_Ne13; createMatrix(N_y + 2, N_x + 2, F_Ne13); double** F_Ne14; createMatrix(N_y + 2, N_x + 2, F_Ne14);
	double** F_Ne15; createMatrix(N_y + 2, N_x + 2, F_Ne15);

	double** F_Ne21; createMatrix(N_y + 2, N_x + 2, F_Ne21); double** F_Ne22; createMatrix(N_y + 2, N_x + 2, F_Ne22);
	double** F_Ne23; createMatrix(N_y + 2, N_x + 2, F_Ne23); double** F_Ne24; createMatrix(N_y + 2, N_x + 2, F_Ne24);
	double** F_Ne25; createMatrix(N_y + 2, N_x + 2, F_Ne25);

	double** F_Ne31; createMatrix(N_y + 2, N_x + 2, F_Ne31); double** F_Ne32; createMatrix(N_y + 2, N_x + 2, F_Ne32);
	double** F_Ne33; createMatrix(N_y + 2, N_x + 2, F_Ne33); double** F_Ne34; createMatrix(N_y + 2, N_x + 2, F_Ne34);
	double** F_Ne35; createMatrix(N_y + 2, N_x + 2, F_Ne35);

	double** F_Ne41; createMatrix(N_y + 2, N_x + 2, F_Ne41); double** F_Ne42; createMatrix(N_y + 2, N_x + 2, F_Ne42);
	double** F_Ne43; createMatrix(N_y + 2, N_x + 2, F_Ne43); double** F_Ne44; createMatrix(N_y + 2, N_x + 2, F_Ne44);
	double** F_Ne45; createMatrix(N_y + 2, N_x + 2, F_Ne45);

	double** F_Ne51; createMatrix(N_y + 2, N_x + 2, F_Ne51); double** F_Ne52; createMatrix(N_y + 2, N_x + 2, F_Ne52);
	double** F_Ne53; createMatrix(N_y + 2, N_x + 2, F_Ne53); double** F_Ne54; createMatrix(N_y + 2, N_x + 2, F_Ne54);
	double** F_Ne55; createMatrix(N_y + 2, N_x + 2, F_Ne55);

	// needed modified
	F_Ne(F_Ne11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ne(F_Ne21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ne(F_Ne31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ne(F_Ne41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ne(F_Ne51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ne(F_Ne55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	double** Ne_GQ_0; createMatrix(N_y + 2, N_x + 2, Ne_GQ_0); double** Ne_GQ_1; createMatrix(N_y + 2, N_x + 2, Ne_GQ_1);
	double** Ne_GQ_2; createMatrix(N_y + 2, N_x + 2, Ne_GQ_2); double** Ne_GQ_3; createMatrix(N_y + 2, N_x + 2, Ne_GQ_3);
	double** Ne_GQ_4; createMatrix(N_y + 2, N_x + 2, Ne_GQ_4); double** Ne_GQ_5; createMatrix(N_y + 2, N_x + 2, Ne_GQ_5);

	get_Gauss_quadrature_2D(Ne_GQ_0, Ne_GQ_1, Ne_GQ_2, Ne_GQ_3, Ne_GQ_4, Ne_GQ_5, F_Ne11, F_Ne21, F_Ne31, F_Ne41, F_Ne51, F_Ne12, F_Ne22, F_Ne32, F_Ne42, F_Ne52, F_Ne13, F_Ne23, F_Ne33, F_Ne43, F_Ne53, F_Ne14, F_Ne24, F_Ne34, F_Ne44, F_Ne54, F_Ne15, F_Ne25, F_Ne35, F_Ne45, F_Ne55, h_x, h_y);

	double** pe_jumpx1; createMatrix(N_y, N_x + 3, pe_jumpx1); double** pe_jumpx2; createMatrix(N_y, N_x + 3, pe_jumpx2);
	double** pe_jumpx3; createMatrix(N_y, N_x + 3, pe_jumpx3); double** pe_jumpx4; createMatrix(N_y, N_x + 3, pe_jumpx4);
	double** pe_jumpx5; createMatrix(N_y, N_x + 3, pe_jumpx5);

	double** pi_jumpx1; createMatrix(N_y, N_x + 3, pi_jumpx1); double** pi_jumpx2; createMatrix(N_y, N_x + 3, pi_jumpx2);
	double** pi_jumpx3; createMatrix(N_y, N_x + 3, pi_jumpx3); double** pi_jumpx4; createMatrix(N_y, N_x + 3, pi_jumpx4);
	double** pi_jumpx5; createMatrix(N_y, N_x + 3, pi_jumpx5);

	double** pr_jumpx1; createMatrix(N_y, N_x + 3, pr_jumpx1); double** pr_jumpx2; createMatrix(N_y, N_x + 3, pr_jumpx2);
	double** pr_jumpx3; createMatrix(N_y, N_x + 3, pr_jumpx3); double** pr_jumpx4; createMatrix(N_y, N_x + 3, pr_jumpx4);
	double** pr_jumpx5; createMatrix(N_y, N_x + 3, pr_jumpx5);

	double** pe_jumpy1; createMatrix(N_y + 3, N_x, pe_jumpy1); double** pe_jumpy2; createMatrix(N_y + 3, N_x, pe_jumpy2);
	double** pe_jumpy3; createMatrix(N_y + 3, N_x, pe_jumpy3); double** pe_jumpy4; createMatrix(N_y + 3, N_x, pe_jumpy4);
	double** pe_jumpy5; createMatrix(N_y + 3, N_x, pe_jumpy5);

	double** pi_jumpy1; createMatrix(N_y + 3, N_x, pi_jumpy1); double** pi_jumpy2; createMatrix(N_y + 3, N_x, pi_jumpy2);
	double** pi_jumpy3; createMatrix(N_y + 3, N_x, pi_jumpy3); double** pi_jumpy4; createMatrix(N_y + 3, N_x, pi_jumpy4);
	double** pi_jumpy5; createMatrix(N_y + 3, N_x, pi_jumpy5);

	double** pr_jumpy1; createMatrix(N_y + 3, N_x, pr_jumpy1); double** pr_jumpy2; createMatrix(N_y + 3, N_x, pr_jumpy2);
	double** pr_jumpy3; createMatrix(N_y + 3, N_x, pr_jumpy3); double** pr_jumpy4; createMatrix(N_y + 3, N_x, pr_jumpy4);
	double** pr_jumpy5; createMatrix(N_y + 3, N_x, pr_jumpy5);

	palpha_jumpx(pe_jumpx1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpx(pi_jumpx1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpx(pr_jumpx1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	palpha_jumpy(pe_jumpy1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpy(pi_jumpy1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpy(pr_jumpy1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	double** GQx_L1; createMatrix(N_y, N_x, GQx_L1); double** GQx_Leta; createMatrix(N_y, N_x, GQx_Leta); double** GQx_Leta2; createMatrix(N_y, N_x, GQx_Leta2);
	double** GQx_R1; createMatrix(N_y, N_x, GQx_R1); double** GQx_Reta; createMatrix(N_y, N_x, GQx_Reta); double** GQx_Reta2; createMatrix(N_y, N_x, GQx_Reta2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			// needed modified
			double tempx_jumpL1 = 2.0 * pe_jumpx1[i][j] - pi_jumpx1[i][j] - pr_jumpx1[i][j];
			double tempx_jumpL2 = 2.0 * pe_jumpx2[i][j] - pi_jumpx2[i][j] - pr_jumpx2[i][j];
			double tempx_jumpL3 = 2.0 * pe_jumpx3[i][j] - pi_jumpx3[i][j] - pr_jumpx3[i][j];
			double tempx_jumpL4 = 2.0 * pe_jumpx4[i][j] - pi_jumpx4[i][j] - pr_jumpx4[i][j];
			double tempx_jumpL5 = 2.0 * pe_jumpx5[i][j] - pi_jumpx5[i][j] - pr_jumpx5[i][j];

			double tempx_jumpR1 = 2.0 * pe_jumpx1[i][j + 1] - pi_jumpx1[i][j + 1] - pr_jumpx1[i][j + 1];
			double tempx_jumpR2 = 2.0 * pe_jumpx2[i][j + 1] - pi_jumpx2[i][j + 1] - pr_jumpx2[i][j + 1];
			double tempx_jumpR3 = 2.0 * pe_jumpx3[i][j + 1] - pi_jumpx3[i][j + 1] - pr_jumpx3[i][j + 1];
			double tempx_jumpR4 = 2.0 * pe_jumpx4[i][j + 1] - pi_jumpx4[i][j + 1] - pr_jumpx4[i][j + 1];
			double tempx_jumpR5 = 2.0 * pe_jumpx5[i][j + 1] - pi_jumpx5[i][j + 1] - pr_jumpx5[i][j + 1];

			double tempx_L1 = coef_L1[i][j - 1]; double tempx_L2 = coef_L2[i][j - 1]; double tempx_L3 = coef_L3[i][j - 1];
			double tempx_L4 = coef_L4[i][j - 1]; double tempx_L5 = coef_L5[i][j - 1];

			double tempx_R1 = coef_R1[i][j - 1]; double tempx_R2 = coef_R2[i][j - 1]; double tempx_R3 = coef_R3[i][j - 1];
			double tempx_R4 = coef_R4[i][j - 1]; double tempx_R5 = coef_R5[i][j - 1];

			GQx_L1[i][j - 1] = weight_1 * tempx_L1 * tempx_jumpL1 + weight_2 * tempx_L2 * tempx_jumpL2 + weight_3 * tempx_L3 * tempx_jumpL3 + weight_4 * tempx_L4 * tempx_jumpL4 + weight_5 * tempx_L5 * tempx_jumpL5;
			GQx_R1[i][j - 1] = weight_1 * tempx_R1 * tempx_jumpR1 + weight_2 * tempx_R2 * tempx_jumpR2 + weight_3 * tempx_R3 * tempx_jumpR3 + weight_4 * tempx_R4 * tempx_jumpR4 + weight_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta[i][j - 1] = weight_1 * GQxi_1 * tempx_L1 * tempx_jumpL1 + weight_2 * GQxi_2 * tempx_L2 * tempx_jumpL2 + weight_3 * GQxi_3 * tempx_L3 * tempx_jumpL3 + weight_4 * GQxi_4 * tempx_L4 * tempx_jumpL4 + weight_5 * GQxi_5 * tempx_L5 * tempx_jumpL5;
			GQx_Reta[i][j - 1] = weight_1 * GQxi_1 * tempx_R1 * tempx_jumpR1 + weight_2 * GQxi_2 * tempx_R2 * tempx_jumpR2 + weight_3 * GQxi_3 * tempx_R3 * tempx_jumpR3 + weight_4 * GQxi_4 * tempx_R4 * tempx_jumpR4 + weight_5 * GQxi_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_L1 * tempx_jumpL1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_L2 * tempx_jumpL2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_L3 * tempx_jumpL3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_L4 * tempx_jumpL4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_L5 * tempx_jumpL5;
			GQx_Reta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_R1 * tempx_jumpR1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_R2 * tempx_jumpR2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_R3 * tempx_jumpR3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_R4 * tempx_jumpR4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_R5 * tempx_jumpR5;

		}
	}

	double** GQy_B1; createMatrix(N_y, N_x, GQy_B1); double** GQy_Bxi; createMatrix(N_y, N_x, GQy_Bxi); double** GQy_Bxi2; createMatrix(N_y, N_x, GQy_Bxi2);
	double** GQy_T1; createMatrix(N_y, N_x, GQy_T1); double** GQy_Txi; createMatrix(N_y, N_x, GQy_Txi); double** GQy_Txi2; createMatrix(N_y, N_x, GQy_Txi2);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 1; i < N_y + 1; i++)
		{
			// needed modified
			double tempy_jumpB1 = 2.0 * pe_jumpy1[i][j] - pi_jumpy1[i][j] - pr_jumpy1[i][j];
			double tempy_jumpB2 = 2.0 * pe_jumpy2[i][j] - pi_jumpy2[i][j] - pr_jumpy2[i][j];
			double tempy_jumpB3 = 2.0 * pe_jumpy3[i][j] - pi_jumpy3[i][j] - pr_jumpy3[i][j];
			double tempy_jumpB4 = 2.0 * pe_jumpy4[i][j] - pi_jumpy4[i][j] - pr_jumpy4[i][j];
			double tempy_jumpB5 = 2.0 * pe_jumpy5[i][j] - pi_jumpy5[i][j] - pr_jumpy5[i][j];

			double tempy_jumpT1 = 2.0 * pe_jumpy1[i + 1][j] - pi_jumpy1[i + 1][j] - pr_jumpy1[i + 1][j];
			double tempy_jumpT2 = 2.0 * pe_jumpy2[i + 1][j] - pi_jumpy2[i + 1][j] - pr_jumpy2[i + 1][j];
			double tempy_jumpT3 = 2.0 * pe_jumpy3[i + 1][j] - pi_jumpy3[i + 1][j] - pr_jumpy3[i + 1][j];
			double tempy_jumpT4 = 2.0 * pe_jumpy4[i + 1][j] - pi_jumpy4[i + 1][j] - pr_jumpy4[i + 1][j];
			double tempy_jumpT5 = 2.0 * pe_jumpy5[i + 1][j] - pi_jumpy5[i + 1][j] - pr_jumpy5[i + 1][j];

			double tempy_B1 = coef_B1[i - 1][j]; double tempy_B2 = coef_B2[i - 1][j]; double tempy_B3 = coef_B3[i - 1][j];
			double tempy_B4 = coef_B4[i - 1][j]; double tempy_B5 = coef_B5[i - 1][j];

			double tempy_T1 = coef_T1[i - 1][j]; double tempy_T2 = coef_T2[i - 1][j]; double tempy_T3 = coef_T3[i - 1][j];
			double tempy_T4 = coef_T4[i - 1][j]; double tempy_T5 = coef_T5[i - 1][j];

			GQy_B1[i - 1][j] = weight_1 * tempy_B1 * tempy_jumpB1 + weight_2 * tempy_B2 * tempy_jumpB2 + weight_3 * tempy_B3 * tempy_jumpB3 + weight_4 * tempy_B4 * tempy_jumpB4 + weight_5 * tempy_B5 * tempy_jumpB5;
			GQy_T1[i - 1][j] = weight_1 * tempy_T1 * tempy_jumpT1 + weight_2 * tempy_T2 * tempy_jumpT2 + weight_3 * tempy_T3 * tempy_jumpT3 + weight_4 * tempy_T4 * tempy_jumpT4 + weight_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi[i - 1][j] = weight_1 * GQxi_1 * tempy_B1 * tempy_jumpB1 + weight_2 * GQxi_2 * tempy_B2 * tempy_jumpB2 + weight_3 * GQxi_3 * tempy_B3 * tempy_jumpB3 + weight_4 * GQxi_4 * tempy_B4 * tempy_jumpB4 + weight_5 * GQxi_5 * tempy_B5 * tempy_jumpB5;
			GQy_Txi[i - 1][j] = weight_1 * GQxi_1 * tempy_T1 * tempy_jumpT1 + weight_2 * GQxi_2 * tempy_T2 * tempy_jumpT2 + weight_3 * GQxi_3 * tempy_T3 * tempy_jumpT3 + weight_4 * GQxi_4 * tempy_T4 * tempy_jumpT4 + weight_5 * GQxi_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_B1 * tempy_jumpB1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_B2 * tempy_jumpB2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_B3 * tempy_jumpB3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_B4 * tempy_jumpB4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_B5 * tempy_jumpB5;
			GQy_Txi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_T1 * tempy_jumpT1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_T2 * tempy_jumpT2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_T3 * tempy_jumpT3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_T4 * tempy_jumpT4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_T5 * tempy_jumpT5;
		}
	}

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i][j] = Ne_GQ_0[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_1[i][j] = Ne_GQ_1[i][j] + h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_2[i][j] = Ne_GQ_2[i][j] - h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_3[i][j] = Ne_GQ_3[i][j] + h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_4[i][j] = Ne_GQ_4[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi2[i - 1][j - 1] - h_x / 2.0 * GQy_Txi2[i - 1][j - 1];
			U_5[i][j] = Ne_GQ_5[i][j] - h_y / 2.0 * GQx_Leta2[i - 1][j - 1] - h_y / 2.0 * GQx_Reta2[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
		}
	}

	deleteMatrix(N_y + 2, F_Ne11); deleteMatrix(N_y + 2, F_Ne12); deleteMatrix(N_y + 2, F_Ne13); deleteMatrix(N_y + 2, F_Ne14); deleteMatrix(N_y + 2, F_Ne15);
	deleteMatrix(N_y + 2, F_Ne21); deleteMatrix(N_y + 2, F_Ne22); deleteMatrix(N_y + 2, F_Ne23); deleteMatrix(N_y + 2, F_Ne24); deleteMatrix(N_y + 2, F_Ne25);
	deleteMatrix(N_y + 2, F_Ne31); deleteMatrix(N_y + 2, F_Ne32); deleteMatrix(N_y + 2, F_Ne33); deleteMatrix(N_y + 2, F_Ne34); deleteMatrix(N_y + 2, F_Ne35);
	deleteMatrix(N_y + 2, F_Ne41); deleteMatrix(N_y + 2, F_Ne42); deleteMatrix(N_y + 2, F_Ne43); deleteMatrix(N_y + 2, F_Ne44); deleteMatrix(N_y + 2, F_Ne45);
	deleteMatrix(N_y + 2, F_Ne51); deleteMatrix(N_y + 2, F_Ne52); deleteMatrix(N_y + 2, F_Ne53); deleteMatrix(N_y + 2, F_Ne54); deleteMatrix(N_y + 2, F_Ne55);

	deleteMatrix(N_y + 2, Ne_GQ_0); deleteMatrix(N_y + 2, Ne_GQ_1); deleteMatrix(N_y + 2, Ne_GQ_2);
	deleteMatrix(N_y + 2, Ne_GQ_3); deleteMatrix(N_y + 2, Ne_GQ_4); deleteMatrix(N_y + 2, Ne_GQ_5);

	deleteMatrix(N_y, pe_jumpx1); deleteMatrix(N_y, pe_jumpx2); deleteMatrix(N_y, pe_jumpx3); deleteMatrix(N_y, pe_jumpx4); deleteMatrix(N_y, pe_jumpx5);
	deleteMatrix(N_y + 3, pe_jumpy1); deleteMatrix(N_y + 3, pe_jumpy2); deleteMatrix(N_y + 3, pe_jumpy3); deleteMatrix(N_y + 3, pe_jumpy4); deleteMatrix(N_y + 3, pe_jumpy5);

	deleteMatrix(N_y, pi_jumpx1); deleteMatrix(N_y, pi_jumpx2); deleteMatrix(N_y, pi_jumpx3); deleteMatrix(N_y, pi_jumpx4); deleteMatrix(N_y, pi_jumpx5);
	deleteMatrix(N_y + 3, pi_jumpy1); deleteMatrix(N_y + 3, pi_jumpy2); deleteMatrix(N_y + 3, pi_jumpy3); deleteMatrix(N_y + 3, pi_jumpy4); deleteMatrix(N_y + 3, pi_jumpy5);

	deleteMatrix(N_y, pr_jumpx1); deleteMatrix(N_y, pr_jumpx2); deleteMatrix(N_y, pr_jumpx3); deleteMatrix(N_y, pr_jumpx4); deleteMatrix(N_y, pr_jumpx5);
	deleteMatrix(N_y + 3, pr_jumpy1); deleteMatrix(N_y + 3, pr_jumpy2); deleteMatrix(N_y + 3, pr_jumpy3); deleteMatrix(N_y + 3, pr_jumpy4); deleteMatrix(N_y + 3, pr_jumpy5);

	deleteMatrix(N_y, GQx_L1); deleteMatrix(N_y, GQx_Leta); deleteMatrix(N_y, GQx_Leta2);
	deleteMatrix(N_y, GQx_R1); deleteMatrix(N_y, GQx_Reta); deleteMatrix(N_y, GQx_Reta2);

	deleteMatrix(N_y, GQy_B1); deleteMatrix(N_y, GQy_Bxi); deleteMatrix(N_y, GQy_Bxi2);
	deleteMatrix(N_y, GQy_T1); deleteMatrix(N_y, GQy_Txi); deleteMatrix(N_y, GQy_Txi2);
}

void Nonconservetive_Ni(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** F_Ne11; createMatrix(N_y + 2, N_x + 2, F_Ne11); double** F_Ne12; createMatrix(N_y + 2, N_x + 2, F_Ne12);
	double** F_Ne13; createMatrix(N_y + 2, N_x + 2, F_Ne13); double** F_Ne14; createMatrix(N_y + 2, N_x + 2, F_Ne14);
	double** F_Ne15; createMatrix(N_y + 2, N_x + 2, F_Ne15);

	double** F_Ne21; createMatrix(N_y + 2, N_x + 2, F_Ne21); double** F_Ne22; createMatrix(N_y + 2, N_x + 2, F_Ne22);
	double** F_Ne23; createMatrix(N_y + 2, N_x + 2, F_Ne23); double** F_Ne24; createMatrix(N_y + 2, N_x + 2, F_Ne24);
	double** F_Ne25; createMatrix(N_y + 2, N_x + 2, F_Ne25);

	double** F_Ne31; createMatrix(N_y + 2, N_x + 2, F_Ne31); double** F_Ne32; createMatrix(N_y + 2, N_x + 2, F_Ne32);
	double** F_Ne33; createMatrix(N_y + 2, N_x + 2, F_Ne33); double** F_Ne34; createMatrix(N_y + 2, N_x + 2, F_Ne34);
	double** F_Ne35; createMatrix(N_y + 2, N_x + 2, F_Ne35);

	double** F_Ne41; createMatrix(N_y + 2, N_x + 2, F_Ne41); double** F_Ne42; createMatrix(N_y + 2, N_x + 2, F_Ne42);
	double** F_Ne43; createMatrix(N_y + 2, N_x + 2, F_Ne43); double** F_Ne44; createMatrix(N_y + 2, N_x + 2, F_Ne44);
	double** F_Ne45; createMatrix(N_y + 2, N_x + 2, F_Ne45);

	double** F_Ne51; createMatrix(N_y + 2, N_x + 2, F_Ne51); double** F_Ne52; createMatrix(N_y + 2, N_x + 2, F_Ne52);
	double** F_Ne53; createMatrix(N_y + 2, N_x + 2, F_Ne53); double** F_Ne54; createMatrix(N_y + 2, N_x + 2, F_Ne54);
	double** F_Ne55; createMatrix(N_y + 2, N_x + 2, F_Ne55);

	// needed modified
	F_Ni(F_Ne11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ni(F_Ne21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ni(F_Ne31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ni(F_Ne41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Ni(F_Ne51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Ni(F_Ne55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	double** Ne_GQ_0; createMatrix(N_y + 2, N_x + 2, Ne_GQ_0); double** Ne_GQ_1; createMatrix(N_y + 2, N_x + 2, Ne_GQ_1);
	double** Ne_GQ_2; createMatrix(N_y + 2, N_x + 2, Ne_GQ_2); double** Ne_GQ_3; createMatrix(N_y + 2, N_x + 2, Ne_GQ_3);
	double** Ne_GQ_4; createMatrix(N_y + 2, N_x + 2, Ne_GQ_4); double** Ne_GQ_5; createMatrix(N_y + 2, N_x + 2, Ne_GQ_5);

	get_Gauss_quadrature_2D(Ne_GQ_0, Ne_GQ_1, Ne_GQ_2, Ne_GQ_3, Ne_GQ_4, Ne_GQ_5, F_Ne11, F_Ne21, F_Ne31, F_Ne41, F_Ne51, F_Ne12, F_Ne22, F_Ne32, F_Ne42, F_Ne52, F_Ne13, F_Ne23, F_Ne33, F_Ne43, F_Ne53, F_Ne14, F_Ne24, F_Ne34, F_Ne44, F_Ne54, F_Ne15, F_Ne25, F_Ne35, F_Ne45, F_Ne55, h_x, h_y);

	double** pe_jumpx1; createMatrix(N_y, N_x + 3, pe_jumpx1); double** pe_jumpx2; createMatrix(N_y, N_x + 3, pe_jumpx2);
	double** pe_jumpx3; createMatrix(N_y, N_x + 3, pe_jumpx3); double** pe_jumpx4; createMatrix(N_y, N_x + 3, pe_jumpx4);
	double** pe_jumpx5; createMatrix(N_y, N_x + 3, pe_jumpx5);

	double** pi_jumpx1; createMatrix(N_y, N_x + 3, pi_jumpx1); double** pi_jumpx2; createMatrix(N_y, N_x + 3, pi_jumpx2);
	double** pi_jumpx3; createMatrix(N_y, N_x + 3, pi_jumpx3); double** pi_jumpx4; createMatrix(N_y, N_x + 3, pi_jumpx4);
	double** pi_jumpx5; createMatrix(N_y, N_x + 3, pi_jumpx5);

	double** pr_jumpx1; createMatrix(N_y, N_x + 3, pr_jumpx1); double** pr_jumpx2; createMatrix(N_y, N_x + 3, pr_jumpx2);
	double** pr_jumpx3; createMatrix(N_y, N_x + 3, pr_jumpx3); double** pr_jumpx4; createMatrix(N_y, N_x + 3, pr_jumpx4);
	double** pr_jumpx5; createMatrix(N_y, N_x + 3, pr_jumpx5);

	double** pe_jumpy1; createMatrix(N_y + 3, N_x, pe_jumpy1); double** pe_jumpy2; createMatrix(N_y + 3, N_x, pe_jumpy2);
	double** pe_jumpy3; createMatrix(N_y + 3, N_x, pe_jumpy3); double** pe_jumpy4; createMatrix(N_y + 3, N_x, pe_jumpy4);
	double** pe_jumpy5; createMatrix(N_y + 3, N_x, pe_jumpy5);

	double** pi_jumpy1; createMatrix(N_y + 3, N_x, pi_jumpy1); double** pi_jumpy2; createMatrix(N_y + 3, N_x, pi_jumpy2);
	double** pi_jumpy3; createMatrix(N_y + 3, N_x, pi_jumpy3); double** pi_jumpy4; createMatrix(N_y + 3, N_x, pi_jumpy4);
	double** pi_jumpy5; createMatrix(N_y + 3, N_x, pi_jumpy5);

	double** pr_jumpy1; createMatrix(N_y + 3, N_x, pr_jumpy1); double** pr_jumpy2; createMatrix(N_y + 3, N_x, pr_jumpy2);
	double** pr_jumpy3; createMatrix(N_y + 3, N_x, pr_jumpy3); double** pr_jumpy4; createMatrix(N_y + 3, N_x, pr_jumpy4);
	double** pr_jumpy5; createMatrix(N_y + 3, N_x, pr_jumpy5);

	palpha_jumpx(pe_jumpx1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpx(pi_jumpx1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpx(pr_jumpx1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	palpha_jumpy(pe_jumpy1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpy(pi_jumpy1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpy(pr_jumpy1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	double** GQx_L1; createMatrix(N_y, N_x, GQx_L1); double** GQx_Leta; createMatrix(N_y, N_x, GQx_Leta); double** GQx_Leta2; createMatrix(N_y, N_x, GQx_Leta2);
	double** GQx_R1; createMatrix(N_y, N_x, GQx_R1); double** GQx_Reta; createMatrix(N_y, N_x, GQx_Reta); double** GQx_Reta2; createMatrix(N_y, N_x, GQx_Reta2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			// needed modified
			double tempx_jumpL1 = 2.0 * pi_jumpx1[i][j] - pe_jumpx1[i][j] - pr_jumpx1[i][j];
			double tempx_jumpL2 = 2.0 * pi_jumpx2[i][j] - pe_jumpx2[i][j] - pr_jumpx2[i][j];
			double tempx_jumpL3 = 2.0 * pi_jumpx3[i][j] - pe_jumpx3[i][j] - pr_jumpx3[i][j];
			double tempx_jumpL4 = 2.0 * pi_jumpx4[i][j] - pe_jumpx4[i][j] - pr_jumpx4[i][j];
			double tempx_jumpL5 = 2.0 * pi_jumpx5[i][j] - pe_jumpx5[i][j] - pr_jumpx5[i][j];

			double tempx_jumpR1 = 2.0 * pi_jumpx1[i][j + 1] - pe_jumpx1[i][j + 1] - pr_jumpx1[i][j + 1];
			double tempx_jumpR2 = 2.0 * pi_jumpx2[i][j + 1] - pe_jumpx2[i][j + 1] - pr_jumpx2[i][j + 1];
			double tempx_jumpR3 = 2.0 * pi_jumpx3[i][j + 1] - pe_jumpx3[i][j + 1] - pr_jumpx3[i][j + 1];
			double tempx_jumpR4 = 2.0 * pi_jumpx4[i][j + 1] - pe_jumpx4[i][j + 1] - pr_jumpx4[i][j + 1];
			double tempx_jumpR5 = 2.0 * pi_jumpx5[i][j + 1] - pe_jumpx5[i][j + 1] - pr_jumpx5[i][j + 1];

			double tempx_L1 = coef_L1[i][j - 1]; double tempx_L2 = coef_L2[i][j - 1]; double tempx_L3 = coef_L3[i][j - 1];
			double tempx_L4 = coef_L4[i][j - 1]; double tempx_L5 = coef_L5[i][j - 1];

			double tempx_R1 = coef_R1[i][j - 1]; double tempx_R2 = coef_R2[i][j - 1]; double tempx_R3 = coef_R3[i][j - 1];
			double tempx_R4 = coef_R4[i][j - 1]; double tempx_R5 = coef_R5[i][j - 1];

			GQx_L1[i][j - 1] = weight_1 * tempx_L1 * tempx_jumpL1 + weight_2 * tempx_L2 * tempx_jumpL2 + weight_3 * tempx_L3 * tempx_jumpL3 + weight_4 * tempx_L4 * tempx_jumpL4 + weight_5 * tempx_L5 * tempx_jumpL5;
			GQx_R1[i][j - 1] = weight_1 * tempx_R1 * tempx_jumpR1 + weight_2 * tempx_R2 * tempx_jumpR2 + weight_3 * tempx_R3 * tempx_jumpR3 + weight_4 * tempx_R4 * tempx_jumpR4 + weight_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta[i][j - 1] = weight_1 * GQxi_1 * tempx_L1 * tempx_jumpL1 + weight_2 * GQxi_2 * tempx_L2 * tempx_jumpL2 + weight_3 * GQxi_3 * tempx_L3 * tempx_jumpL3 + weight_4 * GQxi_4 * tempx_L4 * tempx_jumpL4 + weight_5 * GQxi_5 * tempx_L5 * tempx_jumpL5;
			GQx_Reta[i][j - 1] = weight_1 * GQxi_1 * tempx_R1 * tempx_jumpR1 + weight_2 * GQxi_2 * tempx_R2 * tempx_jumpR2 + weight_3 * GQxi_3 * tempx_R3 * tempx_jumpR3 + weight_4 * GQxi_4 * tempx_R4 * tempx_jumpR4 + weight_5 * GQxi_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_L1 * tempx_jumpL1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_L2 * tempx_jumpL2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_L3 * tempx_jumpL3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_L4 * tempx_jumpL4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_L5 * tempx_jumpL5;
			GQx_Reta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_R1 * tempx_jumpR1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_R2 * tempx_jumpR2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_R3 * tempx_jumpR3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_R4 * tempx_jumpR4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_R5 * tempx_jumpR5;

		}
	}

	double** GQy_B1; createMatrix(N_y, N_x, GQy_B1); double** GQy_Bxi; createMatrix(N_y, N_x, GQy_Bxi); double** GQy_Bxi2; createMatrix(N_y, N_x, GQy_Bxi2);
	double** GQy_T1; createMatrix(N_y, N_x, GQy_T1); double** GQy_Txi; createMatrix(N_y, N_x, GQy_Txi); double** GQy_Txi2; createMatrix(N_y, N_x, GQy_Txi2);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 1; i < N_y + 1; i++)
		{
			// needed modified
			double tempy_jumpB1 = 2.0 * pi_jumpy1[i][j] - pe_jumpy1[i][j] - pr_jumpy1[i][j];
			double tempy_jumpB2 = 2.0 * pi_jumpy2[i][j] - pe_jumpy2[i][j] - pr_jumpy2[i][j];
			double tempy_jumpB3 = 2.0 * pi_jumpy3[i][j] - pe_jumpy3[i][j] - pr_jumpy3[i][j];
			double tempy_jumpB4 = 2.0 * pi_jumpy4[i][j] - pe_jumpy4[i][j] - pr_jumpy4[i][j];
			double tempy_jumpB5 = 2.0 * pi_jumpy5[i][j] - pe_jumpy5[i][j] - pr_jumpy5[i][j];

			double tempy_jumpT1 = 2.0 * pi_jumpy1[i + 1][j] - pe_jumpy1[i + 1][j] - pr_jumpy1[i + 1][j];
			double tempy_jumpT2 = 2.0 * pi_jumpy2[i + 1][j] - pe_jumpy2[i + 1][j] - pr_jumpy2[i + 1][j];
			double tempy_jumpT3 = 2.0 * pi_jumpy3[i + 1][j] - pe_jumpy3[i + 1][j] - pr_jumpy3[i + 1][j];
			double tempy_jumpT4 = 2.0 * pi_jumpy4[i + 1][j] - pe_jumpy4[i + 1][j] - pr_jumpy4[i + 1][j];
			double tempy_jumpT5 = 2.0 * pi_jumpy5[i + 1][j] - pe_jumpy5[i + 1][j] - pr_jumpy5[i + 1][j];

			double tempy_B1 = coef_B1[i - 1][j]; double tempy_B2 = coef_B2[i - 1][j]; double tempy_B3 = coef_B3[i - 1][j];
			double tempy_B4 = coef_B4[i - 1][j]; double tempy_B5 = coef_B5[i - 1][j];

			double tempy_T1 = coef_T1[i - 1][j]; double tempy_T2 = coef_T2[i - 1][j]; double tempy_T3 = coef_T3[i - 1][j];
			double tempy_T4 = coef_T4[i - 1][j]; double tempy_T5 = coef_T5[i - 1][j];

			GQy_B1[i - 1][j] = weight_1 * tempy_B1 * tempy_jumpB1 + weight_2 * tempy_B2 * tempy_jumpB2 + weight_3 * tempy_B3 * tempy_jumpB3 + weight_4 * tempy_B4 * tempy_jumpB4 + weight_5 * tempy_B5 * tempy_jumpB5;
			GQy_T1[i - 1][j] = weight_1 * tempy_T1 * tempy_jumpT1 + weight_2 * tempy_T2 * tempy_jumpT2 + weight_3 * tempy_T3 * tempy_jumpT3 + weight_4 * tempy_T4 * tempy_jumpT4 + weight_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi[i - 1][j] = weight_1 * GQxi_1 * tempy_B1 * tempy_jumpB1 + weight_2 * GQxi_2 * tempy_B2 * tempy_jumpB2 + weight_3 * GQxi_3 * tempy_B3 * tempy_jumpB3 + weight_4 * GQxi_4 * tempy_B4 * tempy_jumpB4 + weight_5 * GQxi_5 * tempy_B5 * tempy_jumpB5;
			GQy_Txi[i - 1][j] = weight_1 * GQxi_1 * tempy_T1 * tempy_jumpT1 + weight_2 * GQxi_2 * tempy_T2 * tempy_jumpT2 + weight_3 * GQxi_3 * tempy_T3 * tempy_jumpT3 + weight_4 * GQxi_4 * tempy_T4 * tempy_jumpT4 + weight_5 * GQxi_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_B1 * tempy_jumpB1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_B2 * tempy_jumpB2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_B3 * tempy_jumpB3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_B4 * tempy_jumpB4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_B5 * tempy_jumpB5;
			GQy_Txi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_T1 * tempy_jumpT1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_T2 * tempy_jumpT2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_T3 * tempy_jumpT3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_T4 * tempy_jumpT4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_T5 * tempy_jumpT5;
		}
	}

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i][j] = Ne_GQ_0[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_1[i][j] = Ne_GQ_1[i][j] + h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_2[i][j] = Ne_GQ_2[i][j] - h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_3[i][j] = Ne_GQ_3[i][j] + h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_4[i][j] = Ne_GQ_4[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi2[i - 1][j - 1] - h_x / 2.0 * GQy_Txi2[i - 1][j - 1];
			U_5[i][j] = Ne_GQ_5[i][j] - h_y / 2.0 * GQx_Leta2[i - 1][j - 1] - h_y / 2.0 * GQx_Reta2[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
		}
	}

	deleteMatrix(N_y + 2, F_Ne11); deleteMatrix(N_y + 2, F_Ne12); deleteMatrix(N_y + 2, F_Ne13); deleteMatrix(N_y + 2, F_Ne14); deleteMatrix(N_y + 2, F_Ne15);
	deleteMatrix(N_y + 2, F_Ne21); deleteMatrix(N_y + 2, F_Ne22); deleteMatrix(N_y + 2, F_Ne23); deleteMatrix(N_y + 2, F_Ne24); deleteMatrix(N_y + 2, F_Ne25);
	deleteMatrix(N_y + 2, F_Ne31); deleteMatrix(N_y + 2, F_Ne32); deleteMatrix(N_y + 2, F_Ne33); deleteMatrix(N_y + 2, F_Ne34); deleteMatrix(N_y + 2, F_Ne35);
	deleteMatrix(N_y + 2, F_Ne41); deleteMatrix(N_y + 2, F_Ne42); deleteMatrix(N_y + 2, F_Ne43); deleteMatrix(N_y + 2, F_Ne44); deleteMatrix(N_y + 2, F_Ne45);
	deleteMatrix(N_y + 2, F_Ne51); deleteMatrix(N_y + 2, F_Ne52); deleteMatrix(N_y + 2, F_Ne53); deleteMatrix(N_y + 2, F_Ne54); deleteMatrix(N_y + 2, F_Ne55);

	deleteMatrix(N_y + 2, Ne_GQ_0); deleteMatrix(N_y + 2, Ne_GQ_1); deleteMatrix(N_y + 2, Ne_GQ_2);
	deleteMatrix(N_y + 2, Ne_GQ_3); deleteMatrix(N_y + 2, Ne_GQ_4); deleteMatrix(N_y + 2, Ne_GQ_5);

	deleteMatrix(N_y, pe_jumpx1); deleteMatrix(N_y, pe_jumpx2); deleteMatrix(N_y, pe_jumpx3); deleteMatrix(N_y, pe_jumpx4); deleteMatrix(N_y, pe_jumpx5);
	deleteMatrix(N_y + 3, pe_jumpy1); deleteMatrix(N_y + 3, pe_jumpy2); deleteMatrix(N_y + 3, pe_jumpy3); deleteMatrix(N_y + 3, pe_jumpy4); deleteMatrix(N_y + 3, pe_jumpy5);

	deleteMatrix(N_y, pi_jumpx1); deleteMatrix(N_y, pi_jumpx2); deleteMatrix(N_y, pi_jumpx3); deleteMatrix(N_y, pi_jumpx4); deleteMatrix(N_y, pi_jumpx5);
	deleteMatrix(N_y + 3, pi_jumpy1); deleteMatrix(N_y + 3, pi_jumpy2); deleteMatrix(N_y + 3, pi_jumpy3); deleteMatrix(N_y + 3, pi_jumpy4); deleteMatrix(N_y + 3, pi_jumpy5);

	deleteMatrix(N_y, pr_jumpx1); deleteMatrix(N_y, pr_jumpx2); deleteMatrix(N_y, pr_jumpx3); deleteMatrix(N_y, pr_jumpx4); deleteMatrix(N_y, pr_jumpx5);
	deleteMatrix(N_y + 3, pr_jumpy1); deleteMatrix(N_y + 3, pr_jumpy2); deleteMatrix(N_y + 3, pr_jumpy3); deleteMatrix(N_y + 3, pr_jumpy4); deleteMatrix(N_y + 3, pr_jumpy5);

	deleteMatrix(N_y, GQx_L1); deleteMatrix(N_y, GQx_Leta); deleteMatrix(N_y, GQx_Leta2);
	deleteMatrix(N_y, GQx_R1); deleteMatrix(N_y, GQx_Reta); deleteMatrix(N_y, GQx_Reta2);

	deleteMatrix(N_y, GQy_B1); deleteMatrix(N_y, GQy_Bxi); deleteMatrix(N_y, GQy_Bxi2);
	deleteMatrix(N_y, GQy_T1); deleteMatrix(N_y, GQy_Txi); deleteMatrix(N_y, GQy_Txi2);
}

void Nonconservetive_Nr(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** F_Ne11; createMatrix(N_y + 2, N_x + 2, F_Ne11); double** F_Ne12; createMatrix(N_y + 2, N_x + 2, F_Ne12);
	double** F_Ne13; createMatrix(N_y + 2, N_x + 2, F_Ne13); double** F_Ne14; createMatrix(N_y + 2, N_x + 2, F_Ne14);
	double** F_Ne15; createMatrix(N_y + 2, N_x + 2, F_Ne15);

	double** F_Ne21; createMatrix(N_y + 2, N_x + 2, F_Ne21); double** F_Ne22; createMatrix(N_y + 2, N_x + 2, F_Ne22);
	double** F_Ne23; createMatrix(N_y + 2, N_x + 2, F_Ne23); double** F_Ne24; createMatrix(N_y + 2, N_x + 2, F_Ne24);
	double** F_Ne25; createMatrix(N_y + 2, N_x + 2, F_Ne25);

	double** F_Ne31; createMatrix(N_y + 2, N_x + 2, F_Ne31); double** F_Ne32; createMatrix(N_y + 2, N_x + 2, F_Ne32);
	double** F_Ne33; createMatrix(N_y + 2, N_x + 2, F_Ne33); double** F_Ne34; createMatrix(N_y + 2, N_x + 2, F_Ne34);
	double** F_Ne35; createMatrix(N_y + 2, N_x + 2, F_Ne35);

	double** F_Ne41; createMatrix(N_y + 2, N_x + 2, F_Ne41); double** F_Ne42; createMatrix(N_y + 2, N_x + 2, F_Ne42);
	double** F_Ne43; createMatrix(N_y + 2, N_x + 2, F_Ne43); double** F_Ne44; createMatrix(N_y + 2, N_x + 2, F_Ne44);
	double** F_Ne45; createMatrix(N_y + 2, N_x + 2, F_Ne45);

	double** F_Ne51; createMatrix(N_y + 2, N_x + 2, F_Ne51); double** F_Ne52; createMatrix(N_y + 2, N_x + 2, F_Ne52);
	double** F_Ne53; createMatrix(N_y + 2, N_x + 2, F_Ne53); double** F_Ne54; createMatrix(N_y + 2, N_x + 2, F_Ne54);
	double** F_Ne55; createMatrix(N_y + 2, N_x + 2, F_Ne55);

	// needed modified
	F_Nr(F_Ne11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Nr(F_Ne21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Nr(F_Ne31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Nr(F_Ne41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	F_Nr(F_Ne51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	F_Nr(F_Ne55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	double** Ne_GQ_0; createMatrix(N_y + 2, N_x + 2, Ne_GQ_0); double** Ne_GQ_1; createMatrix(N_y + 2, N_x + 2, Ne_GQ_1);
	double** Ne_GQ_2; createMatrix(N_y + 2, N_x + 2, Ne_GQ_2); double** Ne_GQ_3; createMatrix(N_y + 2, N_x + 2, Ne_GQ_3);
	double** Ne_GQ_4; createMatrix(N_y + 2, N_x + 2, Ne_GQ_4); double** Ne_GQ_5; createMatrix(N_y + 2, N_x + 2, Ne_GQ_5);

	get_Gauss_quadrature_2D(Ne_GQ_0, Ne_GQ_1, Ne_GQ_2, Ne_GQ_3, Ne_GQ_4, Ne_GQ_5, F_Ne11, F_Ne21, F_Ne31, F_Ne41, F_Ne51, F_Ne12, F_Ne22, F_Ne32, F_Ne42, F_Ne52, F_Ne13, F_Ne23, F_Ne33, F_Ne43, F_Ne53, F_Ne14, F_Ne24, F_Ne34, F_Ne44, F_Ne54, F_Ne15, F_Ne25, F_Ne35, F_Ne45, F_Ne55, h_x, h_y);

	double** pe_jumpx1; createMatrix(N_y, N_x + 3, pe_jumpx1); double** pe_jumpx2; createMatrix(N_y, N_x + 3, pe_jumpx2);
	double** pe_jumpx3; createMatrix(N_y, N_x + 3, pe_jumpx3); double** pe_jumpx4; createMatrix(N_y, N_x + 3, pe_jumpx4);
	double** pe_jumpx5; createMatrix(N_y, N_x + 3, pe_jumpx5);

	double** pi_jumpx1; createMatrix(N_y, N_x + 3, pi_jumpx1); double** pi_jumpx2; createMatrix(N_y, N_x + 3, pi_jumpx2);
	double** pi_jumpx3; createMatrix(N_y, N_x + 3, pi_jumpx3); double** pi_jumpx4; createMatrix(N_y, N_x + 3, pi_jumpx4);
	double** pi_jumpx5; createMatrix(N_y, N_x + 3, pi_jumpx5);

	double** pr_jumpx1; createMatrix(N_y, N_x + 3, pr_jumpx1); double** pr_jumpx2; createMatrix(N_y, N_x + 3, pr_jumpx2);
	double** pr_jumpx3; createMatrix(N_y, N_x + 3, pr_jumpx3); double** pr_jumpx4; createMatrix(N_y, N_x + 3, pr_jumpx4);
	double** pr_jumpx5; createMatrix(N_y, N_x + 3, pr_jumpx5);

	double** pe_jumpy1; createMatrix(N_y + 3, N_x, pe_jumpy1); double** pe_jumpy2; createMatrix(N_y + 3, N_x, pe_jumpy2);
	double** pe_jumpy3; createMatrix(N_y + 3, N_x, pe_jumpy3); double** pe_jumpy4; createMatrix(N_y + 3, N_x, pe_jumpy4);
	double** pe_jumpy5; createMatrix(N_y + 3, N_x, pe_jumpy5);

	double** pi_jumpy1; createMatrix(N_y + 3, N_x, pi_jumpy1); double** pi_jumpy2; createMatrix(N_y + 3, N_x, pi_jumpy2);
	double** pi_jumpy3; createMatrix(N_y + 3, N_x, pi_jumpy3); double** pi_jumpy4; createMatrix(N_y + 3, N_x, pi_jumpy4);
	double** pi_jumpy5; createMatrix(N_y + 3, N_x, pi_jumpy5);

	double** pr_jumpy1; createMatrix(N_y + 3, N_x, pr_jumpy1); double** pr_jumpy2; createMatrix(N_y + 3, N_x, pr_jumpy2);
	double** pr_jumpy3; createMatrix(N_y + 3, N_x, pr_jumpy3); double** pr_jumpy4; createMatrix(N_y + 3, N_x, pr_jumpy4);
	double** pr_jumpy5; createMatrix(N_y + 3, N_x, pr_jumpy5);

	palpha_jumpx(pe_jumpx1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpx(pe_jumpx5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpx(pi_jumpx1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpx(pi_jumpx5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpx(pr_jumpx1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpx(pr_jumpx5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	palpha_jumpy(pe_jumpy1, GQxi_1, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy2, GQxi_2, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy3, GQxi_3, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy4, GQxi_4, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	palpha_jumpy(pe_jumpy5, GQxi_5, gamma_e, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);

	palpha_jumpy(pi_jumpy1, GQxi_1, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy2, GQxi_2, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy3, GQxi_3, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy4, GQxi_4, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	palpha_jumpy(pi_jumpy5, GQxi_5, gamma_i, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);

	palpha_jumpy(pr_jumpy1, GQxi_1, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy2, GQxi_2, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy3, GQxi_3, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy4, GQxi_4, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	palpha_jumpy(pr_jumpy5, GQxi_5, gamma_r, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	double** GQx_L1; createMatrix(N_y, N_x, GQx_L1); double** GQx_Leta; createMatrix(N_y, N_x, GQx_Leta); double** GQx_Leta2; createMatrix(N_y, N_x, GQx_Leta2);
	double** GQx_R1; createMatrix(N_y, N_x, GQx_R1); double** GQx_Reta; createMatrix(N_y, N_x, GQx_Reta); double** GQx_Reta2; createMatrix(N_y, N_x, GQx_Reta2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			// needed modified
			double tempx_jumpL1 = 2.0 * pr_jumpx1[i][j] - pi_jumpx1[i][j] - pe_jumpx1[i][j];
			double tempx_jumpL2 = 2.0 * pr_jumpx2[i][j] - pi_jumpx2[i][j] - pe_jumpx2[i][j];
			double tempx_jumpL3 = 2.0 * pr_jumpx3[i][j] - pi_jumpx3[i][j] - pe_jumpx3[i][j];
			double tempx_jumpL4 = 2.0 * pr_jumpx4[i][j] - pi_jumpx4[i][j] - pe_jumpx4[i][j];
			double tempx_jumpL5 = 2.0 * pr_jumpx5[i][j] - pi_jumpx5[i][j] - pe_jumpx5[i][j];

			double tempx_jumpR1 = 2.0 * pr_jumpx1[i][j + 1] - pi_jumpx1[i][j + 1] - pe_jumpx1[i][j + 1];
			double tempx_jumpR2 = 2.0 * pr_jumpx2[i][j + 1] - pi_jumpx2[i][j + 1] - pe_jumpx2[i][j + 1];
			double tempx_jumpR3 = 2.0 * pr_jumpx3[i][j + 1] - pi_jumpx3[i][j + 1] - pe_jumpx3[i][j + 1];
			double tempx_jumpR4 = 2.0 * pr_jumpx4[i][j + 1] - pi_jumpx4[i][j + 1] - pe_jumpx4[i][j + 1];
			double tempx_jumpR5 = 2.0 * pr_jumpx5[i][j + 1] - pi_jumpx5[i][j + 1] - pe_jumpx5[i][j + 1];

			double tempx_L1 = coef_L1[i][j - 1]; double tempx_L2 = coef_L2[i][j - 1]; double tempx_L3 = coef_L3[i][j - 1];
			double tempx_L4 = coef_L4[i][j - 1]; double tempx_L5 = coef_L5[i][j - 1];

			double tempx_R1 = coef_R1[i][j - 1]; double tempx_R2 = coef_R2[i][j - 1]; double tempx_R3 = coef_R3[i][j - 1];
			double tempx_R4 = coef_R4[i][j - 1]; double tempx_R5 = coef_R5[i][j - 1];

			GQx_L1[i][j - 1] = weight_1 * tempx_L1 * tempx_jumpL1 + weight_2 * tempx_L2 * tempx_jumpL2 + weight_3 * tempx_L3 * tempx_jumpL3 + weight_4 * tempx_L4 * tempx_jumpL4 + weight_5 * tempx_L5 * tempx_jumpL5;
			GQx_R1[i][j - 1] = weight_1 * tempx_R1 * tempx_jumpR1 + weight_2 * tempx_R2 * tempx_jumpR2 + weight_3 * tempx_R3 * tempx_jumpR3 + weight_4 * tempx_R4 * tempx_jumpR4 + weight_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta[i][j - 1] = weight_1 * GQxi_1 * tempx_L1 * tempx_jumpL1 + weight_2 * GQxi_2 * tempx_L2 * tempx_jumpL2 + weight_3 * GQxi_3 * tempx_L3 * tempx_jumpL3 + weight_4 * GQxi_4 * tempx_L4 * tempx_jumpL4 + weight_5 * GQxi_5 * tempx_L5 * tempx_jumpL5;
			GQx_Reta[i][j - 1] = weight_1 * GQxi_1 * tempx_R1 * tempx_jumpR1 + weight_2 * GQxi_2 * tempx_R2 * tempx_jumpR2 + weight_3 * GQxi_3 * tempx_R3 * tempx_jumpR3 + weight_4 * GQxi_4 * tempx_R4 * tempx_jumpR4 + weight_5 * GQxi_5 * tempx_R5 * tempx_jumpR5;

			GQx_Leta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_L1 * tempx_jumpL1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_L2 * tempx_jumpL2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_L3 * tempx_jumpL3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_L4 * tempx_jumpL4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_L5 * tempx_jumpL5;
			GQx_Reta2[i][j - 1] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempx_R1 * tempx_jumpR1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempx_R2 * tempx_jumpR2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempx_R3 * tempx_jumpR3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempx_R4 * tempx_jumpR4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempx_R5 * tempx_jumpR5;

		}
	}

	double** GQy_B1; createMatrix(N_y, N_x, GQy_B1); double** GQy_Bxi; createMatrix(N_y, N_x, GQy_Bxi); double** GQy_Bxi2; createMatrix(N_y, N_x, GQy_Bxi2);
	double** GQy_T1; createMatrix(N_y, N_x, GQy_T1); double** GQy_Txi; createMatrix(N_y, N_x, GQy_Txi); double** GQy_Txi2; createMatrix(N_y, N_x, GQy_Txi2);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 1; i < N_y + 1; i++)
		{
			// needed modified
			double tempy_jumpB1 = 2.0 * pr_jumpy1[i][j] - pi_jumpy1[i][j] - pe_jumpy1[i][j];
			double tempy_jumpB2 = 2.0 * pr_jumpy2[i][j] - pi_jumpy2[i][j] - pe_jumpy2[i][j];
			double tempy_jumpB3 = 2.0 * pr_jumpy3[i][j] - pi_jumpy3[i][j] - pe_jumpy3[i][j];
			double tempy_jumpB4 = 2.0 * pr_jumpy4[i][j] - pi_jumpy4[i][j] - pe_jumpy4[i][j];
			double tempy_jumpB5 = 2.0 * pr_jumpy5[i][j] - pi_jumpy5[i][j] - pe_jumpy5[i][j];

			double tempy_jumpT1 = 2.0 * pr_jumpy1[i + 1][j] - pi_jumpy1[i + 1][j] - pe_jumpy1[i + 1][j];
			double tempy_jumpT2 = 2.0 * pr_jumpy2[i + 1][j] - pi_jumpy2[i + 1][j] - pe_jumpy2[i + 1][j];
			double tempy_jumpT3 = 2.0 * pr_jumpy3[i + 1][j] - pi_jumpy3[i + 1][j] - pe_jumpy3[i + 1][j];
			double tempy_jumpT4 = 2.0 * pr_jumpy4[i + 1][j] - pi_jumpy4[i + 1][j] - pe_jumpy4[i + 1][j];
			double tempy_jumpT5 = 2.0 * pr_jumpy5[i + 1][j] - pi_jumpy5[i + 1][j] - pe_jumpy5[i + 1][j];

			double tempy_B1 = coef_B1[i - 1][j]; double tempy_B2 = coef_B2[i - 1][j]; double tempy_B3 = coef_B3[i - 1][j];
			double tempy_B4 = coef_B4[i - 1][j]; double tempy_B5 = coef_B5[i - 1][j];

			double tempy_T1 = coef_T1[i - 1][j]; double tempy_T2 = coef_T2[i - 1][j]; double tempy_T3 = coef_T3[i - 1][j];
			double tempy_T4 = coef_T4[i - 1][j]; double tempy_T5 = coef_T5[i - 1][j];

			GQy_B1[i - 1][j] = weight_1 * tempy_B1 * tempy_jumpB1 + weight_2 * tempy_B2 * tempy_jumpB2 + weight_3 * tempy_B3 * tempy_jumpB3 + weight_4 * tempy_B4 * tempy_jumpB4 + weight_5 * tempy_B5 * tempy_jumpB5;
			GQy_T1[i - 1][j] = weight_1 * tempy_T1 * tempy_jumpT1 + weight_2 * tempy_T2 * tempy_jumpT2 + weight_3 * tempy_T3 * tempy_jumpT3 + weight_4 * tempy_T4 * tempy_jumpT4 + weight_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi[i - 1][j] = weight_1 * GQxi_1 * tempy_B1 * tempy_jumpB1 + weight_2 * GQxi_2 * tempy_B2 * tempy_jumpB2 + weight_3 * GQxi_3 * tempy_B3 * tempy_jumpB3 + weight_4 * GQxi_4 * tempy_B4 * tempy_jumpB4 + weight_5 * GQxi_5 * tempy_B5 * tempy_jumpB5;
			GQy_Txi[i - 1][j] = weight_1 * GQxi_1 * tempy_T1 * tempy_jumpT1 + weight_2 * GQxi_2 * tempy_T2 * tempy_jumpT2 + weight_3 * GQxi_3 * tempy_T3 * tempy_jumpT3 + weight_4 * GQxi_4 * tempy_T4 * tempy_jumpT4 + weight_5 * GQxi_5 * tempy_T5 * tempy_jumpT5;

			GQy_Bxi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_B1 * tempy_jumpB1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_B2 * tempy_jumpB2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_B3 * tempy_jumpB3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_B4 * tempy_jumpB4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_B5 * tempy_jumpB5;
			GQy_Txi2[i - 1][j] = weight_1 * (3.0 * GQxi_1 * GQxi_1 - 1.0) / 2.0 * tempy_T1 * tempy_jumpT1 + weight_2 * (3.0 * GQxi_2 * GQxi_2 - 1.0) / 2.0 * tempy_T2 * tempy_jumpT2 + weight_3 * (3.0 * GQxi_3 * GQxi_3 - 1.0) / 2.0 * tempy_T3 * tempy_jumpT3 + weight_4 * (3.0 * GQxi_4 * GQxi_4 - 1.0) / 2.0 * tempy_T4 * tempy_jumpT4 + weight_5 * (3.0 * GQxi_5 * GQxi_5 - 1.0) / 2.0 * tempy_T5 * tempy_jumpT5;
		}
	}

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i][j] = Ne_GQ_0[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_1[i][j] = Ne_GQ_1[i][j] + h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_2[i][j] = Ne_GQ_2[i][j] - h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
			U_3[i][j] = Ne_GQ_3[i][j] + h_y / 2.0 * GQx_Leta[i - 1][j - 1] - h_y / 2.0 * GQx_Reta[i - 1][j - 1] + h_x / 2.0 * GQy_Bxi[i - 1][j - 1] - h_x / 2.0 * GQy_Txi[i - 1][j - 1];
			U_4[i][j] = Ne_GQ_4[i][j] - h_y / 2.0 * GQx_L1[i - 1][j - 1] - h_y / 2.0 * GQx_R1[i - 1][j - 1] - h_x / 2.0 * GQy_Bxi2[i - 1][j - 1] - h_x / 2.0 * GQy_Txi2[i - 1][j - 1];
			U_5[i][j] = Ne_GQ_5[i][j] - h_y / 2.0 * GQx_Leta2[i - 1][j - 1] - h_y / 2.0 * GQx_Reta2[i - 1][j - 1] - h_x / 2.0 * GQy_B1[i - 1][j - 1] - h_x / 2.0 * GQy_T1[i - 1][j - 1];
		}
	}

	deleteMatrix(N_y + 2, F_Ne11); deleteMatrix(N_y + 2, F_Ne12); deleteMatrix(N_y + 2, F_Ne13); deleteMatrix(N_y + 2, F_Ne14); deleteMatrix(N_y + 2, F_Ne15);
	deleteMatrix(N_y + 2, F_Ne21); deleteMatrix(N_y + 2, F_Ne22); deleteMatrix(N_y + 2, F_Ne23); deleteMatrix(N_y + 2, F_Ne24); deleteMatrix(N_y + 2, F_Ne25);
	deleteMatrix(N_y + 2, F_Ne31); deleteMatrix(N_y + 2, F_Ne32); deleteMatrix(N_y + 2, F_Ne33); deleteMatrix(N_y + 2, F_Ne34); deleteMatrix(N_y + 2, F_Ne35);
	deleteMatrix(N_y + 2, F_Ne41); deleteMatrix(N_y + 2, F_Ne42); deleteMatrix(N_y + 2, F_Ne43); deleteMatrix(N_y + 2, F_Ne44); deleteMatrix(N_y + 2, F_Ne45);
	deleteMatrix(N_y + 2, F_Ne51); deleteMatrix(N_y + 2, F_Ne52); deleteMatrix(N_y + 2, F_Ne53); deleteMatrix(N_y + 2, F_Ne54); deleteMatrix(N_y + 2, F_Ne55);

	deleteMatrix(N_y + 2, Ne_GQ_0); deleteMatrix(N_y + 2, Ne_GQ_1); deleteMatrix(N_y + 2, Ne_GQ_2);
	deleteMatrix(N_y + 2, Ne_GQ_3); deleteMatrix(N_y + 2, Ne_GQ_4); deleteMatrix(N_y + 2, Ne_GQ_5);

	deleteMatrix(N_y, pe_jumpx1); deleteMatrix(N_y, pe_jumpx2); deleteMatrix(N_y, pe_jumpx3); deleteMatrix(N_y, pe_jumpx4); deleteMatrix(N_y, pe_jumpx5);
	deleteMatrix(N_y + 3, pe_jumpy1); deleteMatrix(N_y + 3, pe_jumpy2); deleteMatrix(N_y + 3, pe_jumpy3); deleteMatrix(N_y + 3, pe_jumpy4); deleteMatrix(N_y + 3, pe_jumpy5);

	deleteMatrix(N_y, pi_jumpx1); deleteMatrix(N_y, pi_jumpx2); deleteMatrix(N_y, pi_jumpx3); deleteMatrix(N_y, pi_jumpx4); deleteMatrix(N_y, pi_jumpx5);
	deleteMatrix(N_y + 3, pi_jumpy1); deleteMatrix(N_y + 3, pi_jumpy2); deleteMatrix(N_y + 3, pi_jumpy3); deleteMatrix(N_y + 3, pi_jumpy4); deleteMatrix(N_y + 3, pi_jumpy5);

	deleteMatrix(N_y, pr_jumpx1); deleteMatrix(N_y, pr_jumpx2); deleteMatrix(N_y, pr_jumpx3); deleteMatrix(N_y, pr_jumpx4); deleteMatrix(N_y, pr_jumpx5);
	deleteMatrix(N_y + 3, pr_jumpy1); deleteMatrix(N_y + 3, pr_jumpy2); deleteMatrix(N_y + 3, pr_jumpy3); deleteMatrix(N_y + 3, pr_jumpy4); deleteMatrix(N_y + 3, pr_jumpy5);

	deleteMatrix(N_y, GQx_L1); deleteMatrix(N_y, GQx_Leta); deleteMatrix(N_y, GQx_Leta2);
	deleteMatrix(N_y, GQx_R1); deleteMatrix(N_y, GQx_Reta); deleteMatrix(N_y, GQx_Reta2);

	deleteMatrix(N_y, GQy_B1); deleteMatrix(N_y, GQy_Bxi); deleteMatrix(N_y, GQy_Bxi2);
	deleteMatrix(N_y, GQy_T1); deleteMatrix(N_y, GQy_Txi); deleteMatrix(N_y, GQy_Txi2);
}


// Diffusion
void compute_Talphaxi(double** u, double xi, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhoxi = rho_1[i][j] + rho_3[i][j] * eta + 3.0 * rho_4[i][j] * xi;
			double temp_Ealphaxi = Ealpha_1[i][j] + Ealpha_3[i][j] * eta + 3.0 * Ealpha_4[i][j] * xi;
			double temp_m1xi = m1_1[i][j] + m1_3[i][j] * eta + 3.0 * m1_4[i][j] * xi;
			double temp_m2xi = m2_1[i][j] + m2_3[i][j] * eta + 3.0 * m2_4[i][j] * xi;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_1 = -1.0 * temp_rhoxi / (C_valpha * temp_rho * temp_rho) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0);
			double temp_2 = 1.0 / (C_valpha * temp_rho) * (temp_Ealphaxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
			u[i][j] = temp_1 + temp_2;
		}
	}
}

void compute_Talphaeta(double** u, double xi, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Ealphaxi = Ealpha_2[i][j] + Ealpha_3[i][j] * xi + 3.0 * Ealpha_5[i][j] * eta;
			double temp_rhoxi = rho_2[i][j] + rho_3[i][j] * xi + 3.0 * rho_5[i][j] * eta;
			double temp_m1xi = m1_2[i][j] + m1_3[i][j] * xi + 3.0 * m1_5[i][j] * eta;
			double temp_m2xi = m2_2[i][j] + m2_3[i][j] * xi + 3.0 * m2_5[i][j] * eta;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_1 = -1.0 * temp_rhoxi / (C_valpha * temp_rho * temp_rho) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0);
			double temp_2 = 1.0 / (C_valpha * temp_rho) * (temp_Ealphaxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
			u[i][j] = temp_1 + temp_2;
		}
	}
}

void compute_Tr4xi(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_Erxi = Er_1[i][j] + Er_3[i][j] * eta + 3.0 * Er_4[i][j] * xi;
			double temp_rhoxi = rho_1[i][j] + rho_3[i][j] * eta + 3.0 * rho_4[i][j] * xi;
			double temp_m1xi = m1_1[i][j] + m1_3[i][j] * eta + 3.0 * m1_4[i][j] * xi;
			double temp_m2xi = m2_1[i][j] + m2_3[i][j] * eta + 3.0 * m2_4[i][j] * xi;

			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			u[i][j] = 1.0 / a * (temp_Erxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
		}
	}
}

void compute_Tr4eta(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_Erxi = Er_2[i][j] + Er_3[i][j] * xi + 3.0 * Er_5[i][j] * eta;
			double temp_rhoxi = rho_2[i][j] + rho_3[i][j] * xi + 3.0 * rho_5[i][j] * eta;
			double temp_m1xi = m1_2[i][j] + m1_3[i][j] * xi + 3.0 * m1_5[i][j] * eta;
			double temp_m2xi = m2_2[i][j] + m2_3[i][j] * xi + 3.0 * m2_5[i][j] * eta;

			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			u[i][j] = 1.0 / a * (temp_Erxi + temp_m1 * temp_m1 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + temp_m2 * temp_m2 * temp_rhoxi / (6.0 * temp_rho * temp_rho) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
		}
	}
}

void compute_Talphaxi2(double** u, double xi, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_rhoxi = rho_1[i][j] + rho_3[i][j] * eta + 3.0 * rho_4[i][j] * xi;
			double temp_Ealphaxi = Ealpha_1[i][j] + Ealpha_3[i][j] * eta + 3.0 * Ealpha_4[i][j] * xi;
			double temp_m1xi = m1_1[i][j] + m1_3[i][j] * eta + 3.0 * m1_4[i][j] * xi;
			double temp_m2xi = m2_1[i][j] + m2_3[i][j] * eta + 3.0 * m2_4[i][j] * xi;

			double temp_rhoxi2 = 3.0 * rho_4[i][j];
			double temp_Ealphaxi2 = 3.0 * Ealpha_4[i][j];
			double temp_m1xi2 = 3.0 * m1_4[i][j];
			double temp_m2xi2 = 3.0 * m2_4[i][j];

			double temp_1 = 1.0 / C_valpha * (2.0 * pow(temp_rhoxi, 2) / pow(temp_rho, 3) - temp_rhoxi2 / pow(temp_rho, 2)) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0);
			double temp_2 = -2.0 * temp_rhoxi / (C_valpha * pow(temp_rho, 2)) * (temp_Ealphaxi + pow(temp_m1, 2) * temp_rhoxi / (6.0 * pow(temp_rho, 2)) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + pow(temp_m2, 2) * temp_rhoxi / (6.0 * pow(temp_rho, 2)) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
			double temp_31 = -1.0 * pow(temp_m1, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m1, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m1 * temp_m1xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m1xi, 2) + temp_m1 * temp_m1xi2) / (3.0 * temp_rho);
			double temp_32 = -1.0 * pow(temp_m2, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m2, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m2 * temp_m2xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m2xi, 2) + temp_m2 * temp_m2xi2) / (3.0 * temp_rho);
			double temp_3 = 1.0 / (C_valpha * temp_rho) * (temp_Ealphaxi2 + temp_31 + temp_32);

			u[i][j] = temp_1 + temp_2 + temp_3;
		}
	}
}

void compute_Talphaeta2(double** u, double xi, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_Ealphaxi = Ealpha_2[i][j] + Ealpha_3[i][j] * xi + 3.0 * Ealpha_5[i][j] * eta;
			double temp_rhoxi = rho_2[i][j] + rho_3[i][j] * xi + 3.0 * rho_5[i][j] * eta;
			double temp_m1xi = m1_2[i][j] + m1_3[i][j] * xi + 3.0 * m1_5[i][j] * eta;
			double temp_m2xi = m2_2[i][j] + m2_3[i][j] * xi + 3.0 * m2_5[i][j] * eta;

			double temp_rhoxi2 = 3.0 * rho_5[i][j];
			double temp_Ealphaxi2 = 3.0 * Ealpha_5[i][j];
			double temp_m1xi2 = 3.0 * m1_5[i][j];
			double temp_m2xi2 = 3.0 * m2_5[i][j];

			double temp_1 = 1.0 / C_valpha * (2.0 * pow(temp_rhoxi, 2) / pow(temp_rho, 3) - temp_rhoxi2 / pow(temp_rho, 2)) * (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0);
			double temp_2 = -2.0 * temp_rhoxi / (C_valpha * pow(temp_rho, 2)) * (temp_Ealphaxi + pow(temp_m1, 2) * temp_rhoxi / (6.0 * pow(temp_rho, 2)) - temp_m1 * temp_m1xi / (3.0 * temp_rho) + pow(temp_m2, 2) * temp_rhoxi / (6.0 * pow(temp_rho, 2)) - temp_m2 * temp_m2xi / (3.0 * temp_rho));
			double temp_31 = -1.0 * pow(temp_m1, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m1, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m1 * temp_m1xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m1xi, 2) + temp_m1 * temp_m1xi2) / (3.0 * temp_rho);
			double temp_32 = -1.0 * pow(temp_m2, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m2, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m2 * temp_m2xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m2xi, 2) + temp_m2 * temp_m2xi2) / (3.0 * temp_rho);
			double temp_3 = 1.0 / (C_valpha * temp_rho) * (temp_Ealphaxi2 + temp_31 + temp_32);

			u[i][j] = temp_1 + temp_2 + temp_3;
		}
	}
}

void compute_Tr4xi2(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_rhoxi = rho_1[i][j] + rho_3[i][j] * eta + 3.0 * rho_4[i][j] * xi;
			double temp_m1xi = m1_1[i][j] + m1_3[i][j] * eta + 3.0 * m1_4[i][j] * xi;
			double temp_m2xi = m2_1[i][j] + m2_3[i][j] * eta + 3.0 * m2_4[i][j] * xi;

			double temp_Erxi2 = 3.0 * Er_4[i][j];
			double temp_rhoxi2 = 3.0 * rho_4[i][j];
			double temp_m1xi2 = 3.0 * m1_4[i][j];
			double temp_m2xi2 = 3.0 * m2_4[i][j];

			double temp_1 = -1.0 * pow(temp_m1, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m1, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m1 * temp_m1xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m1xi, 2) + temp_m1 * temp_m1xi2) / (3.0 * temp_rho);
			double temp_2 = -1.0 * pow(temp_m2, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m2, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m2 * temp_m2xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m2xi, 2) + temp_m2 * temp_m2xi2) / (3.0 * temp_rho);

			u[i][j] = 1.0 / a * (temp_Erxi2 + temp_1 + temp_2);
		}
	}
}

void compute_Tr4eta2(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m1 = m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_m2 = m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;

			double temp_rhoxi = rho_2[i][j] + rho_3[i][j] * xi + 3.0 * rho_5[i][j] * eta;
			double temp_m1xi = m1_2[i][j] + m1_3[i][j] * xi + 3.0 * m1_5[i][j] * eta;
			double temp_m2xi = m2_2[i][j] + m2_3[i][j] * xi + 3.0 * m2_5[i][j] * eta;

			double temp_Erxi2 = 3.0 * Er_5[i][j];
			double temp_rhoxi2 = 3.0 * rho_5[i][j];
			double temp_m1xi2 = 3.0 * m1_5[i][j];
			double temp_m2xi2 = 3.0 * m2_5[i][j];

			double temp_1 = -1.0 * pow(temp_m1, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m1, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m1 * temp_m1xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m1xi, 2) + temp_m1 * temp_m1xi2) / (3.0 * temp_rho);
			double temp_2 = -1.0 * pow(temp_m2, 2) * pow(temp_rhoxi, 2) / (3.0 * pow(temp_rho, 3)) + pow(temp_m2, 2) * temp_rhoxi2 / (6.0 * pow(temp_rho, 2)) + 2.0 * temp_m2 * temp_m2xi * temp_rhoxi / (3.0 * pow(temp_rho, 2)) - (pow(temp_m2xi, 2) + temp_m2 * temp_m2xi2) / (3.0 * temp_rho);

			u[i][j] = 1.0 / a * (temp_Erxi2 + temp_1 + temp_2);
		}
	}
}

void Talpha_jumpx(double** u, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talpha_plus; createMatrix(N_y + 2, N_x + 2, Talpha_plus);
	double** Talpha_minus; createMatrix(N_y + 2, N_x + 2, Talpha_minus);

	compute_Talpha(Talpha_plus, -1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talpha(Talpha_minus, 1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_x(u, Talpha_plus, Talpha_minus);

	deleteMatrix(N_y + 2, Talpha_plus); deleteMatrix(N_y + 2, Talpha_minus);
}

void Talpha_jumpy(double** u, double xi, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talpha_plus; createMatrix(N_y + 2, N_x + 2, Talpha_plus);
	double** Talpha_minus; createMatrix(N_y + 2, N_x + 2, Talpha_minus);

	compute_Talpha(Talpha_plus, xi, -1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talpha(Talpha_minus, xi, 1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_y(u, Talpha_plus, Talpha_minus);

	deleteMatrix(N_y + 2, Talpha_plus); deleteMatrix(N_y + 2, Talpha_minus);
}

void Tr4_jumpx(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talpha_plus; createMatrix(N_y + 2, N_x + 2, Talpha_plus);
	double** Talpha_minus; createMatrix(N_y + 2, N_x + 2, Talpha_minus);

	compute_Tr4(Talpha_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4(Talpha_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_jump_x(u, Talpha_plus, Talpha_minus);

	deleteMatrix(N_y + 2, Talpha_plus); deleteMatrix(N_y + 2, Talpha_minus);
}

void Tr4_jumpy(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talpha_plus; createMatrix(N_y + 2, N_x + 2, Talpha_plus);
	double** Talpha_minus; createMatrix(N_y + 2, N_x + 2, Talpha_minus);

	compute_Tr4(Talpha_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4(Talpha_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_jump_y(u, Talpha_plus, Talpha_minus);

	deleteMatrix(N_y + 2, Talpha_plus); deleteMatrix(N_y + 2, Talpha_minus);
}

void Talphaxi_ave(double** u, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talphaxi_plus; createMatrix(N_y + 2, N_x + 2, Talphaxi_plus);
	double** Talphaxi_minus; createMatrix(N_y + 2, N_x + 2, Talphaxi_minus);

	compute_Talphaxi(Talphaxi_plus, -1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_minus, 1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_ave_x(u, Talphaxi_plus, Talphaxi_minus);

	deleteMatrix(N_y + 2, Talphaxi_plus); deleteMatrix(N_y + 2, Talphaxi_minus);
}

void Talphaeta_ave(double** u, double xi, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talphaeta_plus; createMatrix(N_y + 2, N_x + 2, Talphaeta_plus);
	double** Talphaeta_minus; createMatrix(N_y + 2, N_x + 2, Talphaeta_minus);

	compute_Talphaeta(Talphaeta_plus, xi, -1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_minus, xi, 1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_ave_y(u, Talphaeta_plus, Talphaeta_minus);

	deleteMatrix(N_y + 2, Talphaeta_plus); deleteMatrix(N_y + 2, Talphaeta_minus);
}

void Tr4xi_ave(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talphaxi_plus; createMatrix(N_y + 2, N_x + 2, Talphaxi_plus);
	double** Talphaxi_minus; createMatrix(N_y + 2, N_x + 2, Talphaxi_minus);

	compute_Tr4xi(Talphaxi_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4xi(Talphaxi_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_ave_x(u, Talphaxi_plus, Talphaxi_minus);

	deleteMatrix(N_y + 2, Talphaxi_plus); deleteMatrix(N_y + 2, Talphaxi_minus);
}

void Tr4eta_ave(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talphaeta_plus; createMatrix(N_y + 2, N_x + 2, Talphaeta_plus);
	double** Talphaeta_minus; createMatrix(N_y + 2, N_x + 2, Talphaeta_minus);

	compute_Tr4eta(Talphaeta_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4eta(Talphaeta_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_ave_y(u, Talphaeta_plus, Talphaeta_minus);

	deleteMatrix(N_y + 2, Talphaeta_plus); deleteMatrix(N_y + 2, Talphaeta_minus);
}

void Talphaxi2_jump(double** u, double eta, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talphaxi2_plus; createMatrix(N_y + 2, N_x + 2, Talphaxi2_plus);
	double** Talphaxi2_minus; createMatrix(N_y + 2, N_x + 2, Talphaxi2_minus);

	compute_Talphaxi2(Talphaxi2_plus, -1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi2(Talphaxi2_minus, 1.0, eta, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_x(u, Talphaxi2_plus, Talphaxi2_minus);

	deleteMatrix(N_y + 2, Talphaxi2_plus); deleteMatrix(N_y + 2, Talphaxi2_minus);
}

void Talphaeta2_jump(double** u, double xi, double C_valpha, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** Talphaeta2_plus; createMatrix(N_y + 2, N_x + 2, Talphaeta2_plus);
	double** Talphaeta2_minus; createMatrix(N_y + 2, N_x + 2, Talphaeta2_minus);

	compute_Talphaeta2(Talphaeta2_plus, xi, -1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta2(Talphaeta2_minus, xi, 1.0, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_jump_y(u, Talphaeta2_plus, Talphaeta2_minus);

	deleteMatrix(N_y + 2, Talphaeta2_plus); deleteMatrix(N_y + 2, Talphaeta2_minus);
}

void Tr4xi2_jump(double** u, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talphaxi2_plus; createMatrix(N_y + 2, N_x + 2, Talphaxi2_plus);
	double** Talphaxi2_minus; createMatrix(N_y + 2, N_x + 2, Talphaxi2_minus);

	compute_Tr4xi2(Talphaxi2_plus, -1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4xi2(Talphaxi2_minus, 1.0, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_jump_x(u, Talphaxi2_plus, Talphaxi2_minus);

	deleteMatrix(N_y + 2, Talphaxi2_plus); deleteMatrix(N_y + 2, Talphaxi2_minus);
}

void Tr4eta2_jump(double** u, double xi, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** Talphaeta2_plus; createMatrix(N_y + 2, N_x + 2, Talphaeta2_plus);
	double** Talphaeta2_minus; createMatrix(N_y + 2, N_x + 2, Talphaeta2_minus);

	compute_Tr4eta2(Talphaeta2_plus, xi, -1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
	compute_Tr4eta2(Talphaeta2_minus, xi, 1.0, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	compute_jump_y(u, Talphaeta2_plus, Talphaeta2_minus);

	deleteMatrix(N_y + 2, Talphaeta2_plus); deleteMatrix(N_y + 2, Talphaeta2_minus);
}

void Diffusion_alpha(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double C_valpha, double kappa_alpha, double alpha_d, double beta_d, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double h_x, double h_y)
{
	double** Talpha_jx1; createMatrix(N_y, N_x + 3, Talpha_jx1); double** Talpha_jx2; createMatrix(N_y, N_x + 3, Talpha_jx2);
	double** Talpha_jx3; createMatrix(N_y, N_x + 3, Talpha_jx3); double** Talpha_jx4; createMatrix(N_y, N_x + 3, Talpha_jx4);
	double** Talpha_jx5; createMatrix(N_y, N_x + 3, Talpha_jx5);

	double** Talpha_jy1; createMatrix(N_y + 3, N_x, Talpha_jy1); double** Talpha_jy2; createMatrix(N_y + 3, N_x, Talpha_jy2);
	double** Talpha_jy3; createMatrix(N_y + 3, N_x, Talpha_jy3); double** Talpha_jy4; createMatrix(N_y + 3, N_x, Talpha_jy4);
	double** Talpha_jy5; createMatrix(N_y + 3, N_x, Talpha_jy5);

	double** Talphaxi_a1; createMatrix(N_y, N_x + 3, Talphaxi_a1); double** Talphaxi_a2; createMatrix(N_y, N_x + 3, Talphaxi_a2);
	double** Talphaxi_a3; createMatrix(N_y, N_x + 3, Talphaxi_a3); double** Talphaxi_a4; createMatrix(N_y, N_x + 3, Talphaxi_a4);
	double** Talphaxi_a5; createMatrix(N_y, N_x + 3, Talphaxi_a5);

	double** Talphaeta_a1; createMatrix(N_y + 3, N_x, Talphaeta_a1); double** Talphaeta_a2; createMatrix(N_y + 3, N_x, Talphaeta_a2);
	double** Talphaeta_a3; createMatrix(N_y + 3, N_x, Talphaeta_a3); double** Talphaeta_a4; createMatrix(N_y + 3, N_x, Talphaeta_a4);
	double** Talphaeta_a5; createMatrix(N_y + 3, N_x, Talphaeta_a5);

	double** Talphaxi2_j1; createMatrix(N_y, N_x + 3, Talphaxi2_j1); double** Talphaxi2_j2; createMatrix(N_y, N_x + 3, Talphaxi2_j2);
	double** Talphaxi2_j3; createMatrix(N_y, N_x + 3, Talphaxi2_j3); double** Talphaxi2_j4; createMatrix(N_y, N_x + 3, Talphaxi2_j4);
	double** Talphaxi2_j5; createMatrix(N_y, N_x + 3, Talphaxi2_j5);

	double** Talphaeta2_j1; createMatrix(N_y + 3, N_x, Talphaeta2_j1); double** Talphaeta2_j2; createMatrix(N_y + 3, N_x, Talphaeta2_j2);
	double** Talphaeta2_j3; createMatrix(N_y + 3, N_x, Talphaeta2_j3); double** Talphaeta2_j4; createMatrix(N_y + 3, N_x, Talphaeta2_j4);
	double** Talphaeta2_j5; createMatrix(N_y + 3, N_x, Talphaeta2_j5);

	double** Talphaxi_11; createMatrix(N_y + 2, N_x + 2, Talphaxi_11); double** Talphaxi_12; createMatrix(N_y + 2, N_x + 2, Talphaxi_12);
	double** Talphaxi_13; createMatrix(N_y + 2, N_x + 2, Talphaxi_13); double** Talphaxi_14; createMatrix(N_y + 2, N_x + 2, Talphaxi_14);
	double** Talphaxi_15; createMatrix(N_y + 2, N_x + 2, Talphaxi_15);

	double** Talphaxi_21; createMatrix(N_y + 2, N_x + 2, Talphaxi_21); double** Talphaxi_22; createMatrix(N_y + 2, N_x + 2, Talphaxi_22);
	double** Talphaxi_23; createMatrix(N_y + 2, N_x + 2, Talphaxi_23); double** Talphaxi_24; createMatrix(N_y + 2, N_x + 2, Talphaxi_24);
	double** Talphaxi_25; createMatrix(N_y + 2, N_x + 2, Talphaxi_25);

	double** Talphaxi_31; createMatrix(N_y + 2, N_x + 2, Talphaxi_31); double** Talphaxi_32; createMatrix(N_y + 2, N_x + 2, Talphaxi_32);
	double** Talphaxi_33; createMatrix(N_y + 2, N_x + 2, Talphaxi_33); double** Talphaxi_34; createMatrix(N_y + 2, N_x + 2, Talphaxi_34);
	double** Talphaxi_35; createMatrix(N_y + 2, N_x + 2, Talphaxi_35);

	double** Talphaxi_41; createMatrix(N_y + 2, N_x + 2, Talphaxi_41); double** Talphaxi_42; createMatrix(N_y + 2, N_x + 2, Talphaxi_42);
	double** Talphaxi_43; createMatrix(N_y + 2, N_x + 2, Talphaxi_43); double** Talphaxi_44; createMatrix(N_y + 2, N_x + 2, Talphaxi_44);
	double** Talphaxi_45; createMatrix(N_y + 2, N_x + 2, Talphaxi_45);

	double** Talphaxi_51; createMatrix(N_y + 2, N_x + 2, Talphaxi_51); double** Talphaxi_52; createMatrix(N_y + 2, N_x + 2, Talphaxi_52);
	double** Talphaxi_53; createMatrix(N_y + 2, N_x + 2, Talphaxi_53); double** Talphaxi_54; createMatrix(N_y + 2, N_x + 2, Talphaxi_54);
	double** Talphaxi_55; createMatrix(N_y + 2, N_x + 2, Talphaxi_55);

	double** Talphaeta_11; createMatrix(N_y + 2, N_x + 2, Talphaeta_11); double** Talphaeta_12; createMatrix(N_y + 2, N_x + 2, Talphaeta_12);
	double** Talphaeta_13; createMatrix(N_y + 2, N_x + 2, Talphaeta_13); double** Talphaeta_14; createMatrix(N_y + 2, N_x + 2, Talphaeta_14);
	double** Talphaeta_15; createMatrix(N_y + 2, N_x + 2, Talphaeta_15);

	double** Talphaeta_21; createMatrix(N_y + 2, N_x + 2, Talphaeta_21); double** Talphaeta_22; createMatrix(N_y + 2, N_x + 2, Talphaeta_22);
	double** Talphaeta_23; createMatrix(N_y + 2, N_x + 2, Talphaeta_23); double** Talphaeta_24; createMatrix(N_y + 2, N_x + 2, Talphaeta_24);
	double** Talphaeta_25; createMatrix(N_y + 2, N_x + 2, Talphaeta_25);

	double** Talphaeta_31; createMatrix(N_y + 2, N_x + 2, Talphaeta_31); double** Talphaeta_32; createMatrix(N_y + 2, N_x + 2, Talphaeta_32);
	double** Talphaeta_33; createMatrix(N_y + 2, N_x + 2, Talphaeta_33); double** Talphaeta_34; createMatrix(N_y + 2, N_x + 2, Talphaeta_34);
	double** Talphaeta_35; createMatrix(N_y + 2, N_x + 2, Talphaeta_35);

	double** Talphaeta_41; createMatrix(N_y + 2, N_x + 2, Talphaeta_41); double** Talphaeta_42; createMatrix(N_y + 2, N_x + 2, Talphaeta_42);
	double** Talphaeta_43; createMatrix(N_y + 2, N_x + 2, Talphaeta_43); double** Talphaeta_44; createMatrix(N_y + 2, N_x + 2, Talphaeta_44);
	double** Talphaeta_45; createMatrix(N_y + 2, N_x + 2, Talphaeta_45);

	double** Talphaeta_51; createMatrix(N_y + 2, N_x + 2, Talphaeta_51); double** Talphaeta_52; createMatrix(N_y + 2, N_x + 2, Talphaeta_52);
	double** Talphaeta_53; createMatrix(N_y + 2, N_x + 2, Talphaeta_53); double** Talphaeta_54; createMatrix(N_y + 2, N_x + 2, Talphaeta_54);
	double** Talphaeta_55; createMatrix(N_y + 2, N_x + 2, Talphaeta_55);


	// needed modified
	Talpha_jumpx(Talpha_jx1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpx(Talpha_jx2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpx(Talpha_jx3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpx(Talpha_jx4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpx(Talpha_jx5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Talpha_jumpy(Talpha_jy1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpy(Talpha_jy2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpy(Talpha_jy3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpy(Talpha_jy4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talpha_jumpy(Talpha_jy5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Talphaxi_ave(Talphaxi_a1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi_ave(Talphaxi_a2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi_ave(Talphaxi_a3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi_ave(Talphaxi_a4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi_ave(Talphaxi_a5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Talphaeta_ave(Talphaeta_a1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta_ave(Talphaeta_a2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta_ave(Talphaeta_a3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta_ave(Talphaeta_a4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta_ave(Talphaeta_a5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Talphaxi2_jump(Talphaxi2_j1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi2_jump(Talphaxi2_j2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi2_jump(Talphaxi2_j3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi2_jump(Talphaxi2_j4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaxi2_jump(Talphaxi2_j5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Talphaeta2_jump(Talphaeta2_j1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta2_jump(Talphaeta2_j2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta2_jump(Talphaeta2_j3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta2_jump(Talphaeta2_j4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Talphaeta2_jump(Talphaeta2_j5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaxi(Talphaxi_11, GQxi_1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_12, GQxi_1, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_13, GQxi_1, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_14, GQxi_1, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_15, GQxi_1, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaxi(Talphaxi_21, GQxi_2, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_22, GQxi_2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_23, GQxi_2, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_24, GQxi_2, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_25, GQxi_2, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaxi(Talphaxi_31, GQxi_3, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_32, GQxi_3, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_33, GQxi_3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_34, GQxi_3, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_35, GQxi_3, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaxi(Talphaxi_41, GQxi_4, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_42, GQxi_4, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_43, GQxi_4, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_44, GQxi_4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_45, GQxi_4, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaxi(Talphaxi_51, GQxi_5, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_52, GQxi_5, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_53, GQxi_5, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_54, GQxi_5, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaxi(Talphaxi_55, GQxi_5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaeta(Talphaeta_11, GQxi_1, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_12, GQxi_1, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_13, GQxi_1, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_14, GQxi_1, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_15, GQxi_1, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaeta(Talphaeta_21, GQxi_2, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_22, GQxi_2, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_23, GQxi_2, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_24, GQxi_2, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_25, GQxi_2, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaeta(Talphaeta_31, GQxi_3, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_32, GQxi_3, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_33, GQxi_3, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_34, GQxi_3, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_35, GQxi_3, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaeta(Talphaeta_41, GQxi_4, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_42, GQxi_4, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_43, GQxi_4, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_44, GQxi_4, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_45, GQxi_4, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Talphaeta(Talphaeta_51, GQxi_5, GQxi_1, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_52, GQxi_5, GQxi_2, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_53, GQxi_5, GQxi_3, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_54, GQxi_5, GQxi_4, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Talphaeta(Talphaeta_55, GQxi_5, GQxi_5, C_valpha, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	double** GQ_xi_1; createMatrix(N_y + 2, N_x + 2, GQ_xi_1); double** GQ_xi_eta; createMatrix(N_y + 2, N_x + 2, GQ_xi_eta); double** GQ_xi_xi; createMatrix(N_y + 2, N_x + 2, GQ_xi_xi);
	double** GQ_eta_1; createMatrix(N_y + 2, N_x + 2, GQ_eta_1); double** GQ_eta_xi; createMatrix(N_y + 2, N_x + 2, GQ_eta_xi); double** GQ_eta_eta; createMatrix(N_y + 2, N_x + 2, GQ_eta_eta);

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
			GQ_xi_1[i][j] = -1.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] + w_12 * Talphaxi_12[i][j] + w_13 * Talphaxi_13[i][j] + w_14 * Talphaxi_14[i][j] + w_15 * Talphaxi_15[i][j] + w_21 * Talphaxi_21[i][j] + w_22 * Talphaxi_22[i][j] + w_23 * Talphaxi_23[i][j] + w_24 * Talphaxi_24[i][j] + w_25 * Talphaxi_25[i][j] + w_31 * Talphaxi_31[i][j] + w_32 * Talphaxi_32[i][j] + w_33 * Talphaxi_33[i][j] + w_34 * Talphaxi_34[i][j] + w_35 * Talphaxi_35[i][j] + w_41 * Talphaxi_41[i][j] + w_42 * Talphaxi_42[i][j] + w_43 * Talphaxi_43[i][j] + w_44 * Talphaxi_44[i][j] + w_45 * Talphaxi_45[i][j] + w_51 * Talphaxi_51[i][j] + w_52 * Talphaxi_52[i][j] + w_53 * Talphaxi_53[i][j] + w_54 * Talphaxi_54[i][j] + w_55 * Talphaxi_55[i][j]);
			GQ_xi_eta[i][j] = -1.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] * GQxi_1 + w_12 * Talphaxi_12[i][j] * GQxi_2 + w_13 * Talphaxi_13[i][j] * GQxi_3 + w_14 * Talphaxi_14[i][j] * GQxi_4 + w_15 * Talphaxi_15[i][j] * GQxi_5 + w_21 * Talphaxi_21[i][j] * GQxi_1 + w_22 * Talphaxi_22[i][j] * GQxi_2 + w_23 * Talphaxi_23[i][j] * GQxi_3 + w_24 * Talphaxi_24[i][j] * GQxi_4 + w_25 * Talphaxi_25[i][j] * GQxi_5 + w_31 * Talphaxi_31[i][j] * GQxi_1 + w_32 * Talphaxi_32[i][j] * GQxi_2 + w_33 * Talphaxi_33[i][j] * GQxi_3 + w_34 * Talphaxi_34[i][j] * GQxi_4 + w_35 * Talphaxi_35[i][j] * GQxi_5 + w_41 * Talphaxi_41[i][j] * GQxi_1 + w_42 * Talphaxi_42[i][j] * GQxi_2 + w_43 * Talphaxi_43[i][j] * GQxi_3 + w_44 * Talphaxi_44[i][j] * GQxi_4 + w_45 * Talphaxi_45[i][j] * GQxi_5 + w_51 * Talphaxi_51[i][j] * GQxi_1 + w_52 * Talphaxi_52[i][j] * GQxi_2 + w_53 * Talphaxi_53[i][j] * GQxi_3 + w_54 * Talphaxi_54[i][j] * GQxi_4 + w_55 * Talphaxi_55[i][j] * GQxi_5);
			GQ_xi_xi[i][j] = -3.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] * GQxi_1 + w_12 * Talphaxi_12[i][j] * GQxi_1 + w_13 * Talphaxi_13[i][j] * GQxi_1 + w_14 * Talphaxi_14[i][j] * GQxi_1 + w_15 * Talphaxi_15[i][j] * GQxi_1 + w_21 * Talphaxi_21[i][j] * GQxi_2 + w_22 * Talphaxi_22[i][j] * GQxi_2 + w_23 * Talphaxi_23[i][j] * GQxi_2 + w_24 * Talphaxi_24[i][j] * GQxi_2 + w_25 * Talphaxi_25[i][j] * GQxi_2 + w_31 * Talphaxi_31[i][j] * GQxi_3 + w_32 * Talphaxi_32[i][j] * GQxi_3 + w_33 * Talphaxi_33[i][j] * GQxi_3 + w_34 * Talphaxi_34[i][j] * GQxi_3 + w_35 * Talphaxi_35[i][j] * GQxi_3 + w_41 * Talphaxi_41[i][j] * GQxi_4 + w_42 * Talphaxi_42[i][j] * GQxi_4 + w_43 * Talphaxi_43[i][j] * GQxi_4 + w_44 * Talphaxi_44[i][j] * GQxi_4 + w_45 * Talphaxi_45[i][j] * GQxi_4 + w_51 * Talphaxi_51[i][j] * GQxi_5 + w_52 * Talphaxi_52[i][j] * GQxi_5 + w_53 * Talphaxi_53[i][j] * GQxi_5 + w_54 * Talphaxi_54[i][j] * GQxi_5 + w_55 * Talphaxi_55[i][j] * GQxi_5);


			GQ_eta_1[i][j] = -1.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] + w_12 * Talphaeta_12[i][j] + w_13 * Talphaeta_13[i][j] + w_14 * Talphaeta_14[i][j] + w_15 * Talphaeta_15[i][j] + w_21 * Talphaeta_21[i][j] + w_22 * Talphaeta_22[i][j] + w_23 * Talphaeta_23[i][j] + w_24 * Talphaeta_24[i][j] + w_25 * Talphaeta_25[i][j] + w_31 * Talphaeta_31[i][j] + w_32 * Talphaeta_32[i][j] + w_33 * Talphaeta_33[i][j] + w_34 * Talphaeta_34[i][j] + w_35 * Talphaeta_35[i][j] + w_41 * Talphaeta_41[i][j] + w_42 * Talphaeta_42[i][j] + w_43 * Talphaeta_43[i][j] + w_44 * Talphaeta_44[i][j] + w_45 * Talphaeta_45[i][j] + w_51 * Talphaeta_51[i][j] + w_52 * Talphaeta_52[i][j] + w_53 * Talphaeta_53[i][j] + w_54 * Talphaeta_54[i][j] + w_55 * Talphaeta_55[i][j]);
			GQ_eta_xi[i][j] = -1.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] * GQxi_1 + w_12 * Talphaeta_12[i][j] * GQxi_1 + w_13 * Talphaeta_13[i][j] * GQxi_1 + w_14 * Talphaeta_14[i][j] * GQxi_1 + w_15 * Talphaeta_15[i][j] * GQxi_1 + w_21 * Talphaeta_21[i][j] * GQxi_2 + w_22 * Talphaeta_22[i][j] * GQxi_2 + w_23 * Talphaeta_23[i][j] * GQxi_2 + w_24 * Talphaeta_24[i][j] * GQxi_2 + w_25 * Talphaeta_25[i][j] * GQxi_2 + w_31 * Talphaeta_31[i][j] * GQxi_3 + w_32 * Talphaeta_32[i][j] * GQxi_3 + w_33 * Talphaeta_33[i][j] * GQxi_3 + w_34 * Talphaeta_34[i][j] * GQxi_3 + w_35 * Talphaeta_35[i][j] * GQxi_3 + w_41 * Talphaeta_41[i][j] * GQxi_4 + w_42 * Talphaeta_42[i][j] * GQxi_4 + w_43 * Talphaeta_43[i][j] * GQxi_4 + w_44 * Talphaeta_44[i][j] * GQxi_4 + w_45 * Talphaeta_45[i][j] * GQxi_4 + w_51 * Talphaeta_51[i][j] * GQxi_5 + w_52 * Talphaeta_52[i][j] * GQxi_5 + w_53 * Talphaeta_53[i][j] * GQxi_5 + w_54 * Talphaeta_54[i][j] * GQxi_5 + w_55 * Talphaeta_55[i][j] * GQxi_5);
			GQ_eta_eta[i][j] = -3.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] * GQxi_1 + w_12 * Talphaeta_12[i][j] * GQxi_2 + w_13 * Talphaeta_13[i][j] * GQxi_3 + w_14 * Talphaeta_14[i][j] * GQxi_4 + w_15 * Talphaeta_15[i][j] * GQxi_5 + w_21 * Talphaeta_21[i][j] * GQxi_1 + w_22 * Talphaeta_22[i][j] * GQxi_2 + w_23 * Talphaeta_23[i][j] * GQxi_3 + w_24 * Talphaeta_24[i][j] * GQxi_4 + w_25 * Talphaeta_25[i][j] * GQxi_5 + w_31 * Talphaeta_31[i][j] * GQxi_1 + w_32 * Talphaeta_32[i][j] * GQxi_2 + w_33 * Talphaeta_33[i][j] * GQxi_3 + w_34 * Talphaeta_34[i][j] * GQxi_4 + w_35 * Talphaeta_35[i][j] * GQxi_5 + w_41 * Talphaeta_41[i][j] * GQxi_1 + w_42 * Talphaeta_42[i][j] * GQxi_2 + w_43 * Talphaeta_43[i][j] * GQxi_3 + w_44 * Talphaeta_44[i][j] * GQxi_4 + w_45 * Talphaeta_45[i][j] * GQxi_5 + w_51 * Talphaeta_51[i][j] * GQxi_1 + w_52 * Talphaeta_52[i][j] * GQxi_2 + w_53 * Talphaeta_53[i][j] * GQxi_3 + w_54 * Talphaeta_54[i][j] * GQxi_4 + w_55 * Talphaeta_55[i][j] * GQxi_5);

		}
	}

	double** GQ_Talpha_jx; createMatrix(N_y, N_x + 3, GQ_Talpha_jx); double** GQ_Talphaxi_a; createMatrix(N_y, N_x + 3, GQ_Talphaxi_a); double** GQ_Talphaxi2_j; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_j);
	double** GQ_Talpha_jxeta; createMatrix(N_y, N_x + 3, GQ_Talpha_jxeta); double** GQ_Talphaxi_aeta; createMatrix(N_y, N_x + 3, GQ_Talphaxi_aeta); double** GQ_Talphaxi2_jeta; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_jeta);
	double** GQ_Talpha_jxeta2; createMatrix(N_y, N_x + 3, GQ_Talpha_jxeta2); double** GQ_Talphaxi_aeta2; createMatrix(N_y, N_x + 3, GQ_Talphaxi_aeta2); double** GQ_Talphaxi2_jeta2; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_jeta2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 3; j++)
		{
			GQ_Talpha_jx[i][j] = weight_1 * Talpha_jx1[i][j] + weight_2 * Talpha_jx2[i][j] + weight_3 * Talpha_jx3[i][j] + weight_4 * Talpha_jx4[i][j] + weight_5 * Talpha_jx5[i][j];
			GQ_Talphaxi_a[i][j] = weight_1 * Talphaxi_a1[i][j] + weight_2 * Talphaxi_a2[i][j] + weight_3 * Talphaxi_a3[i][j] + weight_4 * Talphaxi_a4[i][j] + weight_5 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_j[i][j] = weight_1 * Talphaxi2_j1[i][j] + weight_2 * Talphaxi2_j2[i][j] + weight_3 * Talphaxi2_j3[i][j] + weight_4 * Talphaxi2_j4[i][j] + weight_5 * Talphaxi2_j5[i][j];

			GQ_Talpha_jxeta[i][j] = weight_1 * GQxi_1 * Talpha_jx1[i][j] + weight_2 * GQxi_2 * Talpha_jx2[i][j] + weight_3 * GQxi_3 * Talpha_jx3[i][j] + weight_4 * GQxi_4 * Talpha_jx4[i][j] + weight_5 * GQxi_5 * Talpha_jx5[i][j];
			GQ_Talphaxi_aeta[i][j] = weight_1 * GQxi_1 * Talphaxi_a1[i][j] + weight_2 * GQxi_2 * Talphaxi_a2[i][j] + weight_3 * GQxi_3 * Talphaxi_a3[i][j] + weight_4 * GQxi_4 * Talphaxi_a4[i][j] + weight_5 * GQxi_5 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_jeta[i][j] = weight_1 * GQxi_1 * Talphaxi2_j1[i][j] + weight_2 * GQxi_2 * Talphaxi2_j2[i][j] + weight_3 * GQxi_3 * Talphaxi2_j3[i][j] + weight_4 * GQxi_4 * Talphaxi2_j4[i][j] + weight_5 * GQxi_5 * Talphaxi2_j5[i][j];

			GQ_Talpha_jxeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talpha_jx1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talpha_jx2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talpha_jx3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talpha_jx4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talpha_jx5[i][j];
			GQ_Talphaxi_aeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaxi_a1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaxi_a2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaxi_a3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaxi_a4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_jeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaxi2_j1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaxi2_j2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaxi2_j3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaxi2_j4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaxi2_j5[i][j];
		}
	}

	double** GQ_Talpha_jy; createMatrix(N_y + 3, N_x, GQ_Talpha_jy); double** GQ_Talphaeta_a; createMatrix(N_y + 3, N_x, GQ_Talphaeta_a); double** GQ_Talphaeta2_j; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_j);
	double** GQ_Talpha_jyxi; createMatrix(N_y + 3, N_x, GQ_Talpha_jyxi); double** GQ_Talphaeta_axi; createMatrix(N_y + 3, N_x, GQ_Talphaeta_axi); double** GQ_Talphaeta2_jxi; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_jxi);
	double** GQ_Talpha_jyxi2; createMatrix(N_y + 3, N_x, GQ_Talpha_jyxi2); double** GQ_Talphaeta_axi2; createMatrix(N_y + 3, N_x, GQ_Talphaeta_axi2); double** GQ_Talphaeta2_jxi2; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_jxi2);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 3; i++)
		{
			GQ_Talpha_jy[i][j] = weight_1 * Talpha_jy1[i][j] + weight_2 * Talpha_jy2[i][j] + weight_3 * Talpha_jy3[i][j] + weight_4 * Talpha_jy4[i][j] + weight_5 * Talpha_jy5[i][j];
			GQ_Talphaeta_a[i][j] = weight_1 * Talphaeta_a1[i][j] + weight_2 * Talphaeta_a2[i][j] + weight_3 * Talphaeta_a3[i][j] + weight_4 * Talphaeta_a4[i][j] + weight_5 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_j[i][j] = weight_1 * Talphaeta2_j1[i][j] + weight_2 * Talphaeta2_j2[i][j] + weight_3 * Talphaeta2_j3[i][j] + weight_4 * Talphaeta2_j4[i][j] + weight_5 * Talphaeta2_j5[i][j];

			GQ_Talpha_jyxi[i][j] = weight_1 * GQxi_1 * Talpha_jy1[i][j] + weight_2 * GQxi_2 * Talpha_jy2[i][j] + weight_3 * GQxi_3 * Talpha_jy3[i][j] + weight_4 * GQxi_4 * Talpha_jy4[i][j] + weight_5 * GQxi_5 * Talpha_jy5[i][j];
			GQ_Talphaeta_axi[i][j] = weight_1 * GQxi_1 * Talphaeta_a1[i][j] + weight_2 * GQxi_2 * Talphaeta_a2[i][j] + weight_3 * GQxi_3 * Talphaeta_a3[i][j] + weight_4 * GQxi_4 * Talphaeta_a4[i][j] + weight_5 * GQxi_5 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_jxi[i][j] = weight_1 * GQxi_1 * Talphaeta2_j1[i][j] + weight_2 * GQxi_2 * Talphaeta2_j2[i][j] + weight_3 * GQxi_3 * Talphaeta2_j3[i][j] + weight_4 * GQxi_4 * Talphaeta2_j4[i][j] + weight_5 * GQxi_5 * Talphaeta2_j5[i][j];

			GQ_Talpha_jyxi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talpha_jy1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talpha_jy2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talpha_jy3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talpha_jy4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talpha_jy5[i][j];
			GQ_Talphaeta_axi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaeta_a1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaeta_a2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaeta_a3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaeta_a4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_jxi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaeta2_j1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaeta2_j2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaeta2_j3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaeta2_j4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaeta2_j5[i][j];
		}
	}


	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i][j] = 0.0 + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] - GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] - GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] - GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] - GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] - GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] - GQ_Talphaeta2_j[i][j - 1]);
			U_1[i][j] = GQ_xi_1[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] + GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi[i + 1][j - 1] - GQ_Talphaeta_axi[i][j - 1]) + (alpha_d - kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] + GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi[i + 1][j - 1] - GQ_Talpha_jyxi[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] + GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi[i + 1][j - 1] - GQ_Talphaeta2_jxi[i][j - 1]);
			U_2[i][j] = GQ_eta_1[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta[i - 1][j + 1] - GQ_Talphaxi_aeta[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] + GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta[i - 1][j + 1] - GQ_Talpha_jxeta[i - 1][j]) + (alpha_d - kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] + GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta[i - 1][j + 1] - GQ_Talphaxi2_jeta[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] + GQ_Talphaeta2_j[i][j - 1]);
			U_3[i][j] = GQ_xi_eta[i][j] + GQ_eta_xi[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta[i - 1][j + 1] + GQ_Talphaxi_aeta[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi[i + 1][j - 1] + GQ_Talphaeta_axi[i][j - 1]) + (alpha_d - kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta[i - 1][j + 1] + GQ_Talpha_jxeta[i - 1][j]) + (alpha_d - kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi[i + 1][j - 1] + GQ_Talpha_jyxi[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta[i - 1][j + 1] + GQ_Talphaxi2_jeta[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi[i + 1][j - 1] + GQ_Talphaeta2_jxi[i][j - 1]);
			U_4[i][j] = GQ_xi_xi[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] - GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi2[i + 1][j - 1] - GQ_Talphaeta_axi2[i][j - 1]) + (alpha_d - 3.0 * kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] - GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi2[i + 1][j - 1] - GQ_Talpha_jyxi2[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] - GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi2[i + 1][j - 1] - GQ_Talphaeta2_jxi2[i][j - 1]);
			U_5[i][j] = GQ_eta_eta[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta2[i - 1][j + 1] - GQ_Talphaxi_aeta2[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] - GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta2[i - 1][j + 1] - GQ_Talpha_jxeta2[i - 1][j]) + (alpha_d - 3.0 * kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] - GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta2[i - 1][j + 1] - GQ_Talphaxi2_jeta2[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] - GQ_Talphaeta2_j[i][j - 1]);
		}
	}

	deleteMatrix(N_y, Talpha_jx1); deleteMatrix(N_y, Talpha_jx2); deleteMatrix(N_y, Talpha_jx3); deleteMatrix(N_y, Talpha_jx4); deleteMatrix(N_y, Talpha_jx5);
	deleteMatrix(N_y + 3, Talpha_jy1); deleteMatrix(N_y + 3, Talpha_jy2); deleteMatrix(N_y + 3, Talpha_jy3); deleteMatrix(N_y + 3, Talpha_jy4); deleteMatrix(N_y + 3, Talpha_jy5);
	deleteMatrix(N_y, Talphaxi_a1); deleteMatrix(N_y, Talphaxi_a2); deleteMatrix(N_y, Talphaxi_a3); deleteMatrix(N_y, Talphaxi_a4); deleteMatrix(N_y, Talphaxi_a5);
	deleteMatrix(N_y + 3, Talphaeta_a1); deleteMatrix(N_y + 3, Talphaeta_a2); deleteMatrix(N_y + 3, Talphaeta_a3); deleteMatrix(N_y + 3, Talphaeta_a4); deleteMatrix(N_y + 3, Talphaeta_a5);
	deleteMatrix(N_y, Talphaxi2_j1); deleteMatrix(N_y, Talphaxi2_j2); deleteMatrix(N_y, Talphaxi2_j3); deleteMatrix(N_y, Talphaxi2_j4); deleteMatrix(N_y, Talphaxi2_j5);
	deleteMatrix(N_y + 3, Talphaeta2_j1); deleteMatrix(N_y + 3, Talphaeta2_j2); deleteMatrix(N_y + 3, Talphaeta2_j3); deleteMatrix(N_y + 3, Talphaeta2_j4); deleteMatrix(N_y + 3, Talphaeta2_j5);

	deleteMatrix(N_y, GQ_Talpha_jx); deleteMatrix(N_y, GQ_Talphaxi_a); deleteMatrix(N_y, GQ_Talphaxi2_j);
	deleteMatrix(N_y, GQ_Talpha_jxeta); deleteMatrix(N_y, GQ_Talphaxi_aeta); deleteMatrix(N_y, GQ_Talphaxi2_jeta);
	deleteMatrix(N_y, GQ_Talpha_jxeta2); deleteMatrix(N_y, GQ_Talphaxi_aeta2); deleteMatrix(N_y, GQ_Talphaxi2_jeta2);

	deleteMatrix(N_y + 3, GQ_Talpha_jy); deleteMatrix(N_y + 3, GQ_Talphaeta_a); deleteMatrix(N_y + 3, GQ_Talphaeta2_j);
	deleteMatrix(N_y + 3, GQ_Talpha_jyxi); deleteMatrix(N_y + 3, GQ_Talphaeta_axi); deleteMatrix(N_y + 3, GQ_Talphaeta2_jxi);
	deleteMatrix(N_y + 3, GQ_Talpha_jyxi2); deleteMatrix(N_y + 3, GQ_Talphaeta_axi2); deleteMatrix(N_y + 3, GQ_Talphaeta2_jxi2);

	deleteMatrix(N_y + 2, Talphaxi_11); deleteMatrix(N_y + 2, Talphaxi_12); deleteMatrix(N_y + 2, Talphaxi_13); deleteMatrix(N_y + 2, Talphaxi_14); deleteMatrix(N_y + 2, Talphaxi_15);
	deleteMatrix(N_y + 2, Talphaxi_21); deleteMatrix(N_y + 2, Talphaxi_22); deleteMatrix(N_y + 2, Talphaxi_23); deleteMatrix(N_y + 2, Talphaxi_24); deleteMatrix(N_y + 2, Talphaxi_25);
	deleteMatrix(N_y + 2, Talphaxi_31); deleteMatrix(N_y + 2, Talphaxi_32); deleteMatrix(N_y + 2, Talphaxi_33); deleteMatrix(N_y + 2, Talphaxi_34); deleteMatrix(N_y + 2, Talphaxi_35);
	deleteMatrix(N_y + 2, Talphaxi_41); deleteMatrix(N_y + 2, Talphaxi_42); deleteMatrix(N_y + 2, Talphaxi_43); deleteMatrix(N_y + 2, Talphaxi_44); deleteMatrix(N_y + 2, Talphaxi_45);
	deleteMatrix(N_y + 2, Talphaxi_51); deleteMatrix(N_y + 2, Talphaxi_52); deleteMatrix(N_y + 2, Talphaxi_53); deleteMatrix(N_y + 2, Talphaxi_54); deleteMatrix(N_y + 2, Talphaxi_55);

	deleteMatrix(N_y + 2, Talphaeta_11); deleteMatrix(N_y + 2, Talphaeta_12); deleteMatrix(N_y + 2, Talphaeta_13); deleteMatrix(N_y + 2, Talphaeta_14); deleteMatrix(N_y + 2, Talphaeta_15);
	deleteMatrix(N_y + 2, Talphaeta_21); deleteMatrix(N_y + 2, Talphaeta_22); deleteMatrix(N_y + 2, Talphaeta_23); deleteMatrix(N_y + 2, Talphaeta_24); deleteMatrix(N_y + 2, Talphaeta_25);
	deleteMatrix(N_y + 2, Talphaeta_31); deleteMatrix(N_y + 2, Talphaeta_32); deleteMatrix(N_y + 2, Talphaeta_33); deleteMatrix(N_y + 2, Talphaeta_34); deleteMatrix(N_y + 2, Talphaeta_35);
	deleteMatrix(N_y + 2, Talphaeta_41); deleteMatrix(N_y + 2, Talphaeta_42); deleteMatrix(N_y + 2, Talphaeta_43); deleteMatrix(N_y + 2, Talphaeta_44); deleteMatrix(N_y + 2, Talphaeta_45);
	deleteMatrix(N_y + 2, Talphaeta_51); deleteMatrix(N_y + 2, Talphaeta_52); deleteMatrix(N_y + 2, Talphaeta_53); deleteMatrix(N_y + 2, Talphaeta_54); deleteMatrix(N_y + 2, Talphaeta_55);

	deleteMatrix(N_y + 2, GQ_xi_1); deleteMatrix(N_y + 2, GQ_xi_eta); deleteMatrix(N_y + 2, GQ_xi_xi);
	deleteMatrix(N_y + 2, GQ_eta_1); deleteMatrix(N_y + 2, GQ_eta_xi); deleteMatrix(N_y + 2, GQ_eta_eta);
}

void Diffusion_r(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double kappa_alpha, double alpha_d, double beta_d, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5, double h_x, double h_y)
{
	double** Talpha_jx1; createMatrix(N_y, N_x + 3, Talpha_jx1); double** Talpha_jx2; createMatrix(N_y, N_x + 3, Talpha_jx2);
	double** Talpha_jx3; createMatrix(N_y, N_x + 3, Talpha_jx3); double** Talpha_jx4; createMatrix(N_y, N_x + 3, Talpha_jx4);
	double** Talpha_jx5; createMatrix(N_y, N_x + 3, Talpha_jx5);

	double** Talpha_jy1; createMatrix(N_y + 3, N_x, Talpha_jy1); double** Talpha_jy2; createMatrix(N_y + 3, N_x, Talpha_jy2);
	double** Talpha_jy3; createMatrix(N_y + 3, N_x, Talpha_jy3); double** Talpha_jy4; createMatrix(N_y + 3, N_x, Talpha_jy4);
	double** Talpha_jy5; createMatrix(N_y + 3, N_x, Talpha_jy5);

	double** Talphaxi_a1; createMatrix(N_y, N_x + 3, Talphaxi_a1); double** Talphaxi_a2; createMatrix(N_y, N_x + 3, Talphaxi_a2);
	double** Talphaxi_a3; createMatrix(N_y, N_x + 3, Talphaxi_a3); double** Talphaxi_a4; createMatrix(N_y, N_x + 3, Talphaxi_a4);
	double** Talphaxi_a5; createMatrix(N_y, N_x + 3, Talphaxi_a5);

	double** Talphaeta_a1; createMatrix(N_y + 3, N_x, Talphaeta_a1); double** Talphaeta_a2; createMatrix(N_y + 3, N_x, Talphaeta_a2);
	double** Talphaeta_a3; createMatrix(N_y + 3, N_x, Talphaeta_a3); double** Talphaeta_a4; createMatrix(N_y + 3, N_x, Talphaeta_a4);
	double** Talphaeta_a5; createMatrix(N_y + 3, N_x, Talphaeta_a5);

	double** Talphaxi2_j1; createMatrix(N_y, N_x + 3, Talphaxi2_j1); double** Talphaxi2_j2; createMatrix(N_y, N_x + 3, Talphaxi2_j2);
	double** Talphaxi2_j3; createMatrix(N_y, N_x + 3, Talphaxi2_j3); double** Talphaxi2_j4; createMatrix(N_y, N_x + 3, Talphaxi2_j4);
	double** Talphaxi2_j5; createMatrix(N_y, N_x + 3, Talphaxi2_j5);

	double** Talphaeta2_j1; createMatrix(N_y + 3, N_x, Talphaeta2_j1); double** Talphaeta2_j2; createMatrix(N_y + 3, N_x, Talphaeta2_j2);
	double** Talphaeta2_j3; createMatrix(N_y + 3, N_x, Talphaeta2_j3); double** Talphaeta2_j4; createMatrix(N_y + 3, N_x, Talphaeta2_j4);
	double** Talphaeta2_j5; createMatrix(N_y + 3, N_x, Talphaeta2_j5);

	double** Talphaxi_11; createMatrix(N_y + 2, N_x + 2, Talphaxi_11); double** Talphaxi_12; createMatrix(N_y + 2, N_x + 2, Talphaxi_12);
	double** Talphaxi_13; createMatrix(N_y + 2, N_x + 2, Talphaxi_13); double** Talphaxi_14; createMatrix(N_y + 2, N_x + 2, Talphaxi_14);
	double** Talphaxi_15; createMatrix(N_y + 2, N_x + 2, Talphaxi_15);

	double** Talphaxi_21; createMatrix(N_y + 2, N_x + 2, Talphaxi_21); double** Talphaxi_22; createMatrix(N_y + 2, N_x + 2, Talphaxi_22);
	double** Talphaxi_23; createMatrix(N_y + 2, N_x + 2, Talphaxi_23); double** Talphaxi_24; createMatrix(N_y + 2, N_x + 2, Talphaxi_24);
	double** Talphaxi_25; createMatrix(N_y + 2, N_x + 2, Talphaxi_25);

	double** Talphaxi_31; createMatrix(N_y + 2, N_x + 2, Talphaxi_31); double** Talphaxi_32; createMatrix(N_y + 2, N_x + 2, Talphaxi_32);
	double** Talphaxi_33; createMatrix(N_y + 2, N_x + 2, Talphaxi_33); double** Talphaxi_34; createMatrix(N_y + 2, N_x + 2, Talphaxi_34);
	double** Talphaxi_35; createMatrix(N_y + 2, N_x + 2, Talphaxi_35);

	double** Talphaxi_41; createMatrix(N_y + 2, N_x + 2, Talphaxi_41); double** Talphaxi_42; createMatrix(N_y + 2, N_x + 2, Talphaxi_42);
	double** Talphaxi_43; createMatrix(N_y + 2, N_x + 2, Talphaxi_43); double** Talphaxi_44; createMatrix(N_y + 2, N_x + 2, Talphaxi_44);
	double** Talphaxi_45; createMatrix(N_y + 2, N_x + 2, Talphaxi_45);

	double** Talphaxi_51; createMatrix(N_y + 2, N_x + 2, Talphaxi_51); double** Talphaxi_52; createMatrix(N_y + 2, N_x + 2, Talphaxi_52);
	double** Talphaxi_53; createMatrix(N_y + 2, N_x + 2, Talphaxi_53); double** Talphaxi_54; createMatrix(N_y + 2, N_x + 2, Talphaxi_54);
	double** Talphaxi_55; createMatrix(N_y + 2, N_x + 2, Talphaxi_55);

	double** Talphaeta_11; createMatrix(N_y + 2, N_x + 2, Talphaeta_11); double** Talphaeta_12; createMatrix(N_y + 2, N_x + 2, Talphaeta_12);
	double** Talphaeta_13; createMatrix(N_y + 2, N_x + 2, Talphaeta_13); double** Talphaeta_14; createMatrix(N_y + 2, N_x + 2, Talphaeta_14);
	double** Talphaeta_15; createMatrix(N_y + 2, N_x + 2, Talphaeta_15);

	double** Talphaeta_21; createMatrix(N_y + 2, N_x + 2, Talphaeta_21); double** Talphaeta_22; createMatrix(N_y + 2, N_x + 2, Talphaeta_22);
	double** Talphaeta_23; createMatrix(N_y + 2, N_x + 2, Talphaeta_23); double** Talphaeta_24; createMatrix(N_y + 2, N_x + 2, Talphaeta_24);
	double** Talphaeta_25; createMatrix(N_y + 2, N_x + 2, Talphaeta_25);

	double** Talphaeta_31; createMatrix(N_y + 2, N_x + 2, Talphaeta_31); double** Talphaeta_32; createMatrix(N_y + 2, N_x + 2, Talphaeta_32);
	double** Talphaeta_33; createMatrix(N_y + 2, N_x + 2, Talphaeta_33); double** Talphaeta_34; createMatrix(N_y + 2, N_x + 2, Talphaeta_34);
	double** Talphaeta_35; createMatrix(N_y + 2, N_x + 2, Talphaeta_35);

	double** Talphaeta_41; createMatrix(N_y + 2, N_x + 2, Talphaeta_41); double** Talphaeta_42; createMatrix(N_y + 2, N_x + 2, Talphaeta_42);
	double** Talphaeta_43; createMatrix(N_y + 2, N_x + 2, Talphaeta_43); double** Talphaeta_44; createMatrix(N_y + 2, N_x + 2, Talphaeta_44);
	double** Talphaeta_45; createMatrix(N_y + 2, N_x + 2, Talphaeta_45);

	double** Talphaeta_51; createMatrix(N_y + 2, N_x + 2, Talphaeta_51); double** Talphaeta_52; createMatrix(N_y + 2, N_x + 2, Talphaeta_52);
	double** Talphaeta_53; createMatrix(N_y + 2, N_x + 2, Talphaeta_53); double** Talphaeta_54; createMatrix(N_y + 2, N_x + 2, Talphaeta_54);
	double** Talphaeta_55; createMatrix(N_y + 2, N_x + 2, Talphaeta_55);


	// needed modified
	Tr4_jumpx(Talpha_jx1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpx(Talpha_jx2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpx(Talpha_jx3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpx(Talpha_jx4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpx(Talpha_jx5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Tr4_jumpy(Talpha_jy1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpy(Talpha_jy2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpy(Talpha_jy3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpy(Talpha_jy4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4_jumpy(Talpha_jy5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Tr4xi_ave(Talphaxi_a1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi_ave(Talphaxi_a2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi_ave(Talphaxi_a3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi_ave(Talphaxi_a4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi_ave(Talphaxi_a5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Tr4eta_ave(Talphaeta_a1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta_ave(Talphaeta_a2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta_ave(Talphaeta_a3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta_ave(Talphaeta_a4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta_ave(Talphaeta_a5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Tr4xi2_jump(Talphaxi2_j1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi2_jump(Talphaxi2_j2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi2_jump(Talphaxi2_j3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi2_jump(Talphaxi2_j4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4xi2_jump(Talphaxi2_j5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	Tr4eta2_jump(Talphaeta2_j1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta2_jump(Talphaeta2_j2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta2_jump(Talphaeta2_j3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta2_jump(Talphaeta2_j4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	Tr4eta2_jump(Talphaeta2_j5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4xi(Talphaxi_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4xi(Talphaxi_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4xi(Talphaxi_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4xi(Talphaxi_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4xi(Talphaxi_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4xi(Talphaxi_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4eta(Talphaeta_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4eta(Talphaeta_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4eta(Talphaeta_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4eta(Talphaeta_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	compute_Tr4eta(Talphaeta_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	compute_Tr4eta(Talphaeta_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	double** GQ_xi_1; createMatrix(N_y + 2, N_x + 2, GQ_xi_1); double** GQ_xi_eta; createMatrix(N_y + 2, N_x + 2, GQ_xi_eta); double** GQ_xi_xi; createMatrix(N_y + 2, N_x + 2, GQ_xi_xi);
	double** GQ_eta_1; createMatrix(N_y + 2, N_x + 2, GQ_eta_1); double** GQ_eta_xi; createMatrix(N_y + 2, N_x + 2, GQ_eta_xi); double** GQ_eta_eta; createMatrix(N_y + 2, N_x + 2, GQ_eta_eta);

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
			GQ_xi_1[i][j] = -1.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] + w_12 * Talphaxi_12[i][j] + w_13 * Talphaxi_13[i][j] + w_14 * Talphaxi_14[i][j] + w_15 * Talphaxi_15[i][j] + w_21 * Talphaxi_21[i][j] + w_22 * Talphaxi_22[i][j] + w_23 * Talphaxi_23[i][j] + w_24 * Talphaxi_24[i][j] + w_25 * Talphaxi_25[i][j] + w_31 * Talphaxi_31[i][j] + w_32 * Talphaxi_32[i][j] + w_33 * Talphaxi_33[i][j] + w_34 * Talphaxi_34[i][j] + w_35 * Talphaxi_35[i][j] + w_41 * Talphaxi_41[i][j] + w_42 * Talphaxi_42[i][j] + w_43 * Talphaxi_43[i][j] + w_44 * Talphaxi_44[i][j] + w_45 * Talphaxi_45[i][j] + w_51 * Talphaxi_51[i][j] + w_52 * Talphaxi_52[i][j] + w_53 * Talphaxi_53[i][j] + w_54 * Talphaxi_54[i][j] + w_55 * Talphaxi_55[i][j]);
			GQ_xi_eta[i][j] = -1.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] * GQxi_1 + w_12 * Talphaxi_12[i][j] * GQxi_2 + w_13 * Talphaxi_13[i][j] * GQxi_3 + w_14 * Talphaxi_14[i][j] * GQxi_4 + w_15 * Talphaxi_15[i][j] * GQxi_5 + w_21 * Talphaxi_21[i][j] * GQxi_1 + w_22 * Talphaxi_22[i][j] * GQxi_2 + w_23 * Talphaxi_23[i][j] * GQxi_3 + w_24 * Talphaxi_24[i][j] * GQxi_4 + w_25 * Talphaxi_25[i][j] * GQxi_5 + w_31 * Talphaxi_31[i][j] * GQxi_1 + w_32 * Talphaxi_32[i][j] * GQxi_2 + w_33 * Talphaxi_33[i][j] * GQxi_3 + w_34 * Talphaxi_34[i][j] * GQxi_4 + w_35 * Talphaxi_35[i][j] * GQxi_5 + w_41 * Talphaxi_41[i][j] * GQxi_1 + w_42 * Talphaxi_42[i][j] * GQxi_2 + w_43 * Talphaxi_43[i][j] * GQxi_3 + w_44 * Talphaxi_44[i][j] * GQxi_4 + w_45 * Talphaxi_45[i][j] * GQxi_5 + w_51 * Talphaxi_51[i][j] * GQxi_1 + w_52 * Talphaxi_52[i][j] * GQxi_2 + w_53 * Talphaxi_53[i][j] * GQxi_3 + w_54 * Talphaxi_54[i][j] * GQxi_4 + w_55 * Talphaxi_55[i][j] * GQxi_5);
			GQ_xi_xi[i][j] = -3.0 * kappa_alpha * h_y / h_x * (w_11 * Talphaxi_11[i][j] * GQxi_1 + w_12 * Talphaxi_12[i][j] * GQxi_1 + w_13 * Talphaxi_13[i][j] * GQxi_1 + w_14 * Talphaxi_14[i][j] * GQxi_1 + w_15 * Talphaxi_15[i][j] * GQxi_1 + w_21 * Talphaxi_21[i][j] * GQxi_2 + w_22 * Talphaxi_22[i][j] * GQxi_2 + w_23 * Talphaxi_23[i][j] * GQxi_2 + w_24 * Talphaxi_24[i][j] * GQxi_2 + w_25 * Talphaxi_25[i][j] * GQxi_2 + w_31 * Talphaxi_31[i][j] * GQxi_3 + w_32 * Talphaxi_32[i][j] * GQxi_3 + w_33 * Talphaxi_33[i][j] * GQxi_3 + w_34 * Talphaxi_34[i][j] * GQxi_3 + w_35 * Talphaxi_35[i][j] * GQxi_3 + w_41 * Talphaxi_41[i][j] * GQxi_4 + w_42 * Talphaxi_42[i][j] * GQxi_4 + w_43 * Talphaxi_43[i][j] * GQxi_4 + w_44 * Talphaxi_44[i][j] * GQxi_4 + w_45 * Talphaxi_45[i][j] * GQxi_4 + w_51 * Talphaxi_51[i][j] * GQxi_5 + w_52 * Talphaxi_52[i][j] * GQxi_5 + w_53 * Talphaxi_53[i][j] * GQxi_5 + w_54 * Talphaxi_54[i][j] * GQxi_5 + w_55 * Talphaxi_55[i][j] * GQxi_5);


			GQ_eta_1[i][j] = -1.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] + w_12 * Talphaeta_12[i][j] + w_13 * Talphaeta_13[i][j] + w_14 * Talphaeta_14[i][j] + w_15 * Talphaeta_15[i][j] + w_21 * Talphaeta_21[i][j] + w_22 * Talphaeta_22[i][j] + w_23 * Talphaeta_23[i][j] + w_24 * Talphaeta_24[i][j] + w_25 * Talphaeta_25[i][j] + w_31 * Talphaeta_31[i][j] + w_32 * Talphaeta_32[i][j] + w_33 * Talphaeta_33[i][j] + w_34 * Talphaeta_34[i][j] + w_35 * Talphaeta_35[i][j] + w_41 * Talphaeta_41[i][j] + w_42 * Talphaeta_42[i][j] + w_43 * Talphaeta_43[i][j] + w_44 * Talphaeta_44[i][j] + w_45 * Talphaeta_45[i][j] + w_51 * Talphaeta_51[i][j] + w_52 * Talphaeta_52[i][j] + w_53 * Talphaeta_53[i][j] + w_54 * Talphaeta_54[i][j] + w_55 * Talphaeta_55[i][j]);
			GQ_eta_xi[i][j] = -1.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] * GQxi_1 + w_12 * Talphaeta_12[i][j] * GQxi_1 + w_13 * Talphaeta_13[i][j] * GQxi_1 + w_14 * Talphaeta_14[i][j] * GQxi_1 + w_15 * Talphaeta_15[i][j] * GQxi_1 + w_21 * Talphaeta_21[i][j] * GQxi_2 + w_22 * Talphaeta_22[i][j] * GQxi_2 + w_23 * Talphaeta_23[i][j] * GQxi_2 + w_24 * Talphaeta_24[i][j] * GQxi_2 + w_25 * Talphaeta_25[i][j] * GQxi_2 + w_31 * Talphaeta_31[i][j] * GQxi_3 + w_32 * Talphaeta_32[i][j] * GQxi_3 + w_33 * Talphaeta_33[i][j] * GQxi_3 + w_34 * Talphaeta_34[i][j] * GQxi_3 + w_35 * Talphaeta_35[i][j] * GQxi_3 + w_41 * Talphaeta_41[i][j] * GQxi_4 + w_42 * Talphaeta_42[i][j] * GQxi_4 + w_43 * Talphaeta_43[i][j] * GQxi_4 + w_44 * Talphaeta_44[i][j] * GQxi_4 + w_45 * Talphaeta_45[i][j] * GQxi_4 + w_51 * Talphaeta_51[i][j] * GQxi_5 + w_52 * Talphaeta_52[i][j] * GQxi_5 + w_53 * Talphaeta_53[i][j] * GQxi_5 + w_54 * Talphaeta_54[i][j] * GQxi_5 + w_55 * Talphaeta_55[i][j] * GQxi_5);
			GQ_eta_eta[i][j] = -3.0 * kappa_alpha * h_x / h_y * (w_11 * Talphaeta_11[i][j] * GQxi_1 + w_12 * Talphaeta_12[i][j] * GQxi_2 + w_13 * Talphaeta_13[i][j] * GQxi_3 + w_14 * Talphaeta_14[i][j] * GQxi_4 + w_15 * Talphaeta_15[i][j] * GQxi_5 + w_21 * Talphaeta_21[i][j] * GQxi_1 + w_22 * Talphaeta_22[i][j] * GQxi_2 + w_23 * Talphaeta_23[i][j] * GQxi_3 + w_24 * Talphaeta_24[i][j] * GQxi_4 + w_25 * Talphaeta_25[i][j] * GQxi_5 + w_31 * Talphaeta_31[i][j] * GQxi_1 + w_32 * Talphaeta_32[i][j] * GQxi_2 + w_33 * Talphaeta_33[i][j] * GQxi_3 + w_34 * Talphaeta_34[i][j] * GQxi_4 + w_35 * Talphaeta_35[i][j] * GQxi_5 + w_41 * Talphaeta_41[i][j] * GQxi_1 + w_42 * Talphaeta_42[i][j] * GQxi_2 + w_43 * Talphaeta_43[i][j] * GQxi_3 + w_44 * Talphaeta_44[i][j] * GQxi_4 + w_45 * Talphaeta_45[i][j] * GQxi_5 + w_51 * Talphaeta_51[i][j] * GQxi_1 + w_52 * Talphaeta_52[i][j] * GQxi_2 + w_53 * Talphaeta_53[i][j] * GQxi_3 + w_54 * Talphaeta_54[i][j] * GQxi_4 + w_55 * Talphaeta_55[i][j] * GQxi_5);

		}
	}

	double** GQ_Talpha_jx; createMatrix(N_y, N_x + 3, GQ_Talpha_jx); double** GQ_Talphaxi_a; createMatrix(N_y, N_x + 3, GQ_Talphaxi_a); double** GQ_Talphaxi2_j; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_j);
	double** GQ_Talpha_jxeta; createMatrix(N_y, N_x + 3, GQ_Talpha_jxeta); double** GQ_Talphaxi_aeta; createMatrix(N_y, N_x + 3, GQ_Talphaxi_aeta); double** GQ_Talphaxi2_jeta; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_jeta);
	double** GQ_Talpha_jxeta2; createMatrix(N_y, N_x + 3, GQ_Talpha_jxeta2); double** GQ_Talphaxi_aeta2; createMatrix(N_y, N_x + 3, GQ_Talphaxi_aeta2); double** GQ_Talphaxi2_jeta2; createMatrix(N_y, N_x + 3, GQ_Talphaxi2_jeta2);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x + 3; j++)
		{
			GQ_Talpha_jx[i][j] = weight_1 * Talpha_jx1[i][j] + weight_2 * Talpha_jx2[i][j] + weight_3 * Talpha_jx3[i][j] + weight_4 * Talpha_jx4[i][j] + weight_5 * Talpha_jx5[i][j];
			GQ_Talphaxi_a[i][j] = weight_1 * Talphaxi_a1[i][j] + weight_2 * Talphaxi_a2[i][j] + weight_3 * Talphaxi_a3[i][j] + weight_4 * Talphaxi_a4[i][j] + weight_5 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_j[i][j] = weight_1 * Talphaxi2_j1[i][j] + weight_2 * Talphaxi2_j2[i][j] + weight_3 * Talphaxi2_j3[i][j] + weight_4 * Talphaxi2_j4[i][j] + weight_5 * Talphaxi2_j5[i][j];

			GQ_Talpha_jxeta[i][j] = weight_1 * GQxi_1 * Talpha_jx1[i][j] + weight_2 * GQxi_2 * Talpha_jx2[i][j] + weight_3 * GQxi_3 * Talpha_jx3[i][j] + weight_4 * GQxi_4 * Talpha_jx4[i][j] + weight_5 * GQxi_5 * Talpha_jx5[i][j];
			GQ_Talphaxi_aeta[i][j] = weight_1 * GQxi_1 * Talphaxi_a1[i][j] + weight_2 * GQxi_2 * Talphaxi_a2[i][j] + weight_3 * GQxi_3 * Talphaxi_a3[i][j] + weight_4 * GQxi_4 * Talphaxi_a4[i][j] + weight_5 * GQxi_5 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_jeta[i][j] = weight_1 * GQxi_1 * Talphaxi2_j1[i][j] + weight_2 * GQxi_2 * Talphaxi2_j2[i][j] + weight_3 * GQxi_3 * Talphaxi2_j3[i][j] + weight_4 * GQxi_4 * Talphaxi2_j4[i][j] + weight_5 * GQxi_5 * Talphaxi2_j5[i][j];

			GQ_Talpha_jxeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talpha_jx1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talpha_jx2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talpha_jx3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talpha_jx4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talpha_jx5[i][j];
			GQ_Talphaxi_aeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaxi_a1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaxi_a2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaxi_a3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaxi_a4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaxi_a5[i][j];
			GQ_Talphaxi2_jeta2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaxi2_j1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaxi2_j2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaxi2_j3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaxi2_j4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaxi2_j5[i][j];
		}
	}

	double** GQ_Talpha_jy; createMatrix(N_y + 3, N_x, GQ_Talpha_jy); double** GQ_Talphaeta_a; createMatrix(N_y + 3, N_x, GQ_Talphaeta_a); double** GQ_Talphaeta2_j; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_j);
	double** GQ_Talpha_jyxi; createMatrix(N_y + 3, N_x, GQ_Talpha_jyxi); double** GQ_Talphaeta_axi; createMatrix(N_y + 3, N_x, GQ_Talphaeta_axi); double** GQ_Talphaeta2_jxi; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_jxi);
	double** GQ_Talpha_jyxi2; createMatrix(N_y + 3, N_x, GQ_Talpha_jyxi2); double** GQ_Talphaeta_axi2; createMatrix(N_y + 3, N_x, GQ_Talphaeta_axi2); double** GQ_Talphaeta2_jxi2; createMatrix(N_y + 3, N_x, GQ_Talphaeta2_jxi2);

	for (int j = 0; j < N_x; j++)
	{
		for (int i = 0; i < N_y + 3; i++)
		{
			GQ_Talpha_jy[i][j] = weight_1 * Talpha_jy1[i][j] + weight_2 * Talpha_jy2[i][j] + weight_3 * Talpha_jy3[i][j] + weight_4 * Talpha_jy4[i][j] + weight_5 * Talpha_jy5[i][j];
			GQ_Talphaeta_a[i][j] = weight_1 * Talphaeta_a1[i][j] + weight_2 * Talphaeta_a2[i][j] + weight_3 * Talphaeta_a3[i][j] + weight_4 * Talphaeta_a4[i][j] + weight_5 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_j[i][j] = weight_1 * Talphaeta2_j1[i][j] + weight_2 * Talphaeta2_j2[i][j] + weight_3 * Talphaeta2_j3[i][j] + weight_4 * Talphaeta2_j4[i][j] + weight_5 * Talphaeta2_j5[i][j];

			GQ_Talpha_jyxi[i][j] = weight_1 * GQxi_1 * Talpha_jy1[i][j] + weight_2 * GQxi_2 * Talpha_jy2[i][j] + weight_3 * GQxi_3 * Talpha_jy3[i][j] + weight_4 * GQxi_4 * Talpha_jy4[i][j] + weight_5 * GQxi_5 * Talpha_jy5[i][j];
			GQ_Talphaeta_axi[i][j] = weight_1 * GQxi_1 * Talphaeta_a1[i][j] + weight_2 * GQxi_2 * Talphaeta_a2[i][j] + weight_3 * GQxi_3 * Talphaeta_a3[i][j] + weight_4 * GQxi_4 * Talphaeta_a4[i][j] + weight_5 * GQxi_5 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_jxi[i][j] = weight_1 * GQxi_1 * Talphaeta2_j1[i][j] + weight_2 * GQxi_2 * Talphaeta2_j2[i][j] + weight_3 * GQxi_3 * Talphaeta2_j3[i][j] + weight_4 * GQxi_4 * Talphaeta2_j4[i][j] + weight_5 * GQxi_5 * Talphaeta2_j5[i][j];

			GQ_Talpha_jyxi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talpha_jy1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talpha_jy2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talpha_jy3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talpha_jy4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talpha_jy5[i][j];
			GQ_Talphaeta_axi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaeta_a1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaeta_a2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaeta_a3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaeta_a4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaeta_a5[i][j];
			GQ_Talphaeta2_jxi2[i][j] = weight_1 * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 * Talphaeta2_j1[i][j] + weight_2 * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 * Talphaeta2_j2[i][j] + weight_3 * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 * Talphaeta2_j3[i][j] + weight_4 * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 * Talphaeta2_j4[i][j] + weight_5 * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0 * Talphaeta2_j5[i][j];
		}
	}


	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			U_0[i][j] = 0.0 + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] - GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] - GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] - GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] - GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] - GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] - GQ_Talphaeta2_j[i][j - 1]);
			U_1[i][j] = GQ_xi_1[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] + GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi[i + 1][j - 1] - GQ_Talphaeta_axi[i][j - 1]) + (alpha_d - kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] + GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi[i + 1][j - 1] - GQ_Talpha_jyxi[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] + GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi[i + 1][j - 1] - GQ_Talphaeta2_jxi[i][j - 1]);
			U_2[i][j] = GQ_eta_1[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta[i - 1][j + 1] - GQ_Talphaxi_aeta[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] + GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta[i - 1][j + 1] - GQ_Talpha_jxeta[i - 1][j]) + (alpha_d - kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] + GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta[i - 1][j + 1] - GQ_Talphaxi2_jeta[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] + GQ_Talphaeta2_j[i][j - 1]);
			U_3[i][j] = GQ_xi_eta[i][j] + GQ_eta_xi[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta[i - 1][j + 1] + GQ_Talphaxi_aeta[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi[i + 1][j - 1] + GQ_Talphaeta_axi[i][j - 1]) + (alpha_d - kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta[i - 1][j + 1] + GQ_Talpha_jxeta[i - 1][j]) + (alpha_d - kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi[i + 1][j - 1] + GQ_Talpha_jyxi[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta[i - 1][j + 1] + GQ_Talphaxi2_jeta[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi[i + 1][j - 1] + GQ_Talphaeta2_jxi[i][j - 1]);
			U_4[i][j] = GQ_xi_xi[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_a[i - 1][j + 1] - GQ_Talphaxi_a[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_axi2[i + 1][j - 1] - GQ_Talphaeta_axi2[i][j - 1]) + (alpha_d - 3.0 * kappa_alpha) * h_y / (2.0 * h_x) * (GQ_Talpha_jx[i - 1][j + 1] - GQ_Talpha_jx[i - 1][j]) + alpha_d * h_x / (2.0 * h_y) * (GQ_Talpha_jyxi2[i + 1][j - 1] - GQ_Talpha_jyxi2[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_j[i - 1][j + 1] - GQ_Talphaxi2_j[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_jxi2[i + 1][j - 1] - GQ_Talphaeta2_jxi2[i][j - 1]);
			U_5[i][j] = GQ_eta_eta[i][j] + kappa_alpha * h_y / h_x * (GQ_Talphaxi_aeta2[i - 1][j + 1] - GQ_Talphaxi_aeta2[i - 1][j]) + kappa_alpha * h_x / h_y * (GQ_Talphaeta_a[i + 1][j - 1] - GQ_Talphaeta_a[i][j - 1]) + alpha_d * h_y / (2.0 * h_x) * (GQ_Talpha_jxeta2[i - 1][j + 1] - GQ_Talpha_jxeta2[i - 1][j]) + (alpha_d - 3.0 * kappa_alpha) * h_x / (2.0 * h_y) * (GQ_Talpha_jy[i + 1][j - 1] - GQ_Talpha_jy[i][j - 1]) + 2.0 * beta_d * kappa_alpha * h_y / h_x * (GQ_Talphaxi2_jeta2[i - 1][j + 1] - GQ_Talphaxi2_jeta2[i - 1][j]) + 2.0 * beta_d * kappa_alpha * h_x / h_y * (GQ_Talphaeta2_j[i + 1][j - 1] - GQ_Talphaeta2_j[i][j - 1]);
		}
	}

	deleteMatrix(N_y, Talpha_jx1); deleteMatrix(N_y, Talpha_jx2); deleteMatrix(N_y, Talpha_jx3); deleteMatrix(N_y, Talpha_jx4); deleteMatrix(N_y, Talpha_jx5);
	deleteMatrix(N_y + 3, Talpha_jy1); deleteMatrix(N_y + 3, Talpha_jy2); deleteMatrix(N_y + 3, Talpha_jy3); deleteMatrix(N_y + 3, Talpha_jy4); deleteMatrix(N_y + 3, Talpha_jy5);
	deleteMatrix(N_y, Talphaxi_a1); deleteMatrix(N_y, Talphaxi_a2); deleteMatrix(N_y, Talphaxi_a3); deleteMatrix(N_y, Talphaxi_a4); deleteMatrix(N_y, Talphaxi_a5);
	deleteMatrix(N_y + 3, Talphaeta_a1); deleteMatrix(N_y + 3, Talphaeta_a2); deleteMatrix(N_y + 3, Talphaeta_a3); deleteMatrix(N_y + 3, Talphaeta_a4); deleteMatrix(N_y + 3, Talphaeta_a5);
	deleteMatrix(N_y, Talphaxi2_j1); deleteMatrix(N_y, Talphaxi2_j2); deleteMatrix(N_y, Talphaxi2_j3); deleteMatrix(N_y, Talphaxi2_j4); deleteMatrix(N_y, Talphaxi2_j5);
	deleteMatrix(N_y + 3, Talphaeta2_j1); deleteMatrix(N_y + 3, Talphaeta2_j2); deleteMatrix(N_y + 3, Talphaeta2_j3); deleteMatrix(N_y + 3, Talphaeta2_j4); deleteMatrix(N_y + 3, Talphaeta2_j5);

	deleteMatrix(N_y, GQ_Talpha_jx); deleteMatrix(N_y, GQ_Talphaxi_a); deleteMatrix(N_y, GQ_Talphaxi2_j);
	deleteMatrix(N_y, GQ_Talpha_jxeta); deleteMatrix(N_y, GQ_Talphaxi_aeta); deleteMatrix(N_y, GQ_Talphaxi2_jeta);
	deleteMatrix(N_y, GQ_Talpha_jxeta2); deleteMatrix(N_y, GQ_Talphaxi_aeta2); deleteMatrix(N_y, GQ_Talphaxi2_jeta2);

	deleteMatrix(N_y + 3, GQ_Talpha_jy); deleteMatrix(N_y + 3, GQ_Talphaeta_a); deleteMatrix(N_y + 3, GQ_Talphaeta2_j);
	deleteMatrix(N_y + 3, GQ_Talpha_jyxi); deleteMatrix(N_y + 3, GQ_Talphaeta_axi); deleteMatrix(N_y + 3, GQ_Talphaeta2_jxi);
	deleteMatrix(N_y + 3, GQ_Talpha_jyxi2); deleteMatrix(N_y + 3, GQ_Talphaeta_axi2); deleteMatrix(N_y + 3, GQ_Talphaeta2_jxi2);

	deleteMatrix(N_y + 2, Talphaxi_11); deleteMatrix(N_y + 2, Talphaxi_12); deleteMatrix(N_y + 2, Talphaxi_13); deleteMatrix(N_y + 2, Talphaxi_14); deleteMatrix(N_y + 2, Talphaxi_15);
	deleteMatrix(N_y + 2, Talphaxi_21); deleteMatrix(N_y + 2, Talphaxi_22); deleteMatrix(N_y + 2, Talphaxi_23); deleteMatrix(N_y + 2, Talphaxi_24); deleteMatrix(N_y + 2, Talphaxi_25);
	deleteMatrix(N_y + 2, Talphaxi_31); deleteMatrix(N_y + 2, Talphaxi_32); deleteMatrix(N_y + 2, Talphaxi_33); deleteMatrix(N_y + 2, Talphaxi_34); deleteMatrix(N_y + 2, Talphaxi_35);
	deleteMatrix(N_y + 2, Talphaxi_41); deleteMatrix(N_y + 2, Talphaxi_42); deleteMatrix(N_y + 2, Talphaxi_43); deleteMatrix(N_y + 2, Talphaxi_44); deleteMatrix(N_y + 2, Talphaxi_45);
	deleteMatrix(N_y + 2, Talphaxi_51); deleteMatrix(N_y + 2, Talphaxi_52); deleteMatrix(N_y + 2, Talphaxi_53); deleteMatrix(N_y + 2, Talphaxi_54); deleteMatrix(N_y + 2, Talphaxi_55);

	deleteMatrix(N_y + 2, Talphaeta_11); deleteMatrix(N_y + 2, Talphaeta_12); deleteMatrix(N_y + 2, Talphaeta_13); deleteMatrix(N_y + 2, Talphaeta_14); deleteMatrix(N_y + 2, Talphaeta_15);
	deleteMatrix(N_y + 2, Talphaeta_21); deleteMatrix(N_y + 2, Talphaeta_22); deleteMatrix(N_y + 2, Talphaeta_23); deleteMatrix(N_y + 2, Talphaeta_24); deleteMatrix(N_y + 2, Talphaeta_25);
	deleteMatrix(N_y + 2, Talphaeta_31); deleteMatrix(N_y + 2, Talphaeta_32); deleteMatrix(N_y + 2, Talphaeta_33); deleteMatrix(N_y + 2, Talphaeta_34); deleteMatrix(N_y + 2, Talphaeta_35);
	deleteMatrix(N_y + 2, Talphaeta_41); deleteMatrix(N_y + 2, Talphaeta_42); deleteMatrix(N_y + 2, Talphaeta_43); deleteMatrix(N_y + 2, Talphaeta_44); deleteMatrix(N_y + 2, Talphaeta_45);
	deleteMatrix(N_y + 2, Talphaeta_51); deleteMatrix(N_y + 2, Talphaeta_52); deleteMatrix(N_y + 2, Talphaeta_53); deleteMatrix(N_y + 2, Talphaeta_54); deleteMatrix(N_y + 2, Talphaeta_55);

	deleteMatrix(N_y + 2, GQ_xi_1); deleteMatrix(N_y + 2, GQ_xi_eta); deleteMatrix(N_y + 2, GQ_xi_xi);
	deleteMatrix(N_y + 2, GQ_eta_1); deleteMatrix(N_y + 2, GQ_eta_xi); deleteMatrix(N_y + 2, GQ_eta_eta);
}



void space_d(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	convection_d(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac, F2_alphac, h_x, h_y);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j];
			U_1[i][j] = con_1[i][j];
			U_2[i][j] = con_2[i][j];
			U_3[i][j] = con_3[i][j];
			U_4[i][j] = con_4[i][j];
			U_5[i][j] = con_5[i][j];

		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);

}

void space_m1(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	convection_m1(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j];
			U_1[i][j] = con_1[i][j];
			U_2[i][j] = con_2[i][j];
			U_3[i][j] = con_3[i][j];
			U_4[i][j] = con_4[i][j];
			U_5[i][j] = con_5[i][j];

		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);
}

void space_m2(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	convection_m2(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j];
			U_1[i][j] = con_1[i][j];
			U_2[i][j] = con_2[i][j];
			U_3[i][j] = con_3[i][j];
			U_4[i][j] = con_4[i][j];
			U_5[i][j] = con_5[i][j];
		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);
}

void space_e(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t, int index_mode)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	double** Non_0; createMatrix(N_y + 2, N_x + 2, Non_0); double** Non_1; createMatrix(N_y + 2, N_x + 2, Non_1); double** Non_2; createMatrix(N_y + 2, N_x + 2, Non_2);
	double** Non_3; createMatrix(N_y + 2, N_x + 2, Non_3); double** Non_4; createMatrix(N_y + 2, N_x + 2, Non_4); double** Non_5; createMatrix(N_y + 2, N_x + 2, Non_5);

	double** Diff_0; createMatrix(N_y + 2, N_x + 2, Diff_0); double** Diff_1; createMatrix(N_y + 2, N_x + 2, Diff_1); double** Diff_2; createMatrix(N_y + 2, N_x + 2, Diff_2);
	double** Diff_3; createMatrix(N_y + 2, N_x + 2, Diff_3); double** Diff_4; createMatrix(N_y + 2, N_x + 2, Diff_4); double** Diff_5; createMatrix(N_y + 2, N_x + 2, Diff_5);

	double** Source_0; createMatrix(N_y + 2, N_x + 2, Source_0); double** Source_1; createMatrix(N_y + 2, N_x + 2, Source_1); double** Source_2; createMatrix(N_y + 2, N_x + 2, Source_2);
	double** Source_3; createMatrix(N_y + 2, N_x + 2, Source_3); double** Source_4; createMatrix(N_y + 2, N_x + 2, Source_4); double** Source_5; createMatrix(N_y + 2, N_x + 2, Source_5);

	double diff_alpha_e = 0.0;
	if (index_mode == 2)
	{
		diff_alpha_e = Coef_d2 * kappa_e;
	}
	else if (index_mode == 1)
	{
		diff_alpha_e = Coef_d1 * kappa_e;
	}
	else
	{
		diff_alpha_e = kappa_e;
	}

	convection_Ealpha(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, gamma_e, F1_alphac, F2_alphac, h_x, h_y);

	double** Source_i_0; createMatrix(N_y + 2, N_x + 2, Source_i_0); double** Source_i_1; createMatrix(N_y + 2, N_x + 2, Source_i_1);
	double** Source_i_2; createMatrix(N_y + 2, N_x + 2, Source_i_2); double** Source_i_3; createMatrix(N_y + 2, N_x + 2, Source_i_3);
	double** Source_i_4; createMatrix(N_y + 2, N_x + 2, Source_i_4); double** Source_i_5; createMatrix(N_y + 2, N_x + 2, Source_i_5);
	double** Source_r_0; createMatrix(N_y + 2, N_x + 2, Source_r_0); double** Source_r_1; createMatrix(N_y + 2, N_x + 2, Source_r_1);
	double** Source_r_2; createMatrix(N_y + 2, N_x + 2, Source_r_2); double** Source_r_3; createMatrix(N_y + 2, N_x + 2, Source_r_3);
	double** Source_r_4; createMatrix(N_y + 2, N_x + 2, Source_r_4); double** Source_r_5; createMatrix(N_y + 2, N_x + 2, Source_r_5);

	Source_i(Source_i_0, Source_i_1, Source_i_2, Source_i_3, Source_i_4, Source_i_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, h_x, h_y);
	Source_r(Source_r_0, Source_r_1, Source_r_2, Source_r_3, Source_r_4, Source_r_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);
	Source_e(Source_0, Source_1, Source_2, Source_3, Source_4, Source_5, Source_i_0, Source_i_1, Source_i_2, Source_i_3, Source_i_4, Source_i_5, Source_r_0, Source_r_1, Source_r_2, Source_r_3, Source_r_4, Source_r_5);

	deleteMatrix(N_y + 2, Source_i_0); deleteMatrix(N_y + 2, Source_i_1); deleteMatrix(N_y + 2, Source_i_2);
	deleteMatrix(N_y + 2, Source_i_3); deleteMatrix(N_y + 2, Source_i_4); deleteMatrix(N_y + 2, Source_i_5);
	deleteMatrix(N_y + 2, Source_r_0); deleteMatrix(N_y + 2, Source_r_1); deleteMatrix(N_y + 2, Source_r_2);
	deleteMatrix(N_y + 2, Source_r_3); deleteMatrix(N_y + 2, Source_r_4); deleteMatrix(N_y + 2, Source_r_5);

	Nonconservetive_Ne(Non_0, Non_1, Non_2, Non_3, Non_4, Non_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	Diffusion_alpha(Diff_0, Diff_1, Diff_2, Diff_3, Diff_4, Diff_5, C_ve, kappa_e, diff_alpha_e, diff_beta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, h_x, h_y);


	// Accuracy source part
	double** AS_11; createMatrix(N_y + 2, N_x + 2, AS_11); double** AS_12; createMatrix(N_y + 2, N_x + 2, AS_12); double** AS_13; createMatrix(N_y + 2, N_x + 2, AS_13); double** AS_14; createMatrix(N_y + 2, N_x + 2, AS_14); double** AS_15; createMatrix(N_y + 2, N_x + 2, AS_15);
	double** AS_21; createMatrix(N_y + 2, N_x + 2, AS_21); double** AS_22; createMatrix(N_y + 2, N_x + 2, AS_22); double** AS_23; createMatrix(N_y + 2, N_x + 2, AS_23); double** AS_24; createMatrix(N_y + 2, N_x + 2, AS_24); double** AS_25; createMatrix(N_y + 2, N_x + 2, AS_25);
	double** AS_31; createMatrix(N_y + 2, N_x + 2, AS_31); double** AS_32; createMatrix(N_y + 2, N_x + 2, AS_32); double** AS_33; createMatrix(N_y + 2, N_x + 2, AS_33); double** AS_34; createMatrix(N_y + 2, N_x + 2, AS_34); double** AS_35; createMatrix(N_y + 2, N_x + 2, AS_35);
	double** AS_41; createMatrix(N_y + 2, N_x + 2, AS_41); double** AS_42; createMatrix(N_y + 2, N_x + 2, AS_42); double** AS_43; createMatrix(N_y + 2, N_x + 2, AS_43); double** AS_44; createMatrix(N_y + 2, N_x + 2, AS_44); double** AS_45; createMatrix(N_y + 2, N_x + 2, AS_45);
	double** AS_51; createMatrix(N_y + 2, N_x + 2, AS_51); double** AS_52; createMatrix(N_y + 2, N_x + 2, AS_52); double** AS_53; createMatrix(N_y + 2, N_x + 2, AS_53); double** AS_54; createMatrix(N_y + 2, N_x + 2, AS_54); double** AS_55; createMatrix(N_y + 2, N_x + 2, AS_55);

	double** AS_GQ_0; createMatrix(N_y + 2, N_x + 2, AS_GQ_0); double** AS_GQ_1; createMatrix(N_y + 2, N_x + 2, AS_GQ_1);
	double** AS_GQ_2; createMatrix(N_y + 2, N_x + 2, AS_GQ_2); double** AS_GQ_3; createMatrix(N_y + 2, N_x + 2, AS_GQ_3);
	double** AS_GQ_4; createMatrix(N_y + 2, N_x + 2, AS_GQ_4); double** AS_GQ_5; createMatrix(N_y + 2, N_x + 2, AS_GQ_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			AS_11[i][j] = S3_ex(X_1[j], Y_1[i], t); AS_21[i][j] = S3_ex(X_2[j], Y_1[i], t); AS_31[i][j] = S3_ex(X_3[j], Y_1[i], t);
			AS_41[i][j] = S3_ex(X_4[j], Y_1[i], t); AS_51[i][j] = S3_ex(X_5[j], Y_1[i], t);

			AS_12[i][j] = S3_ex(X_1[j], Y_2[i], t); AS_22[i][j] = S3_ex(X_2[j], Y_2[i], t); AS_32[i][j] = S3_ex(X_3[j], Y_2[i], t);
			AS_42[i][j] = S3_ex(X_4[j], Y_2[i], t); AS_52[i][j] = S3_ex(X_5[j], Y_2[i], t);

			AS_13[i][j] = S3_ex(X_1[j], Y_3[i], t); AS_23[i][j] = S3_ex(X_2[j], Y_3[i], t); AS_33[i][j] = S3_ex(X_3[j], Y_3[i], t);
			AS_43[i][j] = S3_ex(X_4[j], Y_3[i], t); AS_53[i][j] = S3_ex(X_5[j], Y_3[i], t);

			AS_14[i][j] = S3_ex(X_1[j], Y_4[i], t); AS_24[i][j] = S3_ex(X_2[j], Y_4[i], t); AS_34[i][j] = S3_ex(X_3[j], Y_4[i], t);
			AS_44[i][j] = S3_ex(X_4[j], Y_4[i], t); AS_54[i][j] = S3_ex(X_5[j], Y_4[i], t);

			AS_15[i][j] = S3_ex(X_1[j], Y_5[i], t); AS_25[i][j] = S3_ex(X_2[j], Y_5[i], t); AS_35[i][j] = S3_ex(X_3[j], Y_5[i], t);
			AS_45[i][j] = S3_ex(X_4[j], Y_5[i], t); AS_55[i][j] = S3_ex(X_5[j], Y_5[i], t);
		}
	}

	get_Gauss_quadrature_2D(AS_GQ_0, AS_GQ_1, AS_GQ_2, AS_GQ_3, AS_GQ_4, AS_GQ_5, AS_11, AS_21, AS_31, AS_41, AS_51, AS_12, AS_22, AS_32, AS_42, AS_52, AS_13, AS_23, AS_33, AS_43, AS_53, AS_14, AS_24, AS_34, AS_44, AS_54, AS_15, AS_25, AS_35, AS_45, AS_55, h_x, h_y);

	deleteMatrix(N_y + 2, AS_11); deleteMatrix(N_y + 2, AS_12); deleteMatrix(N_y + 2, AS_13); deleteMatrix(N_y + 2, AS_14); deleteMatrix(N_y + 2, AS_15);
	deleteMatrix(N_y + 2, AS_21); deleteMatrix(N_y + 2, AS_22); deleteMatrix(N_y + 2, AS_23); deleteMatrix(N_y + 2, AS_24); deleteMatrix(N_y + 2, AS_25);
	deleteMatrix(N_y + 2, AS_31); deleteMatrix(N_y + 2, AS_32); deleteMatrix(N_y + 2, AS_33); deleteMatrix(N_y + 2, AS_34); deleteMatrix(N_y + 2, AS_35);
	deleteMatrix(N_y + 2, AS_41); deleteMatrix(N_y + 2, AS_42); deleteMatrix(N_y + 2, AS_43); deleteMatrix(N_y + 2, AS_44); deleteMatrix(N_y + 2, AS_45);
	deleteMatrix(N_y + 2, AS_51); deleteMatrix(N_y + 2, AS_52); deleteMatrix(N_y + 2, AS_53); deleteMatrix(N_y + 2, AS_54); deleteMatrix(N_y + 2, AS_55);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j] + Non_0[i][j] + Diff_0[i][j] + Source_0[i][j] + AS_GQ_0[i][j];
			U_1[i][j] = con_1[i][j] + Non_1[i][j] + Diff_1[i][j] + Source_1[i][j] + AS_GQ_1[i][j];
			U_2[i][j] = con_2[i][j] + Non_2[i][j] + Diff_2[i][j] + Source_2[i][j] + AS_GQ_2[i][j];
			U_3[i][j] = con_3[i][j] + Non_3[i][j] + Diff_3[i][j] + Source_3[i][j] + AS_GQ_3[i][j];
			U_4[i][j] = con_4[i][j] + Non_4[i][j] + Diff_4[i][j] + Source_4[i][j] + AS_GQ_4[i][j];
			U_5[i][j] = con_5[i][j] + Non_5[i][j] + Diff_5[i][j] + Source_5[i][j] + AS_GQ_5[i][j];

		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);

	deleteMatrix(N_y + 2, Non_0); deleteMatrix(N_y + 2, Non_1); deleteMatrix(N_y + 2, Non_2);
	deleteMatrix(N_y + 2, Non_3); deleteMatrix(N_y + 2, Non_4); deleteMatrix(N_y + 2, Non_5);

	deleteMatrix(N_y + 2, Diff_0); deleteMatrix(N_y + 2, Diff_1); deleteMatrix(N_y + 2, Diff_2);
	deleteMatrix(N_y + 2, Diff_3); deleteMatrix(N_y + 2, Diff_4); deleteMatrix(N_y + 2, Diff_5);

	deleteMatrix(N_y + 2, Source_0); deleteMatrix(N_y + 2, Source_1); deleteMatrix(N_y + 2, Source_2);
	deleteMatrix(N_y + 2, Source_3); deleteMatrix(N_y + 2, Source_4); deleteMatrix(N_y + 2, Source_5);

	deleteMatrix(N_y + 2, AS_GQ_0); deleteMatrix(N_y + 2, AS_GQ_1); deleteMatrix(N_y + 2, AS_GQ_2);
	deleteMatrix(N_y + 2, AS_GQ_3); deleteMatrix(N_y + 2, AS_GQ_4); deleteMatrix(N_y + 2, AS_GQ_5);

}

void space_i(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t, int index_mode)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	double** Non_0; createMatrix(N_y + 2, N_x + 2, Non_0); double** Non_1; createMatrix(N_y + 2, N_x + 2, Non_1); double** Non_2; createMatrix(N_y + 2, N_x + 2, Non_2);
	double** Non_3; createMatrix(N_y + 2, N_x + 2, Non_3); double** Non_4; createMatrix(N_y + 2, N_x + 2, Non_4); double** Non_5; createMatrix(N_y + 2, N_x + 2, Non_5);

	double** Diff_0; createMatrix(N_y + 2, N_x + 2, Diff_0); double** Diff_1; createMatrix(N_y + 2, N_x + 2, Diff_1); double** Diff_2; createMatrix(N_y + 2, N_x + 2, Diff_2);
	double** Diff_3; createMatrix(N_y + 2, N_x + 2, Diff_3); double** Diff_4; createMatrix(N_y + 2, N_x + 2, Diff_4); double** Diff_5; createMatrix(N_y + 2, N_x + 2, Diff_5);

	double** Source_0; createMatrix(N_y + 2, N_x + 2, Source_0); double** Source_1; createMatrix(N_y + 2, N_x + 2, Source_1); double** Source_2; createMatrix(N_y + 2, N_x + 2, Source_2);
	double** Source_3; createMatrix(N_y + 2, N_x + 2, Source_3); double** Source_4; createMatrix(N_y + 2, N_x + 2, Source_4); double** Source_5; createMatrix(N_y + 2, N_x + 2, Source_5);

	double diff_alpha_i = 0.0;
	if (index_mode == 2)
	{
		diff_alpha_i = Coef_d2 * kappa_i;
	}
	else if (index_mode == 1)
	{
		diff_alpha_i = Coef_d1 * kappa_i;
	}
	else
	{
		diff_alpha_i = kappa_i;
	}

	convection_Ealpha(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, gamma_i, F1_alphac, F2_alphac, h_x, h_y);

	Source_i(Source_0, Source_1, Source_2, Source_3, Source_4, Source_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, h_x, h_y);

	Nonconservetive_Ni(Non_0, Non_1, Non_2, Non_3, Non_4, Non_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	Diffusion_alpha(Diff_0, Diff_1, Diff_2, Diff_3, Diff_4, Diff_5, C_vi, kappa_i, diff_alpha_i, diff_beta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, h_x, h_y);

	// Accuracy source part
	double** AS_11; createMatrix(N_y + 2, N_x + 2, AS_11); double** AS_12; createMatrix(N_y + 2, N_x + 2, AS_12); double** AS_13; createMatrix(N_y + 2, N_x + 2, AS_13); double** AS_14; createMatrix(N_y + 2, N_x + 2, AS_14); double** AS_15; createMatrix(N_y + 2, N_x + 2, AS_15);
	double** AS_21; createMatrix(N_y + 2, N_x + 2, AS_21); double** AS_22; createMatrix(N_y + 2, N_x + 2, AS_22); double** AS_23; createMatrix(N_y + 2, N_x + 2, AS_23); double** AS_24; createMatrix(N_y + 2, N_x + 2, AS_24); double** AS_25; createMatrix(N_y + 2, N_x + 2, AS_25);
	double** AS_31; createMatrix(N_y + 2, N_x + 2, AS_31); double** AS_32; createMatrix(N_y + 2, N_x + 2, AS_32); double** AS_33; createMatrix(N_y + 2, N_x + 2, AS_33); double** AS_34; createMatrix(N_y + 2, N_x + 2, AS_34); double** AS_35; createMatrix(N_y + 2, N_x + 2, AS_35);
	double** AS_41; createMatrix(N_y + 2, N_x + 2, AS_41); double** AS_42; createMatrix(N_y + 2, N_x + 2, AS_42); double** AS_43; createMatrix(N_y + 2, N_x + 2, AS_43); double** AS_44; createMatrix(N_y + 2, N_x + 2, AS_44); double** AS_45; createMatrix(N_y + 2, N_x + 2, AS_45);
	double** AS_51; createMatrix(N_y + 2, N_x + 2, AS_51); double** AS_52; createMatrix(N_y + 2, N_x + 2, AS_52); double** AS_53; createMatrix(N_y + 2, N_x + 2, AS_53); double** AS_54; createMatrix(N_y + 2, N_x + 2, AS_54); double** AS_55; createMatrix(N_y + 2, N_x + 2, AS_55);

	double** AS_GQ_0; createMatrix(N_y + 2, N_x + 2, AS_GQ_0); double** AS_GQ_1; createMatrix(N_y + 2, N_x + 2, AS_GQ_1);
	double** AS_GQ_2; createMatrix(N_y + 2, N_x + 2, AS_GQ_2); double** AS_GQ_3; createMatrix(N_y + 2, N_x + 2, AS_GQ_3);
	double** AS_GQ_4; createMatrix(N_y + 2, N_x + 2, AS_GQ_4); double** AS_GQ_5; createMatrix(N_y + 2, N_x + 2, AS_GQ_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			AS_11[i][j] = S4_ex(X_1[j], Y_1[i], t); AS_21[i][j] = S4_ex(X_2[j], Y_1[i], t); AS_31[i][j] = S4_ex(X_3[j], Y_1[i], t);
			AS_41[i][j] = S4_ex(X_4[j], Y_1[i], t); AS_51[i][j] = S4_ex(X_5[j], Y_1[i], t);

			AS_12[i][j] = S4_ex(X_1[j], Y_2[i], t); AS_22[i][j] = S4_ex(X_2[j], Y_2[i], t); AS_32[i][j] = S4_ex(X_3[j], Y_2[i], t);
			AS_42[i][j] = S4_ex(X_4[j], Y_2[i], t); AS_52[i][j] = S4_ex(X_5[j], Y_2[i], t);

			AS_13[i][j] = S4_ex(X_1[j], Y_3[i], t); AS_23[i][j] = S4_ex(X_2[j], Y_3[i], t); AS_33[i][j] = S4_ex(X_3[j], Y_3[i], t);
			AS_43[i][j] = S4_ex(X_4[j], Y_3[i], t); AS_53[i][j] = S4_ex(X_5[j], Y_3[i], t);

			AS_14[i][j] = S4_ex(X_1[j], Y_4[i], t); AS_24[i][j] = S4_ex(X_2[j], Y_4[i], t); AS_34[i][j] = S4_ex(X_3[j], Y_4[i], t);
			AS_44[i][j] = S4_ex(X_4[j], Y_4[i], t); AS_54[i][j] = S4_ex(X_5[j], Y_4[i], t);

			AS_15[i][j] = S4_ex(X_1[j], Y_5[i], t); AS_25[i][j] = S4_ex(X_2[j], Y_5[i], t); AS_35[i][j] = S4_ex(X_3[j], Y_5[i], t);
			AS_45[i][j] = S4_ex(X_4[j], Y_5[i], t); AS_55[i][j] = S4_ex(X_5[j], Y_5[i], t);
		}
	}

	get_Gauss_quadrature_2D(AS_GQ_0, AS_GQ_1, AS_GQ_2, AS_GQ_3, AS_GQ_4, AS_GQ_5, AS_11, AS_21, AS_31, AS_41, AS_51, AS_12, AS_22, AS_32, AS_42, AS_52, AS_13, AS_23, AS_33, AS_43, AS_53, AS_14, AS_24, AS_34, AS_44, AS_54, AS_15, AS_25, AS_35, AS_45, AS_55, h_x, h_y);

	deleteMatrix(N_y + 2, AS_11); deleteMatrix(N_y + 2, AS_12); deleteMatrix(N_y + 2, AS_13); deleteMatrix(N_y + 2, AS_14); deleteMatrix(N_y + 2, AS_15);
	deleteMatrix(N_y + 2, AS_21); deleteMatrix(N_y + 2, AS_22); deleteMatrix(N_y + 2, AS_23); deleteMatrix(N_y + 2, AS_24); deleteMatrix(N_y + 2, AS_25);
	deleteMatrix(N_y + 2, AS_31); deleteMatrix(N_y + 2, AS_32); deleteMatrix(N_y + 2, AS_33); deleteMatrix(N_y + 2, AS_34); deleteMatrix(N_y + 2, AS_35);
	deleteMatrix(N_y + 2, AS_41); deleteMatrix(N_y + 2, AS_42); deleteMatrix(N_y + 2, AS_43); deleteMatrix(N_y + 2, AS_44); deleteMatrix(N_y + 2, AS_45);
	deleteMatrix(N_y + 2, AS_51); deleteMatrix(N_y + 2, AS_52); deleteMatrix(N_y + 2, AS_53); deleteMatrix(N_y + 2, AS_54); deleteMatrix(N_y + 2, AS_55);


	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j] + Non_0[i][j] + Diff_0[i][j] + Source_0[i][j] + AS_GQ_0[i][j];
			U_1[i][j] = con_1[i][j] + Non_1[i][j] + Diff_1[i][j] + Source_1[i][j] + AS_GQ_1[i][j];
			U_2[i][j] = con_2[i][j] + Non_2[i][j] + Diff_2[i][j] + Source_2[i][j] + AS_GQ_2[i][j];
			U_3[i][j] = con_3[i][j] + Non_3[i][j] + Diff_3[i][j] + Source_3[i][j] + AS_GQ_3[i][j];
			U_4[i][j] = con_4[i][j] + Non_4[i][j] + Diff_4[i][j] + Source_4[i][j] + AS_GQ_4[i][j];
			U_5[i][j] = con_5[i][j] + Non_5[i][j] + Diff_5[i][j] + Source_5[i][j] + AS_GQ_5[i][j];

		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);

	deleteMatrix(N_y + 2, Non_0); deleteMatrix(N_y + 2, Non_1); deleteMatrix(N_y + 2, Non_2);
	deleteMatrix(N_y + 2, Non_3); deleteMatrix(N_y + 2, Non_4); deleteMatrix(N_y + 2, Non_5);

	deleteMatrix(N_y + 2, Diff_0); deleteMatrix(N_y + 2, Diff_1); deleteMatrix(N_y + 2, Diff_2);
	deleteMatrix(N_y + 2, Diff_3); deleteMatrix(N_y + 2, Diff_4); deleteMatrix(N_y + 2, Diff_5);

	deleteMatrix(N_y + 2, Source_0); deleteMatrix(N_y + 2, Source_1); deleteMatrix(N_y + 2, Source_2);
	deleteMatrix(N_y + 2, Source_3); deleteMatrix(N_y + 2, Source_4); deleteMatrix(N_y + 2, Source_5);

	deleteMatrix(N_y + 2, AS_GQ_0); deleteMatrix(N_y + 2, AS_GQ_1); deleteMatrix(N_y + 2, AS_GQ_2);
	deleteMatrix(N_y + 2, AS_GQ_3); deleteMatrix(N_y + 2, AS_GQ_4); deleteMatrix(N_y + 2, AS_GQ_5);
}

void space_r(double** U_0, double** U_1, double** U_2, double** U_3, double** U_4, double** U_5, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* Y_1, double* Y_2, double* Y_3, double* Y_4, double* Y_5, double** coef_L1, double** coef_R1, double** coef_L2, double** coef_R2, double** coef_L3, double** coef_R3, double** coef_L4, double** coef_R4, double** coef_L5, double** coef_R5, double** coef_B1, double** coef_T1, double** coef_B2, double** coef_T2, double** coef_B3, double** coef_T3, double** coef_B4, double** coef_T4, double** coef_B5, double** coef_T5, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double** F1_alphac, double** F2_alphac, double h_x, double h_y, double t, int index_mode)
{
	double** con_0; createMatrix(N_y + 2, N_x + 2, con_0); double** con_1; createMatrix(N_y + 2, N_x + 2, con_1); double** con_2; createMatrix(N_y + 2, N_x + 2, con_2);
	double** con_3; createMatrix(N_y + 2, N_x + 2, con_3); double** con_4; createMatrix(N_y + 2, N_x + 2, con_4); double** con_5; createMatrix(N_y + 2, N_x + 2, con_5);

	double** Non_0; createMatrix(N_y + 2, N_x + 2, Non_0); double** Non_1; createMatrix(N_y + 2, N_x + 2, Non_1); double** Non_2; createMatrix(N_y + 2, N_x + 2, Non_2);
	double** Non_3; createMatrix(N_y + 2, N_x + 2, Non_3); double** Non_4; createMatrix(N_y + 2, N_x + 2, Non_4); double** Non_5; createMatrix(N_y + 2, N_x + 2, Non_5);

	double** Diff_0; createMatrix(N_y + 2, N_x + 2, Diff_0); double** Diff_1; createMatrix(N_y + 2, N_x + 2, Diff_1); double** Diff_2; createMatrix(N_y + 2, N_x + 2, Diff_2);
	double** Diff_3; createMatrix(N_y + 2, N_x + 2, Diff_3); double** Diff_4; createMatrix(N_y + 2, N_x + 2, Diff_4); double** Diff_5; createMatrix(N_y + 2, N_x + 2, Diff_5);

	double** Source_0; createMatrix(N_y + 2, N_x + 2, Source_0); double** Source_1; createMatrix(N_y + 2, N_x + 2, Source_1); double** Source_2; createMatrix(N_y + 2, N_x + 2, Source_2);
	double** Source_3; createMatrix(N_y + 2, N_x + 2, Source_3); double** Source_4; createMatrix(N_y + 2, N_x + 2, Source_4); double** Source_5; createMatrix(N_y + 2, N_x + 2, Source_5);

	double diff_alpha_r = 0.0;
	if (index_mode == 2)
	{
		diff_alpha_r = Coef_d2 * kappa_r;
	}
	else if (index_mode == 1)
	{
		diff_alpha_r = Coef_d1 * kappa_r;
	}
	else
	{
		diff_alpha_r = kappa_r;
	}

	convection_Ealpha(con_0, con_1, con_2, con_3, con_4, con_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, gamma_r, F1_alphac, F2_alphac, h_x, h_y);

	Source_r(Source_0, Source_1, Source_2, Source_3, Source_4, Source_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	Nonconservetive_Nr(Non_0, Non_1, Non_2, Non_3, Non_4, Non_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	Diffusion_r(Diff_0, Diff_1, Diff_2, Diff_3, Diff_4, Diff_5, kappa_r, diff_alpha_r, diff_beta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);


	// Accuracy source part
	double** AS_11; createMatrix(N_y + 2, N_x + 2, AS_11); double** AS_12; createMatrix(N_y + 2, N_x + 2, AS_12); double** AS_13; createMatrix(N_y + 2, N_x + 2, AS_13); double** AS_14; createMatrix(N_y + 2, N_x + 2, AS_14); double** AS_15; createMatrix(N_y + 2, N_x + 2, AS_15);
	double** AS_21; createMatrix(N_y + 2, N_x + 2, AS_21); double** AS_22; createMatrix(N_y + 2, N_x + 2, AS_22); double** AS_23; createMatrix(N_y + 2, N_x + 2, AS_23); double** AS_24; createMatrix(N_y + 2, N_x + 2, AS_24); double** AS_25; createMatrix(N_y + 2, N_x + 2, AS_25);
	double** AS_31; createMatrix(N_y + 2, N_x + 2, AS_31); double** AS_32; createMatrix(N_y + 2, N_x + 2, AS_32); double** AS_33; createMatrix(N_y + 2, N_x + 2, AS_33); double** AS_34; createMatrix(N_y + 2, N_x + 2, AS_34); double** AS_35; createMatrix(N_y + 2, N_x + 2, AS_35);
	double** AS_41; createMatrix(N_y + 2, N_x + 2, AS_41); double** AS_42; createMatrix(N_y + 2, N_x + 2, AS_42); double** AS_43; createMatrix(N_y + 2, N_x + 2, AS_43); double** AS_44; createMatrix(N_y + 2, N_x + 2, AS_44); double** AS_45; createMatrix(N_y + 2, N_x + 2, AS_45);
	double** AS_51; createMatrix(N_y + 2, N_x + 2, AS_51); double** AS_52; createMatrix(N_y + 2, N_x + 2, AS_52); double** AS_53; createMatrix(N_y + 2, N_x + 2, AS_53); double** AS_54; createMatrix(N_y + 2, N_x + 2, AS_54); double** AS_55; createMatrix(N_y + 2, N_x + 2, AS_55);

	double** AS_GQ_0; createMatrix(N_y + 2, N_x + 2, AS_GQ_0); double** AS_GQ_1; createMatrix(N_y + 2, N_x + 2, AS_GQ_1);
	double** AS_GQ_2; createMatrix(N_y + 2, N_x + 2, AS_GQ_2); double** AS_GQ_3; createMatrix(N_y + 2, N_x + 2, AS_GQ_3);
	double** AS_GQ_4; createMatrix(N_y + 2, N_x + 2, AS_GQ_4); double** AS_GQ_5; createMatrix(N_y + 2, N_x + 2, AS_GQ_5);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			AS_11[i][j] = S5_ex(X_1[j], Y_1[i], t); AS_21[i][j] = S5_ex(X_2[j], Y_1[i], t); AS_31[i][j] = S5_ex(X_3[j], Y_1[i], t);
			AS_41[i][j] = S5_ex(X_4[j], Y_1[i], t); AS_51[i][j] = S5_ex(X_5[j], Y_1[i], t);

			AS_12[i][j] = S5_ex(X_1[j], Y_2[i], t); AS_22[i][j] = S5_ex(X_2[j], Y_2[i], t); AS_32[i][j] = S5_ex(X_3[j], Y_2[i], t);
			AS_42[i][j] = S5_ex(X_4[j], Y_2[i], t); AS_52[i][j] = S5_ex(X_5[j], Y_2[i], t);

			AS_13[i][j] = S5_ex(X_1[j], Y_3[i], t); AS_23[i][j] = S5_ex(X_2[j], Y_3[i], t); AS_33[i][j] = S5_ex(X_3[j], Y_3[i], t);
			AS_43[i][j] = S5_ex(X_4[j], Y_3[i], t); AS_53[i][j] = S5_ex(X_5[j], Y_3[i], t);

			AS_14[i][j] = S5_ex(X_1[j], Y_4[i], t); AS_24[i][j] = S5_ex(X_2[j], Y_4[i], t); AS_34[i][j] = S5_ex(X_3[j], Y_4[i], t);
			AS_44[i][j] = S5_ex(X_4[j], Y_4[i], t); AS_54[i][j] = S5_ex(X_5[j], Y_4[i], t);

			AS_15[i][j] = S5_ex(X_1[j], Y_5[i], t); AS_25[i][j] = S5_ex(X_2[j], Y_5[i], t); AS_35[i][j] = S5_ex(X_3[j], Y_5[i], t);
			AS_45[i][j] = S5_ex(X_4[j], Y_5[i], t); AS_55[i][j] = S5_ex(X_5[j], Y_5[i], t);
		}
	}

	get_Gauss_quadrature_2D(AS_GQ_0, AS_GQ_1, AS_GQ_2, AS_GQ_3, AS_GQ_4, AS_GQ_5, AS_11, AS_21, AS_31, AS_41, AS_51, AS_12, AS_22, AS_32, AS_42, AS_52, AS_13, AS_23, AS_33, AS_43, AS_53, AS_14, AS_24, AS_34, AS_44, AS_54, AS_15, AS_25, AS_35, AS_45, AS_55, h_x, h_y);

	deleteMatrix(N_y + 2, AS_11); deleteMatrix(N_y + 2, AS_12); deleteMatrix(N_y + 2, AS_13); deleteMatrix(N_y + 2, AS_14); deleteMatrix(N_y + 2, AS_15);
	deleteMatrix(N_y + 2, AS_21); deleteMatrix(N_y + 2, AS_22); deleteMatrix(N_y + 2, AS_23); deleteMatrix(N_y + 2, AS_24); deleteMatrix(N_y + 2, AS_25);
	deleteMatrix(N_y + 2, AS_31); deleteMatrix(N_y + 2, AS_32); deleteMatrix(N_y + 2, AS_33); deleteMatrix(N_y + 2, AS_34); deleteMatrix(N_y + 2, AS_35);
	deleteMatrix(N_y + 2, AS_41); deleteMatrix(N_y + 2, AS_42); deleteMatrix(N_y + 2, AS_43); deleteMatrix(N_y + 2, AS_44); deleteMatrix(N_y + 2, AS_45);
	deleteMatrix(N_y + 2, AS_51); deleteMatrix(N_y + 2, AS_52); deleteMatrix(N_y + 2, AS_53); deleteMatrix(N_y + 2, AS_54); deleteMatrix(N_y + 2, AS_55);

	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			U_0[i][j] = con_0[i][j] + Non_0[i][j] + Diff_0[i][j] + Source_0[i][j] + AS_GQ_0[i][j];
			U_1[i][j] = con_1[i][j] + Non_1[i][j] + Diff_1[i][j] + Source_1[i][j] + AS_GQ_1[i][j];
			U_2[i][j] = con_2[i][j] + Non_2[i][j] + Diff_2[i][j] + Source_2[i][j] + AS_GQ_2[i][j];
			U_3[i][j] = con_3[i][j] + Non_3[i][j] + Diff_3[i][j] + Source_3[i][j] + AS_GQ_3[i][j];
			U_4[i][j] = con_4[i][j] + Non_4[i][j] + Diff_4[i][j] + Source_4[i][j] + AS_GQ_4[i][j];
			U_5[i][j] = con_5[i][j] + Non_5[i][j] + Diff_5[i][j] + Source_5[i][j] + AS_GQ_5[i][j];
		}
	}

	deleteMatrix(N_y + 2, con_0); deleteMatrix(N_y + 2, con_1); deleteMatrix(N_y + 2, con_2);
	deleteMatrix(N_y + 2, con_3); deleteMatrix(N_y + 2, con_4); deleteMatrix(N_y + 2, con_5);

	deleteMatrix(N_y + 2, Non_0); deleteMatrix(N_y + 2, Non_1); deleteMatrix(N_y + 2, Non_2);
	deleteMatrix(N_y + 2, Non_3); deleteMatrix(N_y + 2, Non_4); deleteMatrix(N_y + 2, Non_5);

	deleteMatrix(N_y + 2, Diff_0); deleteMatrix(N_y + 2, Diff_1); deleteMatrix(N_y + 2, Diff_2);
	deleteMatrix(N_y + 2, Diff_3); deleteMatrix(N_y + 2, Diff_4); deleteMatrix(N_y + 2, Diff_5);

	deleteMatrix(N_y + 2, Source_0); deleteMatrix(N_y + 2, Source_1); deleteMatrix(N_y + 2, Source_2);
	deleteMatrix(N_y + 2, Source_3); deleteMatrix(N_y + 2, Source_4); deleteMatrix(N_y + 2, Source_5);

	deleteMatrix(N_y + 2, AS_GQ_0); deleteMatrix(N_y + 2, AS_GQ_1); deleteMatrix(N_y + 2, AS_GQ_2);
	deleteMatrix(N_y + 2, AS_GQ_3); deleteMatrix(N_y + 2, AS_GQ_4); deleteMatrix(N_y + 2, AS_GQ_5);
}
