#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG2D.h"
#include "space.h"
#include "OE.h"


void OE_compute_u(double** u, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m_0, double** m_1, double** m_2, double** m_3, double** m_4, double** m_5)
{
	double** m; createMatrix(N_y + 2, N_x + 2, m); double** rho; createMatrix(N_y + 2, N_x + 2, rho);
	F_d(m, xi, eta, m_0, m_1, m_2, m_3, m_4, m_5); F_d(rho, xi, eta, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5);

	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			u[i - 1][j - 1] = m[i][j] / rho[i][j];
		}
	}
	deleteMatrix(N_y + 2, m); deleteMatrix(N_y + 2, rho);
}

void OE_compute_u_ave(double** u_ave, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m_0, double** m_1, double** m_2, double** m_3, double** m_4, double** m_5)
{
	double** u_11; createMatrix(N_y, N_x, u_11); double** u_12; createMatrix(N_y, N_x, u_12); double** u_13; createMatrix(N_y, N_x, u_13);
	double** u_14; createMatrix(N_y, N_x, u_14); double** u_15; createMatrix(N_y, N_x, u_15);

	double** u_21; createMatrix(N_y, N_x, u_21); double** u_22; createMatrix(N_y, N_x, u_22); double** u_23; createMatrix(N_y, N_x, u_23);
	double** u_24; createMatrix(N_y, N_x, u_24); double** u_25; createMatrix(N_y, N_x, u_25);

	double** u_31; createMatrix(N_y, N_x, u_31); double** u_32; createMatrix(N_y, N_x, u_32); double** u_33; createMatrix(N_y, N_x, u_33);
	double** u_34; createMatrix(N_y, N_x, u_34); double** u_35; createMatrix(N_y, N_x, u_35);

	double** u_41; createMatrix(N_y, N_x, u_41); double** u_42; createMatrix(N_y, N_x, u_42); double** u_43; createMatrix(N_y, N_x, u_43);
	double** u_44; createMatrix(N_y, N_x, u_44); double** u_45; createMatrix(N_y, N_x, u_45);

	double** u_51; createMatrix(N_y, N_x, u_51); double** u_52; createMatrix(N_y, N_x, u_52); double** u_53; createMatrix(N_y, N_x, u_53);
	double** u_54; createMatrix(N_y, N_x, u_54); double** u_55; createMatrix(N_y, N_x, u_55);

	OE_compute_u(u_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);

	OE_compute_u(u_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);

	OE_compute_u(u_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);

	OE_compute_u(u_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);

	OE_compute_u(u_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);
	OE_compute_u(u_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m_0, m_1, m_2, m_3, m_4, m_5);

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

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			u_ave[i][j] = 0.25 * (w_11 * u_11[i][j] + w_12 * u_12[i][j] + w_13 * u_13[i][j] + w_14 * u_14[i][j] + w_15 * u_15[i][j] + w_21 * u_21[i][j] + w_22 * u_22[i][j] + w_23 * u_23[i][j] + w_24 * u_24[i][j] + w_25 * u_25[i][j] + w_31 * u_31[i][j] + w_32 * u_32[i][j] + w_33 * u_33[i][j] + w_34 * u_34[i][j] + w_35 * u_35[i][j] + w_41 * u_41[i][j] + w_42 * u_42[i][j] + w_43 * u_43[i][j] + w_44 * u_44[i][j] + w_45 * u_45[i][j] + w_51 * u_51[i][j] + w_52 * u_52[i][j] + w_53 * u_53[i][j] + w_54 * u_54[i][j] + w_55 * u_55[i][j]);
		}
	}

	deleteMatrix(N_y, u_11); deleteMatrix(N_y, u_12); deleteMatrix(N_y, u_13); deleteMatrix(N_y, u_14); deleteMatrix(N_y, u_15);
	deleteMatrix(N_y, u_21); deleteMatrix(N_y, u_22); deleteMatrix(N_y, u_23); deleteMatrix(N_y, u_24); deleteMatrix(N_y, u_25);
	deleteMatrix(N_y, u_31); deleteMatrix(N_y, u_32); deleteMatrix(N_y, u_33); deleteMatrix(N_y, u_34); deleteMatrix(N_y, u_35);
	deleteMatrix(N_y, u_41); deleteMatrix(N_y, u_42); deleteMatrix(N_y, u_43); deleteMatrix(N_y, u_44); deleteMatrix(N_y, u_45);
	deleteMatrix(N_y, u_51); deleteMatrix(N_y, u_52); deleteMatrix(N_y, u_53); deleteMatrix(N_y, u_54); deleteMatrix(N_y, u_55);
}

void OE_compute_e_alpha(double** e_alpha, double xi, double eta, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			double temp_rho = rho_0[i][j] + rho_1[i][j] * xi + rho_2[i][j] * eta + rho_3[i][j] * xi * eta + rho_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + rho_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			double temp_rhou2 = pow((m1_0[i][j] + m1_1[i][j] * xi + m1_2[i][j] * eta + m1_3[i][j] * xi * eta + m1_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m1_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_rhov2 = pow((m2_0[i][j] + m2_1[i][j] * xi + m2_2[i][j] * eta + m2_3[i][j] * xi * eta + m2_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + m2_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0), 2) / temp_rho;
			double temp_Ealpha = Ealpha_0[i][j] + Ealpha_1[i][j] * xi + Ealpha_2[i][j] * eta + Ealpha_3[i][j] * xi * eta + Ealpha_4[i][j] * (3.0 * pow(xi, 2) - 1.0) / 2.0 + Ealpha_5[i][j] * (3.0 * pow(eta, 2) - 1.0) / 2.0;
			e_alpha[i - 1][j - 1] = (temp_Ealpha - temp_rhou2 / 6.0 - temp_rhov2 / 6.0) / temp_rho;

		}
	}
}

void OE_compute_e_alpha_ave(double** u_ave, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ealpha_0, double** Ealpha_1, double** Ealpha_2, double** Ealpha_3, double** Ealpha_4, double** Ealpha_5)
{
	double** u_11; createMatrix(N_y, N_x, u_11); double** u_12; createMatrix(N_y, N_x, u_12); double** u_13; createMatrix(N_y, N_x, u_13);
	double** u_14; createMatrix(N_y, N_x, u_14); double** u_15; createMatrix(N_y, N_x, u_15);

	double** u_21; createMatrix(N_y, N_x, u_21); double** u_22; createMatrix(N_y, N_x, u_22); double** u_23; createMatrix(N_y, N_x, u_23);
	double** u_24; createMatrix(N_y, N_x, u_24); double** u_25; createMatrix(N_y, N_x, u_25);

	double** u_31; createMatrix(N_y, N_x, u_31); double** u_32; createMatrix(N_y, N_x, u_32); double** u_33; createMatrix(N_y, N_x, u_33);
	double** u_34; createMatrix(N_y, N_x, u_34); double** u_35; createMatrix(N_y, N_x, u_35);

	double** u_41; createMatrix(N_y, N_x, u_41); double** u_42; createMatrix(N_y, N_x, u_42); double** u_43; createMatrix(N_y, N_x, u_43);
	double** u_44; createMatrix(N_y, N_x, u_44); double** u_45; createMatrix(N_y, N_x, u_45);

	double** u_51; createMatrix(N_y, N_x, u_51); double** u_52; createMatrix(N_y, N_x, u_52); double** u_53; createMatrix(N_y, N_x, u_53);
	double** u_54; createMatrix(N_y, N_x, u_54); double** u_55; createMatrix(N_y, N_x, u_55);

	OE_compute_e_alpha(u_11, GQxi_1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_12, GQxi_1, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_13, GQxi_1, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_14, GQxi_1, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_15, GQxi_1, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	OE_compute_e_alpha(u_21, GQxi_2, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_22, GQxi_2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_23, GQxi_2, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_24, GQxi_2, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_25, GQxi_2, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	OE_compute_e_alpha(u_31, GQxi_3, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_32, GQxi_3, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_33, GQxi_3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_34, GQxi_3, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_35, GQxi_3, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	OE_compute_e_alpha(u_41, GQxi_4, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_42, GQxi_4, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_43, GQxi_4, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_44, GQxi_4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_45, GQxi_4, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

	OE_compute_e_alpha(u_51, GQxi_5, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_52, GQxi_5, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_53, GQxi_5, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_54, GQxi_5, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);
	OE_compute_e_alpha(u_55, GQxi_5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ealpha_0, Ealpha_1, Ealpha_2, Ealpha_3, Ealpha_4, Ealpha_5);

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

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			u_ave[i][j] = 0.25 * (w_11 * u_11[i][j] + w_12 * u_12[i][j] + w_13 * u_13[i][j] + w_14 * u_14[i][j] + w_15 * u_15[i][j] + w_21 * u_21[i][j] + w_22 * u_22[i][j] + w_23 * u_23[i][j] + w_24 * u_24[i][j] + w_25 * u_25[i][j] + w_31 * u_31[i][j] + w_32 * u_32[i][j] + w_33 * u_33[i][j] + w_34 * u_34[i][j] + w_35 * u_35[i][j] + w_41 * u_41[i][j] + w_42 * u_42[i][j] + w_43 * u_43[i][j] + w_44 * u_44[i][j] + w_45 * u_45[i][j] + w_51 * u_51[i][j] + w_52 * u_52[i][j] + w_53 * u_53[i][j] + w_54 * u_54[i][j] + w_55 * u_55[i][j]);
		}
	}

	deleteMatrix(N_y, u_11); deleteMatrix(N_y, u_12); deleteMatrix(N_y, u_13); deleteMatrix(N_y, u_14); deleteMatrix(N_y, u_15);
	deleteMatrix(N_y, u_21); deleteMatrix(N_y, u_22); deleteMatrix(N_y, u_23); deleteMatrix(N_y, u_24); deleteMatrix(N_y, u_25);
	deleteMatrix(N_y, u_31); deleteMatrix(N_y, u_32); deleteMatrix(N_y, u_33); deleteMatrix(N_y, u_34); deleteMatrix(N_y, u_35);
	deleteMatrix(N_y, u_41); deleteMatrix(N_y, u_42); deleteMatrix(N_y, u_43); deleteMatrix(N_y, u_44); deleteMatrix(N_y, u_45);
	deleteMatrix(N_y, u_51); deleteMatrix(N_y, u_52); deleteMatrix(N_y, u_53); deleteMatrix(N_y, u_54); deleteMatrix(N_y, u_55);
}

void OE_compute_Cs_ave(double** Cs_ave, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** e_e_ave; createMatrix(N_y, N_x, e_e_ave); double** e_i_ave; createMatrix(N_y, N_x, e_i_ave);
	double** e_r_ave; createMatrix(N_y, N_x, e_r_ave);

	OE_compute_e_alpha_ave(e_e_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5);
	OE_compute_e_alpha_ave(e_i_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5);
	OE_compute_e_alpha_ave(e_r_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			Cs_ave[i][j] = sqrt(gamma_e * (gamma_e - 1.0) * e_e_ave[i][j] + gamma_i * (gamma_i - 1.0) * e_i_ave[i][j] + gamma_r * (gamma_r - 1.0) * e_r_ave[i][j]);

		}
	}
	deleteMatrix(N_y, e_e_ave); deleteMatrix(N_y, e_i_ave); deleteMatrix(N_y, e_r_ave);
}

void OE_compute_beta_ij(double** beta_ij_x, double** beta_ij_y, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5)
{
	double** u_ave; createMatrix(N_y, N_x, u_ave); double** v_ave; createMatrix(N_y, N_x, v_ave);
	double** Cs_ave; createMatrix(N_y, N_x, Cs_ave);

	OE_compute_u_ave(u_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5);
	OE_compute_u_ave(v_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5);
	OE_compute_Cs_ave(Cs_ave, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_max = abs(u_ave[i][j] - Cs_ave[i][j]);
			if (abs(u_ave[i][j]) > temp_max)
			{
				temp_max = abs(u_ave[i][j]);
			}
			if (abs(u_ave[i][j] + Cs_ave[i][j]) > temp_max)
			{
				temp_max = abs(u_ave[i][j] + Cs_ave[i][j]);
			}
			beta_ij_x[i][j] = temp_max;
		}
	}

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double temp_max = abs(v_ave[i][j] - Cs_ave[i][j]);
			if (abs(v_ave[i][j]) > temp_max)
			{
				temp_max = abs(v_ave[i][j]);
			}
			if (abs(v_ave[i][j] + Cs_ave[i][j]) > temp_max)
			{
				temp_max = abs(v_ave[i][j] + Cs_ave[i][j]);
			}
			beta_ij_y[i][j] = temp_max;
		}
	}

	deleteMatrix(N_y, u_ave); deleteMatrix(N_y, v_ave); deleteMatrix(N_y, Cs_ave);
}

void OE_compute_scalar_sigma(double* Linf, double** sigma_0_x, double** sigma_1_x, double** sigma_2_x, double** sigma_0_y, double** sigma_1_y, double** sigma_2_y, double** u_0, double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double h_x, double h_y)
{
	double u_Gave = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			u_Gave = u_Gave + u_0[i][j];
		}
	}
	u_Gave = u_Gave / (N_x * N_y);

	double** u_L1; createMatrix(N_y + 2, N_x + 2, u_L1); double** u_L2; createMatrix(N_y + 2, N_x + 2, u_L2);
	double** u_L3; createMatrix(N_y + 2, N_x + 2, u_L3); double** u_L4; createMatrix(N_y + 2, N_x + 2, u_L4);
	double** u_L5; createMatrix(N_y + 2, N_x + 2, u_L5);

	double** u_R1; createMatrix(N_y + 2, N_x + 2, u_R1); double** u_R2; createMatrix(N_y + 2, N_x + 2, u_R2);
	double** u_R3; createMatrix(N_y + 2, N_x + 2, u_R3); double** u_R4; createMatrix(N_y + 2, N_x + 2, u_R4);
	double** u_R5; createMatrix(N_y + 2, N_x + 2, u_R5);

	double** u_T1; createMatrix(N_y + 2, N_x + 2, u_T1); double** u_T2; createMatrix(N_y + 2, N_x + 2, u_T2);
	double** u_T3; createMatrix(N_y + 2, N_x + 2, u_T3); double** u_T4; createMatrix(N_y + 2, N_x + 2, u_T4);
	double** u_T5; createMatrix(N_y + 2, N_x + 2, u_T5);

	double** u_B1; createMatrix(N_y + 2, N_x + 2, u_B1); double** u_B2; createMatrix(N_y + 2, N_x + 2, u_B2);
	double** u_B3; createMatrix(N_y + 2, N_x + 2, u_B3); double** u_B4; createMatrix(N_y + 2, N_x + 2, u_B4);
	double** u_B5; createMatrix(N_y + 2, N_x + 2, u_B5);

	double** u_GQ11; createMatrix(N_y + 2, N_x + 2, u_GQ11); double** u_GQ12; createMatrix(N_y + 2, N_x + 2, u_GQ12);
	double** u_GQ13; createMatrix(N_y + 2, N_x + 2, u_GQ13); double** u_GQ14; createMatrix(N_y + 2, N_x + 2, u_GQ14);
	double** u_GQ15; createMatrix(N_y + 2, N_x + 2, u_GQ15);

	double** u_GQ21; createMatrix(N_y + 2, N_x + 2, u_GQ21); double** u_GQ22; createMatrix(N_y + 2, N_x + 2, u_GQ22);
	double** u_GQ23; createMatrix(N_y + 2, N_x + 2, u_GQ23); double** u_GQ24; createMatrix(N_y + 2, N_x + 2, u_GQ24);
	double** u_GQ25; createMatrix(N_y + 2, N_x + 2, u_GQ25);

	double** u_GQ31; createMatrix(N_y + 2, N_x + 2, u_GQ31); double** u_GQ32; createMatrix(N_y + 2, N_x + 2, u_GQ32);
	double** u_GQ33; createMatrix(N_y + 2, N_x + 2, u_GQ33); double** u_GQ34; createMatrix(N_y + 2, N_x + 2, u_GQ34);
	double** u_GQ35; createMatrix(N_y + 2, N_x + 2, u_GQ35);

	double** u_GQ41; createMatrix(N_y + 2, N_x + 2, u_GQ41); double** u_GQ42; createMatrix(N_y + 2, N_x + 2, u_GQ42);
	double** u_GQ43; createMatrix(N_y + 2, N_x + 2, u_GQ43); double** u_GQ44; createMatrix(N_y + 2, N_x + 2, u_GQ44);
	double** u_GQ45; createMatrix(N_y + 2, N_x + 2, u_GQ45);

	double** u_GQ51; createMatrix(N_y + 2, N_x + 2, u_GQ51); double** u_GQ52; createMatrix(N_y + 2, N_x + 2, u_GQ52);
	double** u_GQ53; createMatrix(N_y + 2, N_x + 2, u_GQ53); double** u_GQ54; createMatrix(N_y + 2, N_x + 2, u_GQ54);
	double** u_GQ55; createMatrix(N_y + 2, N_x + 2, u_GQ55);

	F_d(u_L1, -1.0, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_L2, -1.0, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_L3, -1.0, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_L4, -1.0, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_L5, -1.0, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_R1, 1.0, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_R2, 1.0, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_R3, 1.0, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_R4, 1.0, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_R5, 1.0, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_B1, GQxi_1, -1.0, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_B2, GQxi_2, -1.0, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_B3, GQxi_3, -1.0, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_B4, GQxi_4, -1.0, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_B5, GQxi_5, -1.0, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_T1, GQxi_1, 1.0, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_T2, GQxi_2, 1.0, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_T3, GQxi_3, 1.0, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_T4, GQxi_4, 1.0, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_T5, GQxi_5, 1.0, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_GQ11, GQxi_1, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ12, GQxi_1, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ13, GQxi_1, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ14, GQxi_1, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ15, GQxi_1, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_GQ21, GQxi_2, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ22, GQxi_2, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ23, GQxi_2, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ24, GQxi_2, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ25, GQxi_2, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_GQ31, GQxi_3, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ32, GQxi_3, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ33, GQxi_3, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ34, GQxi_3, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ35, GQxi_3, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_GQ41, GQxi_4, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ42, GQxi_4, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ43, GQxi_4, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ44, GQxi_4, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ45, GQxi_4, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	F_d(u_GQ51, GQxi_5, GQxi_1, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ52, GQxi_5, GQxi_2, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ53, GQxi_5, GQxi_3, u_0, u_1, u_2, u_3, u_4, u_5); F_d(u_GQ54, GQxi_5, GQxi_4, u_0, u_1, u_2, u_3, u_4, u_5);
	F_d(u_GQ55, GQxi_5, GQxi_5, u_0, u_1, u_2, u_3, u_4, u_5);

	double u_GLinf = 0.0;
	for (int i = 1; i < N_y + 1; i++)
	{
		for (int j = 1; j < N_x + 1; j++)
		{
			if (abs(u_L1[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_L1[i][j] - u_Gave);
			}
			if (abs(u_L2[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_L2[i][j] - u_Gave);
			}
			if (abs(u_L3[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_L3[i][j] - u_Gave);
			}
			if (abs(u_L4[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_L4[i][j] - u_Gave);
			}
			if (abs(u_L5[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_L5[i][j] - u_Gave);
			}

			if (abs(u_R1[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_R1[i][j] - u_Gave);
			}
			if (abs(u_R2[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_R2[i][j] - u_Gave);
			}
			if (abs(u_R3[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_R3[i][j] - u_Gave);
			}
			if (abs(u_R4[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_R4[i][j] - u_Gave);
			}
			if (abs(u_R5[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_R5[i][j] - u_Gave);
			}

			if (abs(u_B1[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_B1[i][j] - u_Gave);
			}
			if (abs(u_B2[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_B2[i][j] - u_Gave);
			}
			if (abs(u_B3[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_B3[i][j] - u_Gave);
			}
			if (abs(u_B4[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_B4[i][j] - u_Gave);
			}
			if (abs(u_B5[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_B5[i][j] - u_Gave);
			}

			if (abs(u_T1[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_T1[i][j] - u_Gave);
			}
			if (abs(u_T2[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_T2[i][j] - u_Gave);
			}
			if (abs(u_T3[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_T3[i][j] - u_Gave);
			}
			if (abs(u_T4[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_T4[i][j] - u_Gave);
			}
			if (abs(u_T5[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_T5[i][j] - u_Gave);
			}

			if (abs(u_GQ11[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ11[i][j] - u_Gave);
			}
			if (abs(u_GQ12[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ12[i][j] - u_Gave);
			}
			if (abs(u_GQ13[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ13[i][j] - u_Gave);
			}
			if (abs(u_GQ14[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ14[i][j] - u_Gave);
			}
			if (abs(u_GQ15[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ15[i][j] - u_Gave);
			}

			if (abs(u_GQ21[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ21[i][j] - u_Gave);
			}
			if (abs(u_GQ22[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ22[i][j] - u_Gave);
			}
			if (abs(u_GQ23[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ23[i][j] - u_Gave);
			}
			if (abs(u_GQ24[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ24[i][j] - u_Gave);
			}
			if (abs(u_GQ25[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ25[i][j] - u_Gave);
			}

			if (abs(u_GQ31[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ31[i][j] - u_Gave);
			}
			if (abs(u_GQ32[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ32[i][j] - u_Gave);
			}
			if (abs(u_GQ33[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ33[i][j] - u_Gave);
			}
			if (abs(u_GQ34[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ34[i][j] - u_Gave);
			}
			if (abs(u_GQ35[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ35[i][j] - u_Gave);
			}

			if (abs(u_GQ41[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ41[i][j] - u_Gave);
			}
			if (abs(u_GQ42[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ42[i][j] - u_Gave);
			}
			if (abs(u_GQ43[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ43[i][j] - u_Gave);
			}
			if (abs(u_GQ44[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ44[i][j] - u_Gave);
			}
			if (abs(u_GQ45[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ45[i][j] - u_Gave);
			}

			if (abs(u_GQ51[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ51[i][j] - u_Gave);
			}
			if (abs(u_GQ52[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ52[i][j] - u_Gave);
			}
			if (abs(u_GQ53[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ53[i][j] - u_Gave);
			}
			if (abs(u_GQ54[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ54[i][j] - u_Gave);
			}
			if (abs(u_GQ55[i][j] - u_Gave) > u_GLinf)
			{
				u_GLinf = abs(u_GQ55[i][j] - u_Gave);
			}
		}
	}

	Linf[0] = u_GLinf;

	//cout << u_Gave << "  " << u_GLinf << endl;

	if (u_GLinf < pow(10, -12))
	{
		for (int i = 0; i < N_y; i++)
		{
			for (int j = 0; j < N_x; j++)
			{
				sigma_0_x[i][j] = 0.0; sigma_1_x[i][j] = 0.0; sigma_2_x[i][j] = 0.0;
				sigma_0_y[i][j] = 0.0; sigma_1_y[i][j] = 0.0; sigma_2_y[i][j] = 0.0;
			}
		}
	}
	else
	{
		double** ux_jump1; createMatrix(N_y, N_x + 3, ux_jump1); double** ux_jump2; createMatrix(N_y, N_x + 3, ux_jump2);
		double** ux_jump3; createMatrix(N_y, N_x + 3, ux_jump3); double** ux_jump4; createMatrix(N_y, N_x + 3, ux_jump4);
		double** ux_jump5; createMatrix(N_y, N_x + 3, ux_jump5); double** sigmax_0; createMatrix(N_y, N_x + 1, sigmax_0);

		double** uy_jump1; createMatrix(N_y + 3, N_x, uy_jump1); double** uy_jump2; createMatrix(N_y + 3, N_x, uy_jump2);
		double** uy_jump3; createMatrix(N_y + 3, N_x, uy_jump3); double** uy_jump4; createMatrix(N_y + 3, N_x, uy_jump4);
		double** uy_jump5; createMatrix(N_y + 3, N_x, uy_jump5); double** sigmay_0; createMatrix(N_y + 1, N_x, sigmay_0);

		compute_jump_x(ux_jump1, u_L1, u_R1); compute_jump_x(ux_jump2, u_L2, u_R2);
		compute_jump_x(ux_jump3, u_L3, u_R3); compute_jump_x(ux_jump4, u_L4, u_R4); compute_jump_x(ux_jump5, u_L5, u_R5);

		compute_jump_y(uy_jump1, u_B1, u_T1); compute_jump_y(uy_jump2, u_B2, u_T2);
		compute_jump_y(uy_jump3, u_B3, u_T3); compute_jump_y(uy_jump4, u_B4, u_T4); compute_jump_y(uy_jump5, u_B5, u_T5);


		double** ux_xi_plus1; createMatrix(N_y + 2, N_x + 2, ux_xi_plus1); double** ux_xi_plus2; createMatrix(N_y + 2, N_x + 2, ux_xi_plus2);
		double** ux_xi_plus3; createMatrix(N_y + 2, N_x + 2, ux_xi_plus3); double** ux_xi_plus4; createMatrix(N_y + 2, N_x + 2, ux_xi_plus4);
		double** ux_xi_plus5; createMatrix(N_y + 2, N_x + 2, ux_xi_plus5);

		double** ux_xi_minus1; createMatrix(N_y + 2, N_x + 2, ux_xi_minus1); double** ux_xi_minus2; createMatrix(N_y + 2, N_x + 2, ux_xi_minus2);
		double** ux_xi_minus3; createMatrix(N_y + 2, N_x + 2, ux_xi_minus3); double** ux_xi_minus4; createMatrix(N_y + 2, N_x + 2, ux_xi_minus4);
		double** ux_xi_minus5; createMatrix(N_y + 2, N_x + 2, ux_xi_minus5);

		double** ux_eta_plus1; createMatrix(N_y + 2, N_x + 2, ux_eta_plus1); double** ux_eta_plus2; createMatrix(N_y + 2, N_x + 2, ux_eta_plus2);
		double** ux_eta_plus3; createMatrix(N_y + 2, N_x + 2, ux_eta_plus3); double** ux_eta_plus4; createMatrix(N_y + 2, N_x + 2, ux_eta_plus4);
		double** ux_eta_plus5; createMatrix(N_y + 2, N_x + 2, ux_eta_plus5);

		double** ux_eta_minus1; createMatrix(N_y + 2, N_x + 2, ux_eta_minus1); double** ux_eta_minus2; createMatrix(N_y + 2, N_x + 2, ux_eta_minus2);
		double** ux_eta_minus3; createMatrix(N_y + 2, N_x + 2, ux_eta_minus3); double** ux_eta_minus4; createMatrix(N_y + 2, N_x + 2, ux_eta_minus4);
		double** ux_eta_minus5; createMatrix(N_y + 2, N_x + 2, ux_eta_minus5);

		compute_upxi(ux_xi_plus1, -1.0, GQxi_1, u_1, u_3, u_4); compute_upxi(ux_xi_plus2, -1.0, GQxi_2, u_1, u_3, u_4);
		compute_upxi(ux_xi_plus3, -1.0, GQxi_3, u_1, u_3, u_4); compute_upxi(ux_xi_plus4, -1.0, GQxi_4, u_1, u_3, u_4);
		compute_upxi(ux_xi_plus5, -1.0, GQxi_5, u_1, u_3, u_4);

		compute_upxi(ux_xi_minus1, 1.0, GQxi_1, u_1, u_3, u_4); compute_upxi(ux_xi_minus2, 1.0, GQxi_2, u_1, u_3, u_4);
		compute_upxi(ux_xi_minus3, 1.0, GQxi_3, u_1, u_3, u_4); compute_upxi(ux_xi_minus4, 1.0, GQxi_4, u_1, u_3, u_4);
		compute_upxi(ux_xi_minus5, 1.0, GQxi_5, u_1, u_3, u_4);

		compute_upeta(ux_eta_plus1, -1.0, GQxi_1, u_2, u_3, u_5); compute_upeta(ux_eta_plus2, -1.0, GQxi_2, u_2, u_3, u_5);
		compute_upeta(ux_eta_plus3, -1.0, GQxi_3, u_2, u_3, u_5); compute_upeta(ux_eta_plus4, -1.0, GQxi_4, u_2, u_3, u_5);
		compute_upeta(ux_eta_plus5, -1.0, GQxi_5, u_2, u_3, u_5);

		compute_upeta(ux_eta_minus1, 1.0, GQxi_1, u_2, u_3, u_5); compute_upeta(ux_eta_minus2, 1.0, GQxi_2, u_2, u_3, u_5);
		compute_upeta(ux_eta_minus3, 1.0, GQxi_3, u_2, u_3, u_5); compute_upeta(ux_eta_minus4, 1.0, GQxi_4, u_2, u_3, u_5);
		compute_upeta(ux_eta_minus5, 1.0, GQxi_5, u_2, u_3, u_5);

		double** ux_xi_jump1; createMatrix(N_y, N_x + 3, ux_xi_jump1);
		double** ux_xi_jump2; createMatrix(N_y, N_x + 3, ux_xi_jump2);
		double** ux_xi_jump3; createMatrix(N_y, N_x + 3, ux_xi_jump3);
		double** ux_xi_jump4; createMatrix(N_y, N_x + 3, ux_xi_jump4);
		double** ux_xi_jump5; createMatrix(N_y, N_x + 3, ux_xi_jump5);

		double** ux_eta_jump1; createMatrix(N_y, N_x + 3, ux_eta_jump1);
		double** ux_eta_jump2; createMatrix(N_y, N_x + 3, ux_eta_jump2);
		double** ux_eta_jump3; createMatrix(N_y, N_x + 3, ux_eta_jump3);
		double** ux_eta_jump4; createMatrix(N_y, N_x + 3, ux_eta_jump4);
		double** ux_eta_jump5; createMatrix(N_y, N_x + 3, ux_eta_jump5);

		compute_jump_x(ux_xi_jump1, ux_xi_plus1, ux_xi_minus1);
		compute_jump_x(ux_xi_jump2, ux_xi_plus2, ux_xi_minus2);
		compute_jump_x(ux_xi_jump3, ux_xi_plus3, ux_xi_minus3);
		compute_jump_x(ux_xi_jump4, ux_xi_plus4, ux_xi_minus4);
		compute_jump_x(ux_xi_jump5, ux_xi_plus5, ux_xi_minus5);

		compute_jump_x(ux_eta_jump1, ux_eta_plus1, ux_eta_minus1);
		compute_jump_x(ux_eta_jump2, ux_eta_plus2, ux_eta_minus2);
		compute_jump_x(ux_eta_jump3, ux_eta_plus3, ux_eta_minus3);
		compute_jump_x(ux_eta_jump4, ux_eta_plus4, ux_eta_minus4);
		compute_jump_x(ux_eta_jump5, ux_eta_plus5, ux_eta_minus5);

		double** uy_xi_plus1; createMatrix(N_y + 2, N_x + 2, uy_xi_plus1); double** uy_xi_plus2; createMatrix(N_y + 2, N_x + 2, uy_xi_plus2);
		double** uy_xi_plus3; createMatrix(N_y + 2, N_x + 2, uy_xi_plus3); double** uy_xi_plus4; createMatrix(N_y + 2, N_x + 2, uy_xi_plus4);
		double** uy_xi_plus5; createMatrix(N_y + 2, N_x + 2, uy_xi_plus5);

		double** uy_xi_minus1; createMatrix(N_y + 2, N_x + 2, uy_xi_minus1); double** uy_xi_minus2; createMatrix(N_y + 2, N_x + 2, uy_xi_minus2);
		double** uy_xi_minus3; createMatrix(N_y + 2, N_x + 2, uy_xi_minus3); double** uy_xi_minus4; createMatrix(N_y + 2, N_x + 2, uy_xi_minus4);
		double** uy_xi_minus5; createMatrix(N_y + 2, N_x + 2, uy_xi_minus5);

		double** uy_eta_plus1; createMatrix(N_y + 2, N_x + 2, uy_eta_plus1); double** uy_eta_plus2; createMatrix(N_y + 2, N_x + 2, uy_eta_plus2);
		double** uy_eta_plus3; createMatrix(N_y + 2, N_x + 2, uy_eta_plus3); double** uy_eta_plus4; createMatrix(N_y + 2, N_x + 2, uy_eta_plus4);
		double** uy_eta_plus5; createMatrix(N_y + 2, N_x + 2, uy_eta_plus5);

		double** uy_eta_minus1; createMatrix(N_y + 2, N_x + 2, uy_eta_minus1); double** uy_eta_minus2; createMatrix(N_y + 2, N_x + 2, uy_eta_minus2);
		double** uy_eta_minus3; createMatrix(N_y + 2, N_x + 2, uy_eta_minus3); double** uy_eta_minus4; createMatrix(N_y + 2, N_x + 2, uy_eta_minus4);
		double** uy_eta_minus5; createMatrix(N_y + 2, N_x + 2, uy_eta_minus5);

		compute_upxi(uy_xi_plus1, GQxi_1, -1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_plus2, GQxi_2, -1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_plus3, GQxi_3, -1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_plus4, GQxi_4, -1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_plus5, GQxi_5, -1.0, u_1, u_3, u_4);

		compute_upxi(uy_xi_minus1, GQxi_1, 1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_minus2, GQxi_2, 1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_minus3, GQxi_3, 1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_minus4, GQxi_4, 1.0, u_1, u_3, u_4);
		compute_upxi(uy_xi_minus5, GQxi_5, 1.0, u_1, u_3, u_4);

		compute_upeta(uy_eta_plus1, GQxi_1, -1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_plus2, GQxi_2, -1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_plus3, GQxi_3, -1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_plus4, GQxi_4, -1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_plus5, GQxi_5, -1.0, u_2, u_3, u_5);

		compute_upeta(uy_eta_minus1, GQxi_1, 1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_minus2, GQxi_2, 1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_minus3, GQxi_3, 1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_minus4, GQxi_4, 1.0, u_2, u_3, u_5);
		compute_upeta(uy_eta_minus5, GQxi_5, 1.0, u_2, u_3, u_5);

		double** uy_xi_jump1; createMatrix(N_y + 3, N_x, uy_xi_jump1);
		double** uy_xi_jump2; createMatrix(N_y + 3, N_x, uy_xi_jump2);
		double** uy_xi_jump3; createMatrix(N_y + 3, N_x, uy_xi_jump3);
		double** uy_xi_jump4; createMatrix(N_y + 3, N_x, uy_xi_jump4);
		double** uy_xi_jump5; createMatrix(N_y + 3, N_x, uy_xi_jump5);

		double** uy_eta_jump1; createMatrix(N_y + 3, N_x, uy_eta_jump1);
		double** uy_eta_jump2; createMatrix(N_y + 3, N_x, uy_eta_jump2);
		double** uy_eta_jump3; createMatrix(N_y + 3, N_x, uy_eta_jump3);
		double** uy_eta_jump4; createMatrix(N_y + 3, N_x, uy_eta_jump4);
		double** uy_eta_jump5; createMatrix(N_y + 3, N_x, uy_eta_jump5);

		compute_jump_y(uy_xi_jump1, uy_xi_plus1, uy_xi_minus1);
		compute_jump_y(uy_xi_jump2, uy_xi_plus2, uy_xi_minus2);
		compute_jump_y(uy_xi_jump3, uy_xi_plus3, uy_xi_minus3);
		compute_jump_y(uy_xi_jump4, uy_xi_plus4, uy_xi_minus4);
		compute_jump_y(uy_xi_jump5, uy_xi_plus5, uy_xi_minus5);

		compute_jump_y(uy_eta_jump1, uy_eta_plus1, uy_eta_minus1);
		compute_jump_y(uy_eta_jump2, uy_eta_plus2, uy_eta_minus2);
		compute_jump_y(uy_eta_jump3, uy_eta_plus3, uy_eta_minus3);
		compute_jump_y(uy_eta_jump4, uy_eta_plus4, uy_eta_minus4);
		compute_jump_y(uy_eta_jump5, uy_eta_plus5, uy_eta_minus5);

		double** sigmax_1; createMatrix(N_y, N_x + 1, sigmax_1);
		double** sigmay_1; createMatrix(N_y + 1, N_x, sigmay_1);

		double** u_xi2_plus; createMatrix(N_y + 2, N_x + 2, u_xi2_plus);
		double** u_xi2_minus; createMatrix(N_y + 2, N_x + 2, u_xi2_minus);
		double** u_xieta_plus; createMatrix(N_y + 2, N_x + 2, u_xieta_plus);
		double** u_xieta_minus; createMatrix(N_y + 2, N_x + 2, u_xieta_minus);
		double** u_eta2_plus; createMatrix(N_y + 2, N_x + 2, u_eta2_plus);
		double** u_eta2_minus; createMatrix(N_y + 2, N_x + 2, u_eta2_minus);

		for (int i = 0; i < N_y + 2; i++)
		{
			for (int j = 0; j < N_x + 2; j++)
			{
				u_xi2_plus[i][j] = 3.0 * u_4[i][j]; u_xi2_minus[i][j] = u_xi2_plus[i][j];
				u_xieta_plus[i][j] = u_3[i][j]; u_xieta_minus[i][j] = u_xieta_plus[i][j];
				u_eta2_plus[i][j] = 3.0 * u_5[i][j]; u_eta2_minus[i][j] = u_eta2_plus[i][j];
			}
		}

		double** ux_xi2_jump; createMatrix(N_y, N_x + 3, ux_xi2_jump);
		double** ux_xieta_jump; createMatrix(N_y, N_x + 3, ux_xieta_jump);
		double** ux_eta2_jump; createMatrix(N_y, N_x + 3, ux_eta2_jump);

		double** uy_xi2_jump; createMatrix(N_y + 3, N_x, uy_xi2_jump);
		double** uy_xieta_jump; createMatrix(N_y + 3, N_x, uy_xieta_jump);
		double** uy_eta2_jump; createMatrix(N_y + 3, N_x, uy_eta2_jump);

		compute_jump_x(ux_xi2_jump, u_xi2_plus, u_xi2_minus);
		compute_jump_x(ux_xieta_jump, u_xieta_plus, u_xieta_minus);
		compute_jump_x(ux_eta2_jump, u_eta2_plus, u_eta2_minus);

		compute_jump_y(uy_xi2_jump, u_xi2_plus, u_xi2_minus);
		compute_jump_y(uy_xieta_jump, u_xieta_plus, u_xieta_minus);
		compute_jump_y(uy_eta2_jump, u_eta2_plus, u_eta2_minus);

		double** sigmax_2; createMatrix(N_y, N_x + 1, sigmax_2);
		double** sigmay_2; createMatrix(N_y + 1, N_x, sigmay_2);

		for (int i = 0; i < N_y; i++)
		{
			for (int j = 0; j < N_x + 1; j++)
			{
				sigmax_0[i][j] = 1.0 / 12.0 * (weight_1 * abs(ux_jump1[i][j + 1]) + weight_2 * abs(ux_jump2[i][j + 1]) + weight_3 * abs(ux_jump3[i][j + 1]) + weight_4 * abs(ux_jump4[i][j + 1]) + weight_5 * abs(ux_jump5[i][j + 1]));

				double temp_11 = weight_1 * abs(ux_xi_jump1[i][j + 1]) + weight_2 * abs(ux_xi_jump2[i][j + 1]) + weight_3 * abs(ux_xi_jump3[i][j + 1]) + weight_4 * abs(ux_xi_jump4[i][j + 1]) + weight_5 * abs(ux_xi_jump5[i][j + 1]);
				double temp_12 = weight_1 * abs(ux_eta_jump1[i][j + 1]) + weight_2 * abs(ux_eta_jump2[i][j + 1]) + weight_3 * abs(ux_eta_jump3[i][j + 1]) + weight_4 * abs(ux_eta_jump4[i][j + 1]) + weight_5 * abs(ux_eta_jump5[i][j + 1]);
				sigmax_1[i][j] = 0.5 * temp_11 + 0.5 * h_x / h_y * temp_12;

				sigmax_2[i][j] = 5.0 / 3.0 * abs(ux_xi2_jump[i][j + 1]) + 5.0 * h_x / (3.0 * h_y) * abs(ux_xieta_jump[i][j + 1]) + 5.0 * h_x * h_x / (3.0 * h_y * h_y) * abs(ux_eta2_jump[i][j + 1]);
			}
		}

		for (int j = 0; j < N_x; j++)
		{
			for (int i = 0; i < N_y + 1; i++)
			{
				sigmay_0[i][j] = 1.0 / 12.0 * (weight_1 * abs(uy_jump1[i + 1][j]) + weight_2 * abs(uy_jump2[i + 1][j]) + weight_3 * abs(uy_jump3[i + 1][j]) + weight_4 * abs(uy_jump4[i + 1][j]) + weight_5 * abs(uy_jump5[i + 1][j]));

				double temp_11 = weight_1 * abs(uy_xi_jump1[i + 1][j]) + weight_2 * abs(uy_xi_jump2[i + 1][j]) + weight_3 * abs(uy_xi_jump3[i + 1][j]) + weight_4 * abs(uy_xi_jump4[i + 1][j]) + weight_5 * abs(uy_xi_jump5[i + 1][j]);
				double temp_12 = weight_1 * abs(uy_eta_jump1[i + 1][j]) + weight_2 * abs(uy_eta_jump2[i + 1][j]) + weight_3 * abs(uy_eta_jump3[i + 1][j]) + weight_4 * abs(uy_eta_jump4[i + 1][j]) + weight_5 * abs(uy_eta_jump5[i + 1][j]);
				sigmay_1[i][j] = 0.5 * h_y / h_x * temp_11 + 0.5 * temp_12;

				sigmay_2[i][j] = 5.0 * h_y * h_y / (3.0 * h_x * h_x) * abs(uy_xi2_jump[i + 1][j]) + 5.0 * h_y / (3.0 * h_x) * abs(uy_xieta_jump[i + 1][j]) + 5.0 / 3.0 * abs(uy_eta2_jump[i + 1][j]);

			}
		}

		for (int i = 0; i < N_y; i++)
		{
			for (int j = 0; j < N_x; j++)
			{
				sigma_0_x[i][j] = Omega * (sigmax_0[i][j] + sigmax_0[i][j + 1]);
				sigma_1_x[i][j] = Omega * (sigmax_1[i][j] + sigmax_1[i][j + 1]);
				sigma_2_x[i][j] = Omega * (sigmax_2[i][j] + sigmax_2[i][j + 1]);
			}
		}

		for (int j = 0; j < N_x; j++)
		{
			for (int i = 0; i < N_y; i++)
			{
				sigma_0_y[i][j] = Omega * (sigmay_0[i][j] + sigmay_0[i + 1][j]);
				sigma_1_y[i][j] = Omega * (sigmay_1[i][j] + sigmay_1[i + 1][j]);
				sigma_2_y[i][j] = Omega * (sigmay_2[i][j] + sigmay_2[i + 1][j]);
			}
		}

		deleteMatrix(N_y, ux_jump1); deleteMatrix(N_y, ux_jump2); deleteMatrix(N_y, ux_jump3);
		deleteMatrix(N_y, ux_jump4); deleteMatrix(N_y, ux_jump5); deleteMatrix(N_y, sigmax_0);
		deleteMatrix(N_y + 3, uy_jump1); deleteMatrix(N_y + 3, uy_jump2); deleteMatrix(N_y + 3, uy_jump3);
		deleteMatrix(N_y + 3, uy_jump4); deleteMatrix(N_y + 3, uy_jump5); deleteMatrix(N_y + 1, sigmay_0);

		deleteMatrix(N_y + 2, ux_xi_plus1); deleteMatrix(N_y + 2, ux_xi_plus2); deleteMatrix(N_y + 2, ux_xi_plus3);
		deleteMatrix(N_y + 2, ux_xi_plus4); deleteMatrix(N_y + 2, ux_xi_plus5);
		deleteMatrix(N_y + 2, ux_xi_minus1); deleteMatrix(N_y + 2, ux_xi_minus2); deleteMatrix(N_y + 2, ux_xi_minus3);
		deleteMatrix(N_y + 2, ux_xi_minus4); deleteMatrix(N_y + 2, ux_xi_minus5);
		deleteMatrix(N_y + 2, ux_eta_plus1); deleteMatrix(N_y + 2, ux_eta_plus2); deleteMatrix(N_y + 2, ux_eta_plus3);
		deleteMatrix(N_y + 2, ux_eta_plus4); deleteMatrix(N_y + 2, ux_eta_plus5);
		deleteMatrix(N_y + 2, ux_eta_minus1); deleteMatrix(N_y + 2, ux_eta_minus2); deleteMatrix(N_y + 2, ux_eta_minus3);
		deleteMatrix(N_y + 2, ux_eta_minus4); deleteMatrix(N_y + 2, ux_eta_minus5);

		deleteMatrix(N_y, ux_xi_jump1); deleteMatrix(N_y, ux_xi_jump2); deleteMatrix(N_y, ux_xi_jump3);
		deleteMatrix(N_y, ux_xi_jump4); deleteMatrix(N_y, ux_xi_jump5);

		deleteMatrix(N_y, ux_eta_jump1); deleteMatrix(N_y, ux_eta_jump2); deleteMatrix(N_y, ux_eta_jump3);
		deleteMatrix(N_y, ux_eta_jump4); deleteMatrix(N_y, ux_eta_jump5);

		deleteMatrix(N_y + 2, uy_xi_plus1); deleteMatrix(N_y + 2, uy_xi_plus2); deleteMatrix(N_y + 2, uy_xi_plus3);
		deleteMatrix(N_y + 2, uy_xi_plus4); deleteMatrix(N_y + 2, uy_xi_plus5);
		deleteMatrix(N_y + 2, uy_xi_minus1); deleteMatrix(N_y + 2, uy_xi_minus2); deleteMatrix(N_y + 2, uy_xi_minus3);
		deleteMatrix(N_y + 2, uy_xi_minus4); deleteMatrix(N_y + 2, uy_xi_minus5);
		deleteMatrix(N_y + 2, uy_eta_plus1); deleteMatrix(N_y + 2, uy_eta_plus2); deleteMatrix(N_y + 2, uy_eta_plus3);
		deleteMatrix(N_y + 2, uy_eta_plus4); deleteMatrix(N_y + 2, uy_eta_plus5);
		deleteMatrix(N_y + 2, uy_eta_minus1); deleteMatrix(N_y + 2, uy_eta_minus2); deleteMatrix(N_y + 2, uy_eta_minus3);
		deleteMatrix(N_y + 2, uy_eta_minus4); deleteMatrix(N_y + 2, uy_eta_minus5);

		deleteMatrix(N_y + 3, uy_xi_jump1); deleteMatrix(N_y + 3, uy_xi_jump2); deleteMatrix(N_y + 3, uy_xi_jump3);
		deleteMatrix(N_y + 3, uy_xi_jump4); deleteMatrix(N_y + 3, uy_xi_jump5);

		deleteMatrix(N_y + 3, uy_eta_jump1); deleteMatrix(N_y + 3, uy_eta_jump2); deleteMatrix(N_y + 3, uy_eta_jump3);
		deleteMatrix(N_y + 3, uy_eta_jump4); deleteMatrix(N_y + 3, uy_eta_jump5);

		deleteMatrix(N_y, sigmax_1); deleteMatrix(N_y + 1, sigmay_1);

		deleteMatrix(N_y + 2, u_xi2_plus); deleteMatrix(N_y + 2, u_xi2_minus);
		deleteMatrix(N_y + 2, u_xieta_plus); deleteMatrix(N_y + 2, u_xieta_minus);
		deleteMatrix(N_y + 2, u_eta2_plus); deleteMatrix(N_y + 2, u_eta2_minus);

		deleteMatrix(N_y, ux_xi2_jump); deleteMatrix(N_y, ux_xieta_jump); deleteMatrix(N_y, ux_eta2_jump);
		deleteMatrix(N_y + 3, uy_xi2_jump); deleteMatrix(N_y + 3, uy_xieta_jump); deleteMatrix(N_y + 3, uy_eta2_jump);

		deleteMatrix(N_y, sigmax_2); deleteMatrix(N_y + 1, sigmay_2);
	}

	deleteMatrix(N_y + 2, u_GQ11); deleteMatrix(N_y + 2, u_GQ12); deleteMatrix(N_y + 2, u_GQ13); deleteMatrix(N_y + 2, u_GQ14); deleteMatrix(N_y + 2, u_GQ15);
	deleteMatrix(N_y + 2, u_GQ21); deleteMatrix(N_y + 2, u_GQ22); deleteMatrix(N_y + 2, u_GQ23); deleteMatrix(N_y + 2, u_GQ24); deleteMatrix(N_y + 2, u_GQ25);
	deleteMatrix(N_y + 2, u_GQ31); deleteMatrix(N_y + 2, u_GQ32); deleteMatrix(N_y + 2, u_GQ33); deleteMatrix(N_y + 2, u_GQ34); deleteMatrix(N_y + 2, u_GQ35);
	deleteMatrix(N_y + 2, u_GQ41); deleteMatrix(N_y + 2, u_GQ42); deleteMatrix(N_y + 2, u_GQ43); deleteMatrix(N_y + 2, u_GQ44); deleteMatrix(N_y + 2, u_GQ45);
	deleteMatrix(N_y + 2, u_GQ51); deleteMatrix(N_y + 2, u_GQ52); deleteMatrix(N_y + 2, u_GQ53); deleteMatrix(N_y + 2, u_GQ54); deleteMatrix(N_y + 2, u_GQ55);

	deleteMatrix(N_y + 2, u_L1); deleteMatrix(N_y + 2, u_L2); deleteMatrix(N_y + 2, u_L3); deleteMatrix(N_y + 2, u_L4); deleteMatrix(N_y + 2, u_L5);
	deleteMatrix(N_y + 2, u_R1); deleteMatrix(N_y + 2, u_R2); deleteMatrix(N_y + 2, u_R3); deleteMatrix(N_y + 2, u_R4); deleteMatrix(N_y + 2, u_R5);
	deleteMatrix(N_y + 2, u_B1); deleteMatrix(N_y + 2, u_B2); deleteMatrix(N_y + 2, u_B3); deleteMatrix(N_y + 2, u_B4); deleteMatrix(N_y + 2, u_B5);
	deleteMatrix(N_y + 2, u_T1); deleteMatrix(N_y + 2, u_T2); deleteMatrix(N_y + 2, u_T3); deleteMatrix(N_y + 2, u_T4); deleteMatrix(N_y + 2, u_T5);
}

void OE_compute_system_sigma(double** sigma_0, double** sigma_1, double** sigma_2, double** rho_0, double** rho_1, double** rho_2, double** rho_3, double** rho_4, double** rho_5, double** m1_0, double** m1_1, double** m1_2, double** m1_3, double** m1_4, double** m1_5, double** m2_0, double** m2_1, double** m2_2, double** m2_3, double** m2_4, double** m2_5, double** Ee_0, double** Ee_1, double** Ee_2, double** Ee_3, double** Ee_4, double** Ee_5, double** Ei_0, double** Ei_1, double** Ei_2, double** Ei_3, double** Ei_4, double** Ei_5, double** Er_0, double** Er_1, double** Er_2, double** Er_3, double** Er_4, double** Er_5, double h_x, double h_y)
{
	double** rho_sigma_0; createMatrix(N_y, N_x, rho_sigma_0); double** rho_sigma_1; createMatrix(N_y, N_x, rho_sigma_1);
	double** rho_sigma_2; createMatrix(N_y, N_x, rho_sigma_2);
	double** m1_sigma_0; createMatrix(N_y, N_x, m1_sigma_0); double** m1_sigma_1; createMatrix(N_y, N_x, m1_sigma_1);
	double** m1_sigma_2; createMatrix(N_y, N_x, m1_sigma_2);
	double** m2_sigma_0; createMatrix(N_y, N_x, m2_sigma_0); double** m2_sigma_1; createMatrix(N_y, N_x, m2_sigma_1);
	double** m2_sigma_2; createMatrix(N_y, N_x, m2_sigma_2);
	double** Ee_sigma_0; createMatrix(N_y, N_x, Ee_sigma_0); double** Ee_sigma_1; createMatrix(N_y, N_x, Ee_sigma_1);
	double** Ee_sigma_2; createMatrix(N_y, N_x, Ee_sigma_2);
	double** Ei_sigma_0; createMatrix(N_y, N_x, Ei_sigma_0); double** Ei_sigma_1; createMatrix(N_y, N_x, Ei_sigma_1);
	double** Ei_sigma_2; createMatrix(N_y, N_x, Ei_sigma_2);
	double** Er_sigma_0; createMatrix(N_y, N_x, Er_sigma_0); double** Er_sigma_1; createMatrix(N_y, N_x, Er_sigma_1);
	double** Er_sigma_2; createMatrix(N_y, N_x, Er_sigma_2);

	double** rho_sigma_0_x; createMatrix(N_y, N_x, rho_sigma_0_x); double** rho_sigma_1_x; createMatrix(N_y, N_x, rho_sigma_1_x);
	double** rho_sigma_2_x; createMatrix(N_y, N_x, rho_sigma_2_x); double** rho_sigma_0_y; createMatrix(N_y, N_x, rho_sigma_0_y);
	double** rho_sigma_1_y; createMatrix(N_y, N_x, rho_sigma_1_y); double** rho_sigma_2_y; createMatrix(N_y, N_x, rho_sigma_2_y);

	double** m1_sigma_0_x; createMatrix(N_y, N_x, m1_sigma_0_x); double** m1_sigma_1_x; createMatrix(N_y, N_x, m1_sigma_1_x);
	double** m1_sigma_2_x; createMatrix(N_y, N_x, m1_sigma_2_x); double** m1_sigma_0_y; createMatrix(N_y, N_x, m1_sigma_0_y);
	double** m1_sigma_1_y; createMatrix(N_y, N_x, m1_sigma_1_y); double** m1_sigma_2_y; createMatrix(N_y, N_x, m1_sigma_2_y);

	double** m2_sigma_0_x; createMatrix(N_y, N_x, m2_sigma_0_x); double** m2_sigma_1_x; createMatrix(N_y, N_x, m2_sigma_1_x);
	double** m2_sigma_2_x; createMatrix(N_y, N_x, m2_sigma_2_x); double** m2_sigma_0_y; createMatrix(N_y, N_x, m2_sigma_0_y);
	double** m2_sigma_1_y; createMatrix(N_y, N_x, m2_sigma_1_y); double** m2_sigma_2_y; createMatrix(N_y, N_x, m2_sigma_2_y);

	double** Ee_sigma_0_x; createMatrix(N_y, N_x, Ee_sigma_0_x); double** Ee_sigma_1_x; createMatrix(N_y, N_x, Ee_sigma_1_x);
	double** Ee_sigma_2_x; createMatrix(N_y, N_x, Ee_sigma_2_x); double** Ee_sigma_0_y; createMatrix(N_y, N_x, Ee_sigma_0_y);
	double** Ee_sigma_1_y; createMatrix(N_y, N_x, Ee_sigma_1_y); double** Ee_sigma_2_y; createMatrix(N_y, N_x, Ee_sigma_2_y);

	double** Ei_sigma_0_x; createMatrix(N_y, N_x, Ei_sigma_0_x); double** Ei_sigma_1_x; createMatrix(N_y, N_x, Ei_sigma_1_x);
	double** Ei_sigma_2_x; createMatrix(N_y, N_x, Ei_sigma_2_x); double** Ei_sigma_0_y; createMatrix(N_y, N_x, Ei_sigma_0_y);
	double** Ei_sigma_1_y; createMatrix(N_y, N_x, Ei_sigma_1_y); double** Ei_sigma_2_y; createMatrix(N_y, N_x, Ei_sigma_2_y);

	double** Er_sigma_0_x; createMatrix(N_y, N_x, Er_sigma_0_x); double** Er_sigma_1_x; createMatrix(N_y, N_x, Er_sigma_1_x);
	double** Er_sigma_2_x; createMatrix(N_y, N_x, Er_sigma_2_x); double** Er_sigma_0_y; createMatrix(N_y, N_x, Er_sigma_0_y);
	double** Er_sigma_1_y; createMatrix(N_y, N_x, Er_sigma_1_y); double** Er_sigma_2_y; createMatrix(N_y, N_x, Er_sigma_2_y);

	double** beta_ij_x; createMatrix(N_y, N_x, beta_ij_x); double** beta_ij_y; createMatrix(N_y, N_x, beta_ij_y);

	double* rhoL_inf = new double[1]; double* m1L_inf = new double[1]; double* m2L_inf = new double[1];
	double* EeL_inf = new double[1]; double* EiL_inf = new double[1]; double* ErL_inf = new double[1];

	OE_compute_scalar_sigma(rhoL_inf, rho_sigma_0_x, rho_sigma_1_x, rho_sigma_2_x, rho_sigma_0_y, rho_sigma_1_y, rho_sigma_2_y, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, h_x, h_y);
	OE_compute_scalar_sigma(m1L_inf, m1_sigma_0_x, m1_sigma_1_x, m1_sigma_2_x, m1_sigma_0_y, m1_sigma_1_y, m1_sigma_2_y, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x, h_y);
	OE_compute_scalar_sigma(m2L_inf, m2_sigma_0_x, m2_sigma_1_x, m2_sigma_2_x, m2_sigma_0_y, m2_sigma_1_y, m2_sigma_2_y, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_x, h_y);
	OE_compute_scalar_sigma(EeL_inf, Ee_sigma_0_x, Ee_sigma_1_x, Ee_sigma_2_x, Ee_sigma_0_y, Ee_sigma_1_y, Ee_sigma_2_y, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, h_x, h_y);
	OE_compute_scalar_sigma(EiL_inf, Ei_sigma_0_x, Ei_sigma_1_x, Ei_sigma_2_x, Ei_sigma_0_y, Ei_sigma_1_y, Ei_sigma_2_y, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, h_x, h_y);
	OE_compute_scalar_sigma(ErL_inf, Er_sigma_0_x, Er_sigma_1_x, Er_sigma_2_x, Er_sigma_0_y, Er_sigma_1_y, Er_sigma_2_y, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, h_x, h_y);

	double L_inf = rhoL_inf[0];
	if (m1L_inf[0] > L_inf)
	{
		L_inf = m1L_inf[0];
	}
	if (m2L_inf[0] > L_inf)
	{
		L_inf = m2L_inf[0];
	}
	if (EeL_inf[0] > L_inf)
	{
		L_inf = EeL_inf[0];
	}
	if (EiL_inf[0] > L_inf)
	{
		L_inf = EiL_inf[0];
	}
	if (ErL_inf[0] > L_inf)
	{
		L_inf = ErL_inf[0];
	}

	OE_compute_beta_ij(beta_ij_x, beta_ij_y, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			rho_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * rho_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * rho_sigma_0_y[i][j]) / L_inf;
			rho_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * rho_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * rho_sigma_1_y[i][j]) / L_inf;
			rho_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * rho_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * rho_sigma_2_y[i][j]) / L_inf;

			m1_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * m1_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * m1_sigma_0_y[i][j]) / L_inf;
			m1_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * m1_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * m1_sigma_1_y[i][j]) / L_inf;
			m1_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * m1_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * m1_sigma_2_y[i][j]) / L_inf;

			m2_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * m2_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * m2_sigma_0_y[i][j]) / L_inf;
			m2_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * m2_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * m2_sigma_1_y[i][j]) / L_inf;
			m2_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * m2_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * m2_sigma_2_y[i][j]) / L_inf;

			Ee_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * Ee_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * Ee_sigma_0_y[i][j]) / L_inf;
			Ee_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * Ee_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * Ee_sigma_1_y[i][j]) / L_inf;
			Ee_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * Ee_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * Ee_sigma_2_y[i][j]) / L_inf;

			Ei_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * Ei_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * Ei_sigma_0_y[i][j]) / L_inf;
			Ei_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * Ei_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * Ei_sigma_1_y[i][j]) / L_inf;
			Ei_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * Ei_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * Ei_sigma_2_y[i][j]) / L_inf;

			Er_sigma_0[i][j] = (beta_ij_x[i][j] / h_x * Er_sigma_0_x[i][j] + beta_ij_y[i][j] / h_y * Er_sigma_0_y[i][j]) / L_inf;
			Er_sigma_1[i][j] = (beta_ij_x[i][j] / h_x * Er_sigma_1_x[i][j] + beta_ij_y[i][j] / h_y * Er_sigma_1_y[i][j]) / L_inf;
			Er_sigma_2[i][j] = (beta_ij_x[i][j] / h_x * Er_sigma_2_x[i][j] + beta_ij_y[i][j] / h_y * Er_sigma_2_y[i][j]) / L_inf;

			double temp_max0 = 0.0;
			if (rho_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = rho_sigma_0[i][j];
			}
			if (m1_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = m1_sigma_0[i][j];
			}
			if (m2_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = m2_sigma_0[i][j];
			}
			if (Ee_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = Ee_sigma_0[i][j];
			}
			if (Ei_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = Ei_sigma_0[i][j];
			}
			if (Er_sigma_0[i][j] > temp_max0)
			{
				temp_max0 = Er_sigma_0[i][j];
			}
			sigma_0[i][j] = temp_max0;

			double temp_max1 = 0.0;
			if (rho_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = rho_sigma_1[i][j];
			}
			if (m1_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = m1_sigma_1[i][j];
			}
			if (m2_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = m2_sigma_1[i][j];
			}
			if (Ee_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = Ee_sigma_1[i][j];
			}
			if (Ei_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = Ei_sigma_1[i][j];
			}
			if (Er_sigma_1[i][j] > temp_max1)
			{
				temp_max1 = Er_sigma_1[i][j];
			}
			sigma_1[i][j] = temp_max1;

			double temp_max2 = 0.0;
			if (rho_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = rho_sigma_2[i][j];
			}
			if (m1_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = m1_sigma_2[i][j];
			}
			if (m2_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = m2_sigma_2[i][j];
			}
			if (Ee_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = Ee_sigma_2[i][j];
			}
			if (Ei_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = Ei_sigma_2[i][j];
			}
			if (Er_sigma_2[i][j] > temp_max2)
			{
				temp_max2 = Er_sigma_2[i][j];
			}
			sigma_2[i][j] = temp_max2;
		}
	}

	deleteMatrix(N_y, rho_sigma_0); deleteMatrix(N_y, rho_sigma_1); deleteMatrix(N_y, rho_sigma_2);
	deleteMatrix(N_y, rho_sigma_0_x); deleteMatrix(N_y, rho_sigma_1_x); deleteMatrix(N_y, rho_sigma_2_x);
	deleteMatrix(N_y, rho_sigma_0_y); deleteMatrix(N_y, rho_sigma_1_y); deleteMatrix(N_y, rho_sigma_2_y);
	deleteMatrix(N_y, m1_sigma_0); deleteMatrix(N_y, m1_sigma_1); deleteMatrix(N_y, m1_sigma_2);
	deleteMatrix(N_y, m1_sigma_0_x); deleteMatrix(N_y, m1_sigma_1_x); deleteMatrix(N_y, m1_sigma_2_x);
	deleteMatrix(N_y, m1_sigma_0_y); deleteMatrix(N_y, m1_sigma_1_y); deleteMatrix(N_y, m1_sigma_2_y);
	deleteMatrix(N_y, m2_sigma_0); deleteMatrix(N_y, m2_sigma_1); deleteMatrix(N_y, m2_sigma_2);
	deleteMatrix(N_y, m2_sigma_0_x); deleteMatrix(N_y, m2_sigma_1_x); deleteMatrix(N_y, m2_sigma_2_x);
	deleteMatrix(N_y, m2_sigma_0_y); deleteMatrix(N_y, m2_sigma_1_y); deleteMatrix(N_y, m2_sigma_2_y);
	deleteMatrix(N_y, Ee_sigma_0); deleteMatrix(N_y, Ee_sigma_1); deleteMatrix(N_y, Ee_sigma_2);
	deleteMatrix(N_y, Ee_sigma_0_x); deleteMatrix(N_y, Ee_sigma_1_x); deleteMatrix(N_y, Ee_sigma_2_x);
	deleteMatrix(N_y, Ee_sigma_0_y); deleteMatrix(N_y, Ee_sigma_1_y); deleteMatrix(N_y, Ee_sigma_2_y);
	deleteMatrix(N_y, Ei_sigma_0); deleteMatrix(N_y, Ei_sigma_1); deleteMatrix(N_y, Ei_sigma_2);
	deleteMatrix(N_y, Ei_sigma_0_x); deleteMatrix(N_y, Ei_sigma_1_x); deleteMatrix(N_y, Ei_sigma_2_x);
	deleteMatrix(N_y, Ei_sigma_0_y); deleteMatrix(N_y, Ei_sigma_1_y); deleteMatrix(N_y, Ei_sigma_2_y);
	deleteMatrix(N_y, Er_sigma_0); deleteMatrix(N_y, Er_sigma_1); deleteMatrix(N_y, Er_sigma_2);
	deleteMatrix(N_y, Er_sigma_0_x); deleteMatrix(N_y, Er_sigma_1_x); deleteMatrix(N_y, Er_sigma_2_x);
	deleteMatrix(N_y, Er_sigma_0_y); deleteMatrix(N_y, Er_sigma_1_y); deleteMatrix(N_y, Er_sigma_2_y);
	deleteMatrix(N_y, beta_ij_x); deleteMatrix(N_y, beta_ij_y);
	delete[] rhoL_inf; delete[] m1L_inf; delete[] m2L_inf;
	delete[] EeL_inf; delete[] EiL_inf; delete[] ErL_inf;
}

void OE_procedure(double** u_1, double** u_2, double** u_3, double** u_4, double** u_5, double** sigma_0, double** sigma_1, double** sigma_2, double tau)
{
	for (int i = 0; i < N_y; i++)
	{
		for (int j = 0; j < N_x; j++)
		{
			double C1 = exp(-1.0 * tau * (sigma_0[i][j] + sigma_1[i][j]));
			double C2 = exp(-1.0 * tau * (sigma_0[i][j] + sigma_1[i][j] + sigma_2[i][j]));
			u_1[i + 1][j + 1] = C1 * u_1[i + 1][j + 1];
			u_2[i + 1][j + 1] = C1 * u_2[i + 1][j + 1];

			u_3[i + 1][j + 1] = C2 * u_3[i + 1][j + 1];
			u_4[i + 1][j + 1] = C2 * u_4[i + 1][j + 1];
			u_5[i + 1][j + 1] = C2 * u_5[i + 1][j + 1];
		}

	}
}
