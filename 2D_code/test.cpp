#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include <string.h>
#include "modalDG2D.h"
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
	double x_min = 0.0, x_max = 2 * pi; double h_x = (x_max - x_min) / N_x;
	double y_min = 0.0, y_max = 2 * pi; double h_y = (y_max - y_min) / N_y;

	double* x = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		x[i] = x_min - h_x / 2 + i * h_x;
	}
	double* y = new double[N_y + 2];
	for (int j = 0; j < N_y + 2; j++)
	{
		y[j] = y_min - h_y / 2 + j * h_y;
	}

	// temporal discretization
	double t_ini = 0.0;

	// initialization
	double t = t_ini;
	//L2-projection
	double* X_1 = new double[N_x + 2]; double* X_2 = new double[N_x + 2]; double* X_3 = new double[N_x + 2];
	double* X_4 = new double[N_x + 2]; double* X_5 = new double[N_x + 2];
	get_Gaussian_points_x(X_1, X_2, X_3, X_4, X_5, x, h_x);
	double* Y_1 = new double[N_y + 2]; double* Y_2 = new double[N_y + 2]; double* Y_3 = new double[N_y + 2];
	double* Y_4 = new double[N_y + 2]; double* Y_5 = new double[N_y + 2];
	get_Gaussian_points_y(Y_1, Y_2, Y_3, Y_4, Y_5, y, h_y);

	int mrows = N_y + 2; int mcols = N_x + 2;
	double** rho_11; createMatrix(mrows, mcols, rho_11); double** rho_12; createMatrix(mrows, mcols, rho_12); double** rho_13; createMatrix(mrows, mcols, rho_13); double** rho_14; createMatrix(mrows, mcols, rho_14); double** rho_15; createMatrix(mrows, mcols, rho_15);
	double** rho_21; createMatrix(mrows, mcols, rho_21); double** rho_22; createMatrix(mrows, mcols, rho_22); double** rho_23; createMatrix(mrows, mcols, rho_23); double** rho_24; createMatrix(mrows, mcols, rho_24); double** rho_25; createMatrix(mrows, mcols, rho_25);
	double** rho_31; createMatrix(mrows, mcols, rho_31); double** rho_32; createMatrix(mrows, mcols, rho_32); double** rho_33; createMatrix(mrows, mcols, rho_33); double** rho_34; createMatrix(mrows, mcols, rho_34); double** rho_35; createMatrix(mrows, mcols, rho_35);
	double** rho_41; createMatrix(mrows, mcols, rho_41); double** rho_42; createMatrix(mrows, mcols, rho_42); double** rho_43; createMatrix(mrows, mcols, rho_43); double** rho_44; createMatrix(mrows, mcols, rho_44); double** rho_45; createMatrix(mrows, mcols, rho_45);
	double** rho_51; createMatrix(mrows, mcols, rho_51); double** rho_52; createMatrix(mrows, mcols, rho_52); double** rho_53; createMatrix(mrows, mcols, rho_53); double** rho_54; createMatrix(mrows, mcols, rho_54); double** rho_55; createMatrix(mrows, mcols, rho_55);

	double** m1_11; createMatrix(mrows, mcols, m1_11); double** m1_12; createMatrix(mrows, mcols, m1_12); double** m1_13; createMatrix(mrows, mcols, m1_13); double** m1_14; createMatrix(mrows, mcols, m1_14); double** m1_15; createMatrix(mrows, mcols, m1_15);
	double** m1_21; createMatrix(mrows, mcols, m1_21); double** m1_22; createMatrix(mrows, mcols, m1_22); double** m1_23; createMatrix(mrows, mcols, m1_23); double** m1_24; createMatrix(mrows, mcols, m1_24); double** m1_25; createMatrix(mrows, mcols, m1_25);
	double** m1_31; createMatrix(mrows, mcols, m1_31); double** m1_32; createMatrix(mrows, mcols, m1_32); double** m1_33; createMatrix(mrows, mcols, m1_33); double** m1_34; createMatrix(mrows, mcols, m1_34); double** m1_35; createMatrix(mrows, mcols, m1_35);
	double** m1_41; createMatrix(mrows, mcols, m1_41); double** m1_42; createMatrix(mrows, mcols, m1_42); double** m1_43; createMatrix(mrows, mcols, m1_43); double** m1_44; createMatrix(mrows, mcols, m1_44); double** m1_45; createMatrix(mrows, mcols, m1_45);
	double** m1_51; createMatrix(mrows, mcols, m1_51); double** m1_52; createMatrix(mrows, mcols, m1_52); double** m1_53; createMatrix(mrows, mcols, m1_53); double** m1_54; createMatrix(mrows, mcols, m1_54); double** m1_55; createMatrix(mrows, mcols, m1_55);

	double** m2_11; createMatrix(mrows, mcols, m2_11); double** m2_12; createMatrix(mrows, mcols, m2_12); double** m2_13; createMatrix(mrows, mcols, m2_13); double** m2_14; createMatrix(mrows, mcols, m2_14); double** m2_15; createMatrix(mrows, mcols, m2_15);
	double** m2_21; createMatrix(mrows, mcols, m2_21); double** m2_22; createMatrix(mrows, mcols, m2_22); double** m2_23; createMatrix(mrows, mcols, m2_23); double** m2_24; createMatrix(mrows, mcols, m2_24); double** m2_25; createMatrix(mrows, mcols, m2_25);
	double** m2_31; createMatrix(mrows, mcols, m2_31); double** m2_32; createMatrix(mrows, mcols, m2_32); double** m2_33; createMatrix(mrows, mcols, m2_33); double** m2_34; createMatrix(mrows, mcols, m2_34); double** m2_35; createMatrix(mrows, mcols, m2_35);
	double** m2_41; createMatrix(mrows, mcols, m2_41); double** m2_42; createMatrix(mrows, mcols, m2_42); double** m2_43; createMatrix(mrows, mcols, m2_43); double** m2_44; createMatrix(mrows, mcols, m2_44); double** m2_45; createMatrix(mrows, mcols, m2_45);
	double** m2_51; createMatrix(mrows, mcols, m2_51); double** m2_52; createMatrix(mrows, mcols, m2_52); double** m2_53; createMatrix(mrows, mcols, m2_53); double** m2_54; createMatrix(mrows, mcols, m2_54); double** m2_55; createMatrix(mrows, mcols, m2_55);

	double** Ee_11; createMatrix(mrows, mcols, Ee_11); double** Ee_12; createMatrix(mrows, mcols, Ee_12); double** Ee_13; createMatrix(mrows, mcols, Ee_13); double** Ee_14; createMatrix(mrows, mcols, Ee_14); double** Ee_15; createMatrix(mrows, mcols, Ee_15);
	double** Ee_21; createMatrix(mrows, mcols, Ee_21); double** Ee_22; createMatrix(mrows, mcols, Ee_22); double** Ee_23; createMatrix(mrows, mcols, Ee_23); double** Ee_24; createMatrix(mrows, mcols, Ee_24); double** Ee_25; createMatrix(mrows, mcols, Ee_25);
	double** Ee_31; createMatrix(mrows, mcols, Ee_31); double** Ee_32; createMatrix(mrows, mcols, Ee_32); double** Ee_33; createMatrix(mrows, mcols, Ee_33); double** Ee_34; createMatrix(mrows, mcols, Ee_34); double** Ee_35; createMatrix(mrows, mcols, Ee_35);
	double** Ee_41; createMatrix(mrows, mcols, Ee_41); double** Ee_42; createMatrix(mrows, mcols, Ee_42); double** Ee_43; createMatrix(mrows, mcols, Ee_43); double** Ee_44; createMatrix(mrows, mcols, Ee_44); double** Ee_45; createMatrix(mrows, mcols, Ee_45);
	double** Ee_51; createMatrix(mrows, mcols, Ee_51); double** Ee_52; createMatrix(mrows, mcols, Ee_52); double** Ee_53; createMatrix(mrows, mcols, Ee_53); double** Ee_54; createMatrix(mrows, mcols, Ee_54); double** Ee_55; createMatrix(mrows, mcols, Ee_55);

	double** Ei_11; createMatrix(mrows, mcols, Ei_11); double** Ei_12; createMatrix(mrows, mcols, Ei_12); double** Ei_13; createMatrix(mrows, mcols, Ei_13); double** Ei_14; createMatrix(mrows, mcols, Ei_14); double** Ei_15; createMatrix(mrows, mcols, Ei_15);
	double** Ei_21; createMatrix(mrows, mcols, Ei_21); double** Ei_22; createMatrix(mrows, mcols, Ei_22); double** Ei_23; createMatrix(mrows, mcols, Ei_23); double** Ei_24; createMatrix(mrows, mcols, Ei_24); double** Ei_25; createMatrix(mrows, mcols, Ei_25);
	double** Ei_31; createMatrix(mrows, mcols, Ei_31); double** Ei_32; createMatrix(mrows, mcols, Ei_32); double** Ei_33; createMatrix(mrows, mcols, Ei_33); double** Ei_34; createMatrix(mrows, mcols, Ei_34); double** Ei_35; createMatrix(mrows, mcols, Ei_35);
	double** Ei_41; createMatrix(mrows, mcols, Ei_41); double** Ei_42; createMatrix(mrows, mcols, Ei_42); double** Ei_43; createMatrix(mrows, mcols, Ei_43); double** Ei_44; createMatrix(mrows, mcols, Ei_44); double** Ei_45; createMatrix(mrows, mcols, Ei_45);
	double** Ei_51; createMatrix(mrows, mcols, Ei_51); double** Ei_52; createMatrix(mrows, mcols, Ei_52); double** Ei_53; createMatrix(mrows, mcols, Ei_53); double** Ei_54; createMatrix(mrows, mcols, Ei_54); double** Ei_55; createMatrix(mrows, mcols, Ei_55);

	double** Er_11; createMatrix(mrows, mcols, Er_11); double** Er_12; createMatrix(mrows, mcols, Er_12); double** Er_13; createMatrix(mrows, mcols, Er_13); double** Er_14; createMatrix(mrows, mcols, Er_14); double** Er_15; createMatrix(mrows, mcols, Er_15);
	double** Er_21; createMatrix(mrows, mcols, Er_21); double** Er_22; createMatrix(mrows, mcols, Er_22); double** Er_23; createMatrix(mrows, mcols, Er_23); double** Er_24; createMatrix(mrows, mcols, Er_24); double** Er_25; createMatrix(mrows, mcols, Er_25);
	double** Er_31; createMatrix(mrows, mcols, Er_31); double** Er_32; createMatrix(mrows, mcols, Er_32); double** Er_33; createMatrix(mrows, mcols, Er_33); double** Er_34; createMatrix(mrows, mcols, Er_34); double** Er_35; createMatrix(mrows, mcols, Er_35);
	double** Er_41; createMatrix(mrows, mcols, Er_41); double** Er_42; createMatrix(mrows, mcols, Er_42); double** Er_43; createMatrix(mrows, mcols, Er_43); double** Er_44; createMatrix(mrows, mcols, Er_44); double** Er_45; createMatrix(mrows, mcols, Er_45);
	double** Er_51; createMatrix(mrows, mcols, Er_51); double** Er_52; createMatrix(mrows, mcols, Er_52); double** Er_53; createMatrix(mrows, mcols, Er_53); double** Er_54; createMatrix(mrows, mcols, Er_54); double** Er_55; createMatrix(mrows, mcols, Er_55);
	
	for (int i = 0; i < N_y + 2; i++)
	{
		for (int j = 0; j < N_x + 2; j++)
		{
			rho_11[i][j] = rho_ex(X_1[j], Y_1[i], t); rho_21[i][j] = rho_ex(X_2[j], Y_1[i], t); rho_31[i][j] = rho_ex(X_3[j], Y_1[i], t);
			rho_41[i][j] = rho_ex(X_4[j], Y_1[i], t); rho_51[i][j] = rho_ex(X_5[j], Y_1[i], t);

			rho_12[i][j] = rho_ex(X_1[j], Y_2[i], t); rho_22[i][j] = rho_ex(X_2[j], Y_2[i], t); rho_32[i][j] = rho_ex(X_3[j], Y_2[i], t);
			rho_42[i][j] = rho_ex(X_4[j], Y_2[i], t); rho_52[i][j] = rho_ex(X_5[j], Y_2[i], t);

			rho_13[i][j] = rho_ex(X_1[j], Y_3[i], t); rho_23[i][j] = rho_ex(X_2[j], Y_3[i], t); rho_33[i][j] = rho_ex(X_3[j], Y_3[i], t);
			rho_43[i][j] = rho_ex(X_4[j], Y_3[i], t); rho_53[i][j] = rho_ex(X_5[j], Y_3[i], t);

			rho_14[i][j] = rho_ex(X_1[j], Y_4[i], t); rho_24[i][j] = rho_ex(X_2[j], Y_4[i], t); rho_34[i][j] = rho_ex(X_3[j], Y_4[i], t);
			rho_44[i][j] = rho_ex(X_4[j], Y_4[i], t); rho_54[i][j] = rho_ex(X_5[j], Y_4[i], t);

			rho_15[i][j] = rho_ex(X_1[j], Y_5[i], t); rho_25[i][j] = rho_ex(X_2[j], Y_5[i], t); rho_35[i][j] = rho_ex(X_3[j], Y_5[i], t);
			rho_45[i][j] = rho_ex(X_4[j], Y_5[i], t); rho_55[i][j] = rho_ex(X_5[j], Y_5[i], t);

			m1_11[i][j] = m1_ex(X_1[j], Y_1[i], t); m1_21[i][j] = m1_ex(X_2[j], Y_1[i], t); m1_31[i][j] = m1_ex(X_3[j], Y_1[i], t);
			m1_41[i][j] = m1_ex(X_4[j], Y_1[i], t); m1_51[i][j] = m1_ex(X_5[j], Y_1[i], t);

			m1_12[i][j] = m1_ex(X_1[j], Y_2[i], t); m1_22[i][j] = m1_ex(X_2[j], Y_2[i], t); m1_32[i][j] = m1_ex(X_3[j], Y_2[i], t);
			m1_42[i][j] = m1_ex(X_4[j], Y_2[i], t); m1_52[i][j] = m1_ex(X_5[j], Y_2[i], t);

			m1_13[i][j] = m1_ex(X_1[j], Y_3[i], t); m1_23[i][j] = m1_ex(X_2[j], Y_3[i], t); m1_33[i][j] = m1_ex(X_3[j], Y_3[i], t);
			m1_43[i][j] = m1_ex(X_4[j], Y_3[i], t); m1_53[i][j] = m1_ex(X_5[j], Y_3[i], t);

			m1_14[i][j] = m1_ex(X_1[j], Y_4[i], t); m1_24[i][j] = m1_ex(X_2[j], Y_4[i], t); m1_34[i][j] = m1_ex(X_3[j], Y_4[i], t);
			m1_44[i][j] = m1_ex(X_4[j], Y_4[i], t); m1_54[i][j] = m1_ex(X_5[j], Y_4[i], t);

			m1_15[i][j] = m1_ex(X_1[j], Y_5[i], t); m1_25[i][j] = m1_ex(X_2[j], Y_5[i], t); m1_35[i][j] = m1_ex(X_3[j], Y_5[i], t);
			m1_45[i][j] = m1_ex(X_4[j], Y_5[i], t); m1_55[i][j] = m1_ex(X_5[j], Y_5[i], t);

			m2_11[i][j] = m2_ex(X_1[j], Y_1[i], t); m2_21[i][j] = m2_ex(X_2[j], Y_1[i], t); m2_31[i][j] = m2_ex(X_3[j], Y_1[i], t);
			m2_41[i][j] = m2_ex(X_4[j], Y_1[i], t); m2_51[i][j] = m2_ex(X_5[j], Y_1[i], t);

			m2_12[i][j] = m2_ex(X_1[j], Y_2[i], t); m2_22[i][j] = m2_ex(X_2[j], Y_2[i], t); m2_32[i][j] = m2_ex(X_3[j], Y_2[i], t);
			m2_42[i][j] = m2_ex(X_4[j], Y_2[i], t); m2_52[i][j] = m2_ex(X_5[j], Y_2[i], t);

			m2_13[i][j] = m2_ex(X_1[j], Y_3[i], t); m2_23[i][j] = m2_ex(X_2[j], Y_3[i], t); m2_33[i][j] = m2_ex(X_3[j], Y_3[i], t);
			m2_43[i][j] = m2_ex(X_4[j], Y_3[i], t); m2_53[i][j] = m2_ex(X_5[j], Y_3[i], t);

			m2_14[i][j] = m2_ex(X_1[j], Y_4[i], t); m2_24[i][j] = m2_ex(X_2[j], Y_4[i], t); m2_34[i][j] = m2_ex(X_3[j], Y_4[i], t);
			m2_44[i][j] = m2_ex(X_4[j], Y_4[i], t); m2_54[i][j] = m2_ex(X_5[j], Y_4[i], t);

			m2_15[i][j] = m2_ex(X_1[j], Y_5[i], t); m2_25[i][j] = m2_ex(X_2[j], Y_5[i], t); m2_35[i][j] = m2_ex(X_3[j], Y_5[i], t);
			m2_45[i][j] = m2_ex(X_4[j], Y_5[i], t); m2_55[i][j] = m2_ex(X_5[j], Y_5[i], t);

			Ee_11[i][j] = Ee_ex(X_1[j], Y_1[i], t); Ee_21[i][j] = Ee_ex(X_2[j], Y_1[i], t); Ee_31[i][j] = Ee_ex(X_3[j], Y_1[i], t);
			Ee_41[i][j] = Ee_ex(X_4[j], Y_1[i], t); Ee_51[i][j] = Ee_ex(X_5[j], Y_1[i], t);

			Ee_12[i][j] = Ee_ex(X_1[j], Y_2[i], t); Ee_22[i][j] = Ee_ex(X_2[j], Y_2[i], t); Ee_32[i][j] = Ee_ex(X_3[j], Y_2[i], t);
			Ee_42[i][j] = Ee_ex(X_4[j], Y_2[i], t); Ee_52[i][j] = Ee_ex(X_5[j], Y_2[i], t);

			Ee_13[i][j] = Ee_ex(X_1[j], Y_3[i], t); Ee_23[i][j] = Ee_ex(X_2[j], Y_3[i], t); Ee_33[i][j] = Ee_ex(X_3[j], Y_3[i], t);
			Ee_43[i][j] = Ee_ex(X_4[j], Y_3[i], t); Ee_53[i][j] = Ee_ex(X_5[j], Y_3[i], t);

			Ee_14[i][j] = Ee_ex(X_1[j], Y_4[i], t); Ee_24[i][j] = Ee_ex(X_2[j], Y_4[i], t); Ee_34[i][j] = Ee_ex(X_3[j], Y_4[i], t);
			Ee_44[i][j] = Ee_ex(X_4[j], Y_4[i], t); Ee_54[i][j] = Ee_ex(X_5[j], Y_4[i], t);

			Ee_15[i][j] = Ee_ex(X_1[j], Y_5[i], t); Ee_25[i][j] = Ee_ex(X_2[j], Y_5[i], t); Ee_35[i][j] = Ee_ex(X_3[j], Y_5[i], t);
			Ee_45[i][j] = Ee_ex(X_4[j], Y_5[i], t); Ee_55[i][j] = Ee_ex(X_5[j], Y_5[i], t);

			Ei_11[i][j] = Ei_ex(X_1[j], Y_1[i], t); Ei_21[i][j] = Ei_ex(X_2[j], Y_1[i], t); Ei_31[i][j] = Ei_ex(X_3[j], Y_1[i], t);
			Ei_41[i][j] = Ei_ex(X_4[j], Y_1[i], t); Ei_51[i][j] = Ei_ex(X_5[j], Y_1[i], t);

			Ei_12[i][j] = Ei_ex(X_1[j], Y_2[i], t); Ei_22[i][j] = Ei_ex(X_2[j], Y_2[i], t); Ei_32[i][j] = Ei_ex(X_3[j], Y_2[i], t);
			Ei_42[i][j] = Ei_ex(X_4[j], Y_2[i], t); Ei_52[i][j] = Ei_ex(X_5[j], Y_2[i], t);

			Ei_13[i][j] = Ei_ex(X_1[j], Y_3[i], t); Ei_23[i][j] = Ei_ex(X_2[j], Y_3[i], t); Ei_33[i][j] = Ei_ex(X_3[j], Y_3[i], t);
			Ei_43[i][j] = Ei_ex(X_4[j], Y_3[i], t); Ei_53[i][j] = Ei_ex(X_5[j], Y_3[i], t);

			Ei_14[i][j] = Ei_ex(X_1[j], Y_4[i], t); Ei_24[i][j] = Ei_ex(X_2[j], Y_4[i], t); Ei_34[i][j] = Ei_ex(X_3[j], Y_4[i], t);
			Ei_44[i][j] = Ei_ex(X_4[j], Y_4[i], t); Ei_54[i][j] = Ei_ex(X_5[j], Y_4[i], t);

			Ei_15[i][j] = Ei_ex(X_1[j], Y_5[i], t); Ei_25[i][j] = Ei_ex(X_2[j], Y_5[i], t); Ei_35[i][j] = Ei_ex(X_3[j], Y_5[i], t);
			Ei_45[i][j] = Ei_ex(X_4[j], Y_5[i], t); Ei_55[i][j] = Ei_ex(X_5[j], Y_5[i], t);

			Er_11[i][j] = Er_ex(X_1[j], Y_1[i], t); Er_21[i][j] = Er_ex(X_2[j], Y_1[i], t); Er_31[i][j] = Er_ex(X_3[j], Y_1[i], t);
			Er_41[i][j] = Er_ex(X_4[j], Y_1[i], t); Er_51[i][j] = Er_ex(X_5[j], Y_1[i], t);

			Er_12[i][j] = Er_ex(X_1[j], Y_2[i], t); Er_22[i][j] = Er_ex(X_2[j], Y_2[i], t); Er_32[i][j] = Er_ex(X_3[j], Y_2[i], t);
			Er_42[i][j] = Er_ex(X_4[j], Y_2[i], t); Er_52[i][j] = Er_ex(X_5[j], Y_2[i], t);

			Er_13[i][j] = Er_ex(X_1[j], Y_3[i], t); Er_23[i][j] = Er_ex(X_2[j], Y_3[i], t); Er_33[i][j] = Er_ex(X_3[j], Y_3[i], t);
			Er_43[i][j] = Er_ex(X_4[j], Y_3[i], t); Er_53[i][j] = Er_ex(X_5[j], Y_3[i], t);

			Er_14[i][j] = Er_ex(X_1[j], Y_4[i], t); Er_24[i][j] = Er_ex(X_2[j], Y_4[i], t); Er_34[i][j] = Er_ex(X_3[j], Y_4[i], t);
			Er_44[i][j] = Er_ex(X_4[j], Y_4[i], t); Er_54[i][j] = Er_ex(X_5[j], Y_4[i], t);

			Er_15[i][j] = Er_ex(X_1[j], Y_5[i], t); Er_25[i][j] = Er_ex(X_2[j], Y_5[i], t); Er_35[i][j] = Er_ex(X_3[j], Y_5[i], t);
			Er_45[i][j] = Er_ex(X_4[j], Y_5[i], t); Er_55[i][j] = Er_ex(X_5[j], Y_5[i], t);
		}
	}

	double** GQ_0; createMatrix(mrows, mcols, GQ_0); double** GQ_1; createMatrix(mrows, mcols, GQ_1); double** GQ_2; createMatrix(mrows, mcols, GQ_2);
	double** GQ_3; createMatrix(mrows, mcols, GQ_3); double** GQ_4; createMatrix(mrows, mcols, GQ_4); double** GQ_5; createMatrix(mrows, mcols, GQ_5);
	double** rho_0; createMatrix(mrows, mcols, rho_0); double** rho_1; createMatrix(mrows, mcols, rho_1); double** rho_2; createMatrix(mrows, mcols, rho_2);
	double** rho_3; createMatrix(mrows, mcols, rho_3); double** rho_4; createMatrix(mrows, mcols, rho_4); double** rho_5; createMatrix(mrows, mcols, rho_5);
	double** m1_0; createMatrix(mrows, mcols, m1_0); double** m1_1; createMatrix(mrows, mcols, m1_1); double** m1_2; createMatrix(mrows, mcols, m1_2);
	double** m1_3; createMatrix(mrows, mcols, m1_3); double** m1_4; createMatrix(mrows, mcols, m1_4); double** m1_5; createMatrix(mrows, mcols, m1_5);
	double** m2_0; createMatrix(mrows, mcols, m2_0); double** m2_1; createMatrix(mrows, mcols, m2_1); double** m2_2; createMatrix(mrows, mcols, m2_2);
	double** m2_3; createMatrix(mrows, mcols, m2_3); double** m2_4; createMatrix(mrows, mcols, m2_4); double** m2_5; createMatrix(mrows, mcols, m2_5);
	double** Ee_0; createMatrix(mrows, mcols, Ee_0); double** Ee_1; createMatrix(mrows, mcols, Ee_1); double** Ee_2; createMatrix(mrows, mcols, Ee_2);
	double** Ee_3; createMatrix(mrows, mcols, Ee_3); double** Ee_4; createMatrix(mrows, mcols, Ee_4); double** Ee_5; createMatrix(mrows, mcols, Ee_5);
	double** Ei_0; createMatrix(mrows, mcols, Ei_0); double** Ei_1; createMatrix(mrows, mcols, Ei_1); double** Ei_2; createMatrix(mrows, mcols, Ei_2);
	double** Ei_3; createMatrix(mrows, mcols, Ei_3); double** Ei_4; createMatrix(mrows, mcols, Ei_4); double** Ei_5; createMatrix(mrows, mcols, Ei_5);
	double** Er_0; createMatrix(mrows, mcols, Er_0); double** Er_1; createMatrix(mrows, mcols, Er_1); double** Er_2; createMatrix(mrows, mcols, Er_2);
	double** Er_3; createMatrix(mrows, mcols, Er_3); double** Er_4; createMatrix(mrows, mcols, Er_4); double** Er_5; createMatrix(mrows, mcols, Er_5);

	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, rho_11, rho_21, rho_31, rho_41, rho_51, rho_12, rho_22, rho_32, rho_42, rho_52, rho_13, rho_23, rho_33, rho_43, rho_53, rho_14, rho_24, rho_34, rho_44, rho_54, rho_15, rho_25, rho_35, rho_45, rho_55, h_x, h_y);
	get_L2_projection_2D(rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);

	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, m1_11, m1_21, m1_31, m1_41, m1_51, m1_12, m1_22, m1_32, m1_42, m1_52, m1_13, m1_23, m1_33, m1_43, m1_53, m1_14, m1_24, m1_34, m1_44, m1_54, m1_15, m1_25, m1_35, m1_45, m1_55, h_x, h_y);
	get_L2_projection_2D(m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);

	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, m2_11, m2_21, m2_31, m2_41, m2_51, m2_12, m2_22, m2_32, m2_42, m2_52, m2_13, m2_23, m2_33, m2_43, m2_53, m2_14, m2_24, m2_34, m2_44, m2_54, m2_15, m2_25, m2_35, m2_45, m2_55, h_x, h_y);
	get_L2_projection_2D(m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);
	
	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, Ee_11, Ee_21, Ee_31, Ee_41, Ee_51, Ee_12, Ee_22, Ee_32, Ee_42, Ee_52, Ee_13, Ee_23, Ee_33, Ee_43, Ee_53, Ee_14, Ee_24, Ee_34, Ee_44, Ee_54, Ee_15, Ee_25, Ee_35, Ee_45, Ee_55, h_x, h_y);
	get_L2_projection_2D(Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);

	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, Ei_11, Ei_21, Ei_31, Ei_41, Ei_51, Ei_12, Ei_22, Ei_32, Ei_42, Ei_52, Ei_13, Ei_23, Ei_33, Ei_43, Ei_53, Ei_14, Ei_24, Ei_34, Ei_44, Ei_54, Ei_15, Ei_25, Ei_35, Ei_45, Ei_55, h_x, h_y);
	get_L2_projection_2D(Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);

	get_Gauss_quadrature_2D(GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, Er_11, Er_21, Er_31, Er_41, Er_51, Er_12, Er_22, Er_32, Er_42, Er_52, Er_13, Er_23, Er_33, Er_43, Er_53, Er_14, Er_24, Er_34, Er_44, Er_54, Er_15, Er_25, Er_35, Er_45, Er_55, h_x, h_y);
	get_L2_projection_2D(Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, GQ_0, GQ_1, GQ_2, GQ_3, GQ_4, GQ_5, h_x, h_y);

	deleteMatrix(mrows, rho_11); deleteMatrix(mrows, rho_12); deleteMatrix(mrows, rho_13); deleteMatrix(mrows, rho_14); deleteMatrix(mrows, rho_15);
	deleteMatrix(mrows, rho_21); deleteMatrix(mrows, rho_22); deleteMatrix(mrows, rho_23); deleteMatrix(mrows, rho_24); deleteMatrix(mrows, rho_25);
	deleteMatrix(mrows, rho_31); deleteMatrix(mrows, rho_32); deleteMatrix(mrows, rho_33); deleteMatrix(mrows, rho_34); deleteMatrix(mrows, rho_35);
	deleteMatrix(mrows, rho_41); deleteMatrix(mrows, rho_42); deleteMatrix(mrows, rho_43); deleteMatrix(mrows, rho_44); deleteMatrix(mrows, rho_45);
	deleteMatrix(mrows, rho_51); deleteMatrix(mrows, rho_52); deleteMatrix(mrows, rho_53); deleteMatrix(mrows, rho_54); deleteMatrix(mrows, rho_55);

	deleteMatrix(mrows, m1_11); deleteMatrix(mrows, m1_12); deleteMatrix(mrows, m1_13); deleteMatrix(mrows, m1_14); deleteMatrix(mrows, m1_15);
	deleteMatrix(mrows, m1_21); deleteMatrix(mrows, m1_22); deleteMatrix(mrows, m1_23); deleteMatrix(mrows, m1_24); deleteMatrix(mrows, m1_25);
	deleteMatrix(mrows, m1_31); deleteMatrix(mrows, m1_32); deleteMatrix(mrows, m1_33); deleteMatrix(mrows, m1_34); deleteMatrix(mrows, m1_35);
	deleteMatrix(mrows, m1_41); deleteMatrix(mrows, m1_42); deleteMatrix(mrows, m1_43); deleteMatrix(mrows, m1_44); deleteMatrix(mrows, m1_45);
	deleteMatrix(mrows, m1_51); deleteMatrix(mrows, m1_52); deleteMatrix(mrows, m1_53); deleteMatrix(mrows, m1_54); deleteMatrix(mrows, m1_55);

	deleteMatrix(mrows, m2_11); deleteMatrix(mrows, m2_12); deleteMatrix(mrows, m2_13); deleteMatrix(mrows, m2_14); deleteMatrix(mrows, m2_15);
	deleteMatrix(mrows, m2_21); deleteMatrix(mrows, m2_22); deleteMatrix(mrows, m2_23); deleteMatrix(mrows, m2_24); deleteMatrix(mrows, m2_25);
	deleteMatrix(mrows, m2_31); deleteMatrix(mrows, m2_32); deleteMatrix(mrows, m2_33); deleteMatrix(mrows, m2_34); deleteMatrix(mrows, m2_35);
	deleteMatrix(mrows, m2_41); deleteMatrix(mrows, m2_42); deleteMatrix(mrows, m2_43); deleteMatrix(mrows, m2_44); deleteMatrix(mrows, m2_45);
	deleteMatrix(mrows, m2_51); deleteMatrix(mrows, m2_52); deleteMatrix(mrows, m2_53); deleteMatrix(mrows, m2_54); deleteMatrix(mrows, m2_55);

	deleteMatrix(mrows, Ee_11); deleteMatrix(mrows, Ee_12); deleteMatrix(mrows, Ee_13); deleteMatrix(mrows, Ee_14); deleteMatrix(mrows, Ee_15);
	deleteMatrix(mrows, Ee_21); deleteMatrix(mrows, Ee_22); deleteMatrix(mrows, Ee_23); deleteMatrix(mrows, Ee_24); deleteMatrix(mrows, Ee_25);
	deleteMatrix(mrows, Ee_31); deleteMatrix(mrows, Ee_32); deleteMatrix(mrows, Ee_33); deleteMatrix(mrows, Ee_34); deleteMatrix(mrows, Ee_35);
	deleteMatrix(mrows, Ee_41); deleteMatrix(mrows, Ee_42); deleteMatrix(mrows, Ee_43); deleteMatrix(mrows, Ee_44); deleteMatrix(mrows, Ee_45);
	deleteMatrix(mrows, Ee_51); deleteMatrix(mrows, Ee_52); deleteMatrix(mrows, Ee_53); deleteMatrix(mrows, Ee_54); deleteMatrix(mrows, Ee_55);

	deleteMatrix(mrows, Ei_11); deleteMatrix(mrows, Ei_12); deleteMatrix(mrows, Ei_13); deleteMatrix(mrows, Ei_14); deleteMatrix(mrows, Ei_15);
	deleteMatrix(mrows, Ei_21); deleteMatrix(mrows, Ei_22); deleteMatrix(mrows, Ei_23); deleteMatrix(mrows, Ei_24); deleteMatrix(mrows, Ei_25);
	deleteMatrix(mrows, Ei_31); deleteMatrix(mrows, Ei_32); deleteMatrix(mrows, Ei_33); deleteMatrix(mrows, Ei_34); deleteMatrix(mrows, Ei_35);
	deleteMatrix(mrows, Ei_41); deleteMatrix(mrows, Ei_42); deleteMatrix(mrows, Ei_43); deleteMatrix(mrows, Ei_44); deleteMatrix(mrows, Ei_45);
	deleteMatrix(mrows, Ei_51); deleteMatrix(mrows, Ei_52); deleteMatrix(mrows, Ei_53); deleteMatrix(mrows, Ei_54); deleteMatrix(mrows, Ei_55);

	deleteMatrix(mrows, Er_11); deleteMatrix(mrows, Er_12); deleteMatrix(mrows, Er_13); deleteMatrix(mrows, Er_14); deleteMatrix(mrows, Er_15);
	deleteMatrix(mrows, Er_21); deleteMatrix(mrows, Er_22); deleteMatrix(mrows, Er_23); deleteMatrix(mrows, Er_24); deleteMatrix(mrows, Er_25);
	deleteMatrix(mrows, Er_31); deleteMatrix(mrows, Er_32); deleteMatrix(mrows, Er_33); deleteMatrix(mrows, Er_34); deleteMatrix(mrows, Er_35);
	deleteMatrix(mrows, Er_41); deleteMatrix(mrows, Er_42); deleteMatrix(mrows, Er_43); deleteMatrix(mrows, Er_44); deleteMatrix(mrows, Er_45);
	deleteMatrix(mrows, Er_51); deleteMatrix(mrows, Er_52); deleteMatrix(mrows, Er_53); deleteMatrix(mrows, Er_54); deleteMatrix(mrows, Er_55);

	deleteMatrix(mrows, GQ_0); deleteMatrix(mrows, GQ_1); deleteMatrix(mrows, GQ_2);
	deleteMatrix(mrows, GQ_3); deleteMatrix(mrows, GQ_4); deleteMatrix(mrows, GQ_5);

	// determining time step
	double** s_j; createMatrix(N_y, N_x, s_j);
	compute_st(s_j, rho_0, m1_0, m2_0, Ee_0);
	double s_max = compute_sj_max(s_j);
	double tau_s = CFL / s_max;

	double** beta_ij_x; createMatrix(N_y, N_x, beta_ij_x); double** beta_ij_y; createMatrix(N_y, N_x, beta_ij_y);
	OE_compute_beta_ij(beta_ij_x, beta_ij_y, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
    double betax_max = compute_sj_max(beta_ij_x); double betay_max = compute_sj_max(beta_ij_y);
	double tau_a = CFL / (betax_max / h_x + betay_max / h_y);

	double tau = tau_a;
	if (tau_s < tau)
	{
		tau = tau_s;
	}
	double t_fin = T_final;

	cout << tau_a << endl;
	cout << tau_s << endl;
	cout << tau << endl;

	// IMEX3 
	double b1 = -3.0 / 2.0 * pow(gamma, 2) + 4.0 * gamma - 1.0 / 4.0;
	double b2 = 3.0 / 2.0 * pow(gamma, 2) - 5.0 * gamma + 5.0 / 4.0;
	double a1 = -0.35;
	double a2 = (1.0 / 3.0 - 2.0 * pow(gamma, 2) - 2.0 * b2 * a1 * gamma) / (gamma * (1.0 - gamma));

	double rho_min = comput_rho_min(rho_0);
	double a0 = compute_dj(rho_min);

	// mass matrix
	double alpha_0 = Coef_d2; double beta_0 = diff_beta;

	int m_rows = (mrows - 2) * (mcols - 2); int m_cols = m_rows;
	int M_rows = 6 * m_rows; int M_cols = M_rows;
	
	Eigen::SparseMatrix<double> M(M_rows, M_cols);
	piece_Linear_coef_matrix(M, M_rows, M_cols, m_rows, m_cols, alpha_0, beta_0, h_x, h_y, 0.0, 1.0);
	
	double coff = -1.0 * tau * gamma * a0;
	Eigen::SparseMatrix<double> A(M_rows, M_cols);
	piece_Linear_coef_matrix(A, M_rows, M_cols, m_rows, m_cols, alpha_0, beta_0, h_x, h_y, 1.0, coff);

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	double** E; createMatrix(mrows, mcols, E); double** E_new; createMatrix(mrows, mcols, E_new);
	double rho_c = 0.0; double m1_c = 0.0; double m2_c = 0.0; double E_c = 0.0;
	// time stepping
	cout << "Method: Modal DG with P" << Mode << " element. Mesh: " << N_x << "¡Á" << N_y << endl;
	cout << "Start computing the numerical solution ... " << endl;
	int Timestep_index = 0; 
	double Tau = 10.0;

	while (t < t_fin)
	{
		if (t_fin - t < tau)
		{
			tau = t_fin - t;
			coff = -1.0 * tau * gamma * a0;
			piece_Linear_coef_matrix(A, M_rows, M_cols, m_rows, m_cols, alpha_0, beta_0, h_x, h_y, 1.0, coff);
			solver.compute(A);
		}

		Timestep_index = Timestep_index + 1;
		
		if (fmod(Timestep_index, N_timestep) == 0)
		{
			compute_st(s_j, rho_0, m1_0, m2_0, Ee_0);
			s_max = compute_sj_max(s_j);
			tau_s = CFL / s_max;

			OE_compute_beta_ij(beta_ij_x, beta_ij_y, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
			betax_max = compute_sj_max(beta_ij_x); double betay_max = compute_sj_max(beta_ij_y);
			tau_a = CFL / (betax_max / h_x + betay_max / h_y);

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
			//cout << Timestep_index << " " << Tau << " " << tau << endl;
			coff = -1.0 * tau * gamma * a0;
			piece_Linear_coef_matrix(A, M_rows, M_cols, m_rows, m_cols, alpha_0, beta_0, h_x, h_y, 1.0, coff);
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

		for (int i = 0; i < N_y + 2; i++)
		{
			for (int j = 0; j < N_x + 2; j++)
			{
				E[i][j] = Ee_0[i][j] + Ei_0[i][j] + Er_0[i][j];
			}
		}

		double** sigma_0; createMatrix(N_y, N_x, sigma_0); double** sigma_1; createMatrix(N_y, N_x, sigma_1);
		double** sigma_2; createMatrix(N_y, N_x, sigma_2);

		double** coef_L1; createMatrix(N_y, N_x, coef_L1); double** coef_L2; createMatrix(N_y, N_x, coef_L2);
		double** coef_L3; createMatrix(N_y, N_x, coef_L3); double** coef_L4; createMatrix(N_y, N_x, coef_L4);
		double** coef_L5; createMatrix(N_y, N_x, coef_L5);

		double** coef_R1; createMatrix(N_y, N_x, coef_R1); double** coef_R2; createMatrix(N_y, N_x, coef_R2);
		double** coef_R3; createMatrix(N_y, N_x, coef_R3); double** coef_R4; createMatrix(N_y, N_x, coef_R4);
		double** coef_R5; createMatrix(N_y, N_x, coef_R5);

		double** coef_B1; createMatrix(N_y, N_x, coef_B1); double** coef_B2; createMatrix(N_y, N_x, coef_B2);
		double** coef_B3; createMatrix(N_y, N_x, coef_B3); double** coef_B4; createMatrix(N_y, N_x, coef_B4);
		double** coef_B5; createMatrix(N_y, N_x, coef_B5);

		double** coef_T1; createMatrix(N_y, N_x, coef_T1); double** coef_T2; createMatrix(N_y, N_x, coef_T2);
		double** coef_T3; createMatrix(N_y, N_x, coef_T3); double** coef_T4; createMatrix(N_y, N_x, coef_T4);
		double** coef_T5; createMatrix(N_y, N_x, coef_T5);

		double** F1_alphac1; createMatrix(N_y, N_x + 1, F1_alphac1); double** F1_alphac2; createMatrix(N_y, N_x + 1, F1_alphac2);
		double** F1_alphac3; createMatrix(N_y, N_x + 1, F1_alphac3); double** F1_alphac4; createMatrix(N_y, N_x + 1, F1_alphac4);
		double** F1_alphac5; createMatrix(N_y, N_x + 1, F1_alphac5); double** F1_alphac; createMatrix(N_y, N_x + 1, F1_alphac);

		double** F2_alphac1; createMatrix(N_y + 1, N_x, F2_alphac1); double** F2_alphac2; createMatrix(N_y + 1, N_x, F2_alphac2);
		double** F2_alphac3; createMatrix(N_y + 1, N_x, F2_alphac3); double** F2_alphac4; createMatrix(N_y + 1, N_x, F2_alphac4);
		double** F2_alphac5; createMatrix(N_y + 1, N_x, F2_alphac5); double** F2_alphac; createMatrix(N_y + 1, N_x, F2_alphac);

		double** S_0; createMatrix(mrows, mcols, S_0); double** S_1; createMatrix(mrows, mcols, S_1);
		double** S_2; createMatrix(mrows, mcols, S_2); double** S_3; createMatrix(mrows, mcols, S_3);
		double** S_4; createMatrix(mrows, mcols, S_4); double** S_5; createMatrix(mrows, mcols, S_5);

		double t_1 = t; double t_2 = t + gamma * tau;
		double t_3 = t + (1.0 + gamma) / 2.0 * tau; double t_4 = t + tau;

		// First step of IMEX3
		double** rhof_0; createMatrix(mrows, mcols, rhof_0); double** rhof_1; createMatrix(mrows, mcols, rhof_1); double** rhof_2; createMatrix(mrows, mcols, rhof_2);
		double** rhof_3; createMatrix(mrows, mcols, rhof_3); double** rhof_4; createMatrix(mrows, mcols, rhof_4); double** rhof_5; createMatrix(mrows, mcols, rhof_5);
		double** m1f_0; createMatrix(mrows, mcols, m1f_0); double** m1f_1; createMatrix(mrows, mcols, m1f_1); double** m1f_2; createMatrix(mrows, mcols, m1f_2);
		double** m1f_3; createMatrix(mrows, mcols, m1f_3); double** m1f_4; createMatrix(mrows, mcols, m1f_4); double** m1f_5; createMatrix(mrows, mcols, m1f_5);
		double** m2f_0; createMatrix(mrows, mcols, m2f_0); double** m2f_1; createMatrix(mrows, mcols, m2f_1); double** m2f_2; createMatrix(mrows, mcols, m2f_2);
		double** m2f_3; createMatrix(mrows, mcols, m2f_3); double** m2f_4; createMatrix(mrows, mcols, m2f_4); double** m2f_5; createMatrix(mrows, mcols, m2f_5);
		double** Eef_0; createMatrix(mrows, mcols, Eef_0); double** Eef_1; createMatrix(mrows, mcols, Eef_1); double** Eef_2; createMatrix(mrows, mcols, Eef_2);
		double** Eef_3; createMatrix(mrows, mcols, Eef_3); double** Eef_4; createMatrix(mrows, mcols, Eef_4); double** Eef_5; createMatrix(mrows, mcols, Eef_5);
		double** Eif_0; createMatrix(mrows, mcols, Eif_0); double** Eif_1; createMatrix(mrows, mcols, Eif_1); double** Eif_2; createMatrix(mrows, mcols, Eif_2);
		double** Eif_3; createMatrix(mrows, mcols, Eif_3); double** Eif_4; createMatrix(mrows, mcols, Eif_4); double** Eif_5; createMatrix(mrows, mcols, Eif_5);
		double** Erf_0; createMatrix(mrows, mcols, Erf_0); double** Erf_1; createMatrix(mrows, mcols, Erf_1); double** Erf_2; createMatrix(mrows, mcols, Erf_2);
		double** Erf_3; createMatrix(mrows, mcols, Erf_3); double** Erf_4; createMatrix(mrows, mcols, Erf_4); double** Erf_5; createMatrix(mrows, mcols, Erf_5);
		
		compute_coefx(coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, h_x);
		compute_coefy(coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, h_y);

		compute_F1_alphaC(F1_alphac1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F1_alphaC(F1_alphac2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F1_alphaC(F1_alphac3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F1_alphaC(F1_alphac4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F1_alphaC(F1_alphac5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);

		compute_F2_alphaC(F2_alphac1, GQxi_1, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F2_alphaC(F2_alphac2, GQxi_2, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F2_alphaC(F2_alphac3, GQxi_3, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F2_alphaC(F2_alphac4, GQxi_4, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		compute_F2_alphaC(F2_alphac5, GQxi_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5);
		
		compute_F1_alphaC_max(F1_alphac, F1_alphac1, F1_alphac2, F1_alphac3, F1_alphac4, F1_alphac5);
		compute_F2_alphaC_max(F2_alphac, F2_alphac1, F2_alphac2, F2_alphac3, F2_alphac4, F2_alphac5);

		if (t > t_fin - 5 * tau)
		{
			double alphacx_min = compute_alphaCx_min(F1_alphac);
			double alphacx_max = compute_alphaCx_max(F1_alphac);
			double alphacy_min = compute_alphaCy_min(F2_alphac);
			double alphacy_max = compute_alphaCy_max(F2_alphac);
			cout << endl;
			cout << alphacx_min << " , " << alphacx_max << endl;
			cout << alphacy_min << " , " << alphacy_max << endl;
			cout << endl;
		}

		double** rhofTP_0; createMatrix(mrows, mcols, rhofTP_0); double** rhofTP_1; createMatrix(mrows, mcols, rhofTP_1); double** rhofTP_2; createMatrix(mrows, mcols, rhofTP_2);
		double** rhofTP_3; createMatrix(mrows, mcols, rhofTP_3); double** rhofTP_4; createMatrix(mrows, mcols, rhofTP_4); double** rhofTP_5; createMatrix(mrows, mcols, rhofTP_5);
		double** m1fTP_0; createMatrix(mrows, mcols, m1fTP_0); double** m1fTP_1; createMatrix(mrows, mcols, m1fTP_1); double** m1fTP_2; createMatrix(mrows, mcols, m1fTP_2);
		double** m1fTP_3; createMatrix(mrows, mcols, m1fTP_3); double** m1fTP_4; createMatrix(mrows, mcols, m1fTP_4); double** m1fTP_5; createMatrix(mrows, mcols, m1fTP_5);
		double** m2fTP_0; createMatrix(mrows, mcols, m2fTP_0); double** m2fTP_1; createMatrix(mrows, mcols, m2fTP_1); double** m2fTP_2; createMatrix(mrows, mcols, m2fTP_2);
		double** m2fTP_3; createMatrix(mrows, mcols, m2fTP_3); double** m2fTP_4; createMatrix(mrows, mcols, m2fTP_4); double** m2fTP_5; createMatrix(mrows, mcols, m2fTP_5);
		double** EefTP_0; createMatrix(mrows, mcols, EefTP_0); double** EefTP_1; createMatrix(mrows, mcols, EefTP_1); double** EefTP_2; createMatrix(mrows, mcols, EefTP_2);
		double** EefTP_3; createMatrix(mrows, mcols, EefTP_3); double** EefTP_4; createMatrix(mrows, mcols, EefTP_4); double** EefTP_5; createMatrix(mrows, mcols, EefTP_5);
		double** EifTP_0; createMatrix(mrows, mcols, EifTP_0); double** EifTP_1; createMatrix(mrows, mcols, EifTP_1); double** EifTP_2; createMatrix(mrows, mcols, EifTP_2);
		double** EifTP_3; createMatrix(mrows, mcols, EifTP_3); double** EifTP_4; createMatrix(mrows, mcols, EifTP_4); double** EifTP_5; createMatrix(mrows, mcols, EifTP_5);
		double** ErfTP_0; createMatrix(mrows, mcols, ErfTP_0); double** ErfTP_1; createMatrix(mrows, mcols, ErfTP_1); double** ErfTP_2; createMatrix(mrows, mcols, ErfTP_2);
		double** ErfTP_3; createMatrix(mrows, mcols, ErfTP_3); double** ErfTP_4; createMatrix(mrows, mcols, ErfTP_4); double** ErfTP_5; createMatrix(mrows, mcols, ErfTP_5);

		space_d(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, F1_alphac1, F2_alphac, h_x, h_y, t_1);
		get_L2_projection_2D(rhofTP_0, rhofTP_1, rhofTP_2, rhofTP_3, rhofTP_4, rhofTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_f(rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, rhofTP_0, rhofTP_1, rhofTP_2, rhofTP_3, rhofTP_4, rhofTP_5, tau);

		space_m1(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y, t_1);
		get_L2_projection_2D(m1fTP_0, m1fTP_1, m1fTP_2, m1fTP_3, m1fTP_4, m1fTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_f(m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m1fTP_0, m1fTP_1, m1fTP_2, m1fTP_3, m1fTP_4, m1fTP_5, tau);

		space_m2(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y, t_1);
		get_L2_projection_2D(m2fTP_0, m2fTP_1, m2fTP_2, m2fTP_3, m2fTP_4, m2fTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_f(m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, m2fTP_0, m2fTP_1, m2fTP_2, m2fTP_3, m2fTP_4, m2fTP_5, tau);

		space_e(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y, t_1, Mode);
		get_L2_projection_2D(EefTP_0, EefTP_1, EefTP_2, EefTP_3, EefTP_4, EefTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Ee_vec(M_rows); Eigen::VectorXd Eef_vec(M_rows); Eigen::VectorXd EeN_f(M_rows);
		IMEX3_FS(Eef_vec, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Ee_vec, EeN_f, M, m_rows, M_rows, EefTP_0, EefTP_1, EefTP_2, EefTP_3, EefTP_4, EefTP_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, tau, a0, solver);

		space_i(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y, t_1, Mode);
		get_L2_projection_2D(EifTP_0, EifTP_1, EifTP_2, EifTP_3, EifTP_4, EifTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Ei_vec(M_rows); Eigen::VectorXd Eif_vec(M_rows); Eigen::VectorXd EiN_f(M_rows);
		IMEX3_FS(Eif_vec, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Ei_vec, EiN_f, M, m_rows, M_rows, EifTP_0, EifTP_1, EifTP_2, EifTP_3, EifTP_4, EifTP_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, tau, a0, solver);

		space_r(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, F1_alphac, F2_alphac, h_x, h_y, t_1, Mode);
		get_L2_projection_2D(ErfTP_0, ErfTP_1, ErfTP_2, ErfTP_3, ErfTP_4, ErfTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Er_vec(M_rows); Eigen::VectorXd Erf_vec(M_rows); Eigen::VectorXd ErN_f(M_rows);
		IMEX3_FS(Erf_vec, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, Er_vec, ErN_f, M, m_rows, M_rows, ErfTP_0, ErfTP_1, ErfTP_2, ErfTP_3, ErfTP_4, ErfTP_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, tau, a0, solver);

		// update boundary
		update_boundary(rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5); update_boundary(m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5);
		update_boundary(m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5); update_boundary(Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5);
		update_boundary(Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5); update_boundary(Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);

		double* Ee0_vec = new double[m_rows]; double* Ee1_vec = new double[m_rows]; double* Ee2_vec = new double[m_rows];
		double* Ee3_vec = new double[m_rows]; double* Ee4_vec = new double[m_rows]; double* Ee5_vec = new double[m_rows];

		double* Ei0_vec = new double[m_rows]; double* Ei1_vec = new double[m_rows]; double* Ei2_vec = new double[m_rows];
		double* Ei3_vec = new double[m_rows]; double* Ei4_vec = new double[m_rows]; double* Ei5_vec = new double[m_rows];

		double* Er0_vec = new double[m_rows]; double* Er1_vec = new double[m_rows]; double* Er2_vec = new double[m_rows];
		double* Er3_vec = new double[m_rows]; double* Er4_vec = new double[m_rows]; double* Er5_vec = new double[m_rows];

		if (OE_index == 1)
		{
			OE_compute_system_sigma(sigma_0, sigma_1, sigma_2, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, h_x, h_y);
			
			OE_procedure(rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, sigma_0, sigma_1, sigma_2, tau);

			Matrix_to_vector(Ee0_vec, Eef_0, N_y + 2, N_x + 2); Matrix_to_vector(Ee1_vec, Eef_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Ee2_vec, Eef_2, N_y + 2, N_x + 2); Matrix_to_vector(Ee3_vec, Eef_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Ee4_vec, Eef_4, N_y + 2, N_x + 2); Matrix_to_vector(Ee5_vec, Eef_5, N_y + 2, N_x + 2);

			Matrix_to_vector(Ei0_vec, Eif_0, N_y + 2, N_x + 2); Matrix_to_vector(Ei1_vec, Eif_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Ei2_vec, Eif_2, N_y + 2, N_x + 2); Matrix_to_vector(Ei3_vec, Eif_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Ei4_vec, Eif_4, N_y + 2, N_x + 2); Matrix_to_vector(Ei5_vec, Eif_5, N_y + 2, N_x + 2);

			Matrix_to_vector(Er0_vec, Erf_0, N_y + 2, N_x + 2); Matrix_to_vector(Er1_vec, Erf_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Er2_vec, Erf_2, N_y + 2, N_x + 2); Matrix_to_vector(Er3_vec, Erf_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Er4_vec, Erf_4, N_y + 2, N_x + 2); Matrix_to_vector(Er5_vec, Erf_5, N_y + 2, N_x + 2);

			for (int j = 0; j < m_rows; j++)
			{
				Eef_vec(j) = Ee0_vec[j]; Eef_vec(j + m_rows) = Ee1_vec[j]; Eef_vec(j + 2 * m_rows) = Ee2_vec[j];
				Eef_vec(j + 3 * m_rows) = Ee3_vec[j]; Eef_vec(j + 4 * m_rows) = Ee4_vec[j]; Eef_vec(j + 5 * m_rows) = Ee5_vec[j];

				Eif_vec(j) = Ei0_vec[j]; Eif_vec(j + m_rows) = Ei1_vec[j]; Eif_vec(j + 2 * m_rows) = Ei2_vec[j];
				Eif_vec(j + 3 * m_rows) = Ei3_vec[j]; Eif_vec(j + 4 * m_rows) = Ei4_vec[j]; Eif_vec(j + 5 * m_rows) = Ei5_vec[j];

				Erf_vec(j) = Er0_vec[j]; Erf_vec(j + m_rows) = Er1_vec[j]; Erf_vec(j + 2 * m_rows) = Er2_vec[j];
				Erf_vec(j + 3 * m_rows) = Er3_vec[j]; Erf_vec(j + 4 * m_rows) = Er4_vec[j]; Erf_vec(j + 5 * m_rows) = Er5_vec[j];
			}

			// update boundary
			update_boundary(rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5); update_boundary(m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5);
			update_boundary(m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5); update_boundary(Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5);
			update_boundary(Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5); update_boundary(Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		}

		// Second step of IMEX3
		double** rhos_0; createMatrix(mrows, mcols, rhos_0); double** rhos_1; createMatrix(mrows, mcols, rhos_1); double** rhos_2; createMatrix(mrows, mcols, rhos_2);
		double** rhos_3; createMatrix(mrows, mcols, rhos_3); double** rhos_4; createMatrix(mrows, mcols, rhos_4); double** rhos_5; createMatrix(mrows, mcols, rhos_5);
		double** m1s_0; createMatrix(mrows, mcols, m1s_0); double** m1s_1; createMatrix(mrows, mcols, m1s_1); double** m1s_2; createMatrix(mrows, mcols, m1s_2);
		double** m1s_3; createMatrix(mrows, mcols, m1s_3); double** m1s_4; createMatrix(mrows, mcols, m1s_4); double** m1s_5; createMatrix(mrows, mcols, m1s_5);
		double** m2s_0; createMatrix(mrows, mcols, m2s_0); double** m2s_1; createMatrix(mrows, mcols, m2s_1); double** m2s_2; createMatrix(mrows, mcols, m2s_2);
		double** m2s_3; createMatrix(mrows, mcols, m2s_3); double** m2s_4; createMatrix(mrows, mcols, m2s_4); double** m2s_5; createMatrix(mrows, mcols, m2s_5);
		double** Ees_0; createMatrix(mrows, mcols, Ees_0); double** Ees_1; createMatrix(mrows, mcols, Ees_1); double** Ees_2; createMatrix(mrows, mcols, Ees_2);
		double** Ees_3; createMatrix(mrows, mcols, Ees_3); double** Ees_4; createMatrix(mrows, mcols, Ees_4); double** Ees_5; createMatrix(mrows, mcols, Ees_5);
		double** Eis_0; createMatrix(mrows, mcols, Eis_0); double** Eis_1; createMatrix(mrows, mcols, Eis_1); double** Eis_2; createMatrix(mrows, mcols, Eis_2);
		double** Eis_3; createMatrix(mrows, mcols, Eis_3); double** Eis_4; createMatrix(mrows, mcols, Eis_4); double** Eis_5; createMatrix(mrows, mcols, Eis_5);
		double** Ers_0; createMatrix(mrows, mcols, Ers_0); double** Ers_1; createMatrix(mrows, mcols, Ers_1); double** Ers_2; createMatrix(mrows, mcols, Ers_2);
		double** Ers_3; createMatrix(mrows, mcols, Ers_3); double** Ers_4; createMatrix(mrows, mcols, Ers_4); double** Ers_5; createMatrix(mrows, mcols, Ers_5);

		compute_coefx(coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, h_x);
		compute_coefy(coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, h_y);

		compute_F1_alphaC(F1_alphac1, GQxi_1, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F1_alphaC(F1_alphac2, GQxi_2, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F1_alphaC(F1_alphac3, GQxi_3, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F1_alphaC(F1_alphac4, GQxi_4, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F1_alphaC(F1_alphac5, GQxi_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);

		compute_F2_alphaC(F2_alphac1, GQxi_1, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F2_alphaC(F2_alphac2, GQxi_2, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F2_alphaC(F2_alphac3, GQxi_3, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F2_alphaC(F2_alphac4, GQxi_4, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);
		compute_F2_alphaC(F2_alphac5, GQxi_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5);

		compute_F1_alphaC_max(F1_alphac, F1_alphac1, F1_alphac2, F1_alphac3, F1_alphac4, F1_alphac5);
		compute_F2_alphaC_max(F2_alphac, F2_alphac1, F2_alphac2, F2_alphac3, F2_alphac4, F2_alphac5);

		double** rhosTP_0; createMatrix(mrows, mcols, rhosTP_0); double** rhosTP_1; createMatrix(mrows, mcols, rhosTP_1); double** rhosTP_2; createMatrix(mrows, mcols, rhosTP_2);
		double** rhosTP_3; createMatrix(mrows, mcols, rhosTP_3); double** rhosTP_4; createMatrix(mrows, mcols, rhosTP_4); double** rhosTP_5; createMatrix(mrows, mcols, rhosTP_5);
		double** m1sTP_0; createMatrix(mrows, mcols, m1sTP_0); double** m1sTP_1; createMatrix(mrows, mcols, m1sTP_1); double** m1sTP_2; createMatrix(mrows, mcols, m1sTP_2);
		double** m1sTP_3; createMatrix(mrows, mcols, m1sTP_3); double** m1sTP_4; createMatrix(mrows, mcols, m1sTP_4); double** m1sTP_5; createMatrix(mrows, mcols, m1sTP_5);
		double** m2sTP_0; createMatrix(mrows, mcols, m2sTP_0); double** m2sTP_1; createMatrix(mrows, mcols, m2sTP_1); double** m2sTP_2; createMatrix(mrows, mcols, m2sTP_2);
		double** m2sTP_3; createMatrix(mrows, mcols, m2sTP_3); double** m2sTP_4; createMatrix(mrows, mcols, m2sTP_4); double** m2sTP_5; createMatrix(mrows, mcols, m2sTP_5);
		double** EesTP_0; createMatrix(mrows, mcols, EesTP_0); double** EesTP_1; createMatrix(mrows, mcols, EesTP_1); double** EesTP_2; createMatrix(mrows, mcols, EesTP_2);
		double** EesTP_3; createMatrix(mrows, mcols, EesTP_3); double** EesTP_4; createMatrix(mrows, mcols, EesTP_4); double** EesTP_5; createMatrix(mrows, mcols, EesTP_5);
		double** EisTP_0; createMatrix(mrows, mcols, EisTP_0); double** EisTP_1; createMatrix(mrows, mcols, EisTP_1); double** EisTP_2; createMatrix(mrows, mcols, EisTP_2);
		double** EisTP_3; createMatrix(mrows, mcols, EisTP_3); double** EisTP_4; createMatrix(mrows, mcols, EisTP_4); double** EisTP_5; createMatrix(mrows, mcols, EisTP_5);
		double** ErsTP_0; createMatrix(mrows, mcols, ErsTP_0); double** ErsTP_1; createMatrix(mrows, mcols, ErsTP_1); double** ErsTP_2; createMatrix(mrows, mcols, ErsTP_2);
		double** ErsTP_3; createMatrix(mrows, mcols, ErsTP_3); double** ErsTP_4; createMatrix(mrows, mcols, ErsTP_4); double** ErsTP_5; createMatrix(mrows, mcols, ErsTP_5);

		space_d(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, F1_alphac, F2_alphac, h_x, h_y, t_2);
		get_L2_projection_2D(rhosTP_0, rhosTP_1, rhosTP_2, rhosTP_3, rhosTP_4, rhosTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_s(rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, rhofTP_0, rhofTP_1, rhofTP_2, rhofTP_3, rhofTP_4, rhofTP_5, rhosTP_0, rhosTP_1, rhosTP_2, rhosTP_3, rhosTP_4, rhosTP_5, tau, a1);

		space_m1(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, F1_alphac, F2_alphac, h_x, h_y, t_2);
		get_L2_projection_2D(m1sTP_0, m1sTP_1, m1sTP_2, m1sTP_3, m1sTP_4, m1sTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_s(m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m1fTP_0, m1fTP_1, m1fTP_2, m1fTP_3, m1fTP_4, m1fTP_5, m1sTP_0, m1sTP_1, m1sTP_2, m1sTP_3, m1sTP_4, m1sTP_5, tau, a1);

		space_m2(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, F1_alphac, F2_alphac, h_x, h_y, t_2);
		get_L2_projection_2D(m2sTP_0, m2sTP_1, m2sTP_2, m2sTP_3, m2sTP_4, m2sTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_s(m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, m2fTP_0, m2fTP_1, m2fTP_2, m2fTP_3, m2fTP_4, m2fTP_5, m2sTP_0, m2sTP_1, m2sTP_2, m2sTP_3, m2sTP_4, m2sTP_5, tau, a1);

		space_e(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, F1_alphac, F2_alphac, h_x, h_y, t_2, Mode);
		get_L2_projection_2D(EesTP_0, EesTP_1, EesTP_2, EesTP_3, EesTP_4, EesTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Ees_vec(M_rows); Eigen::VectorXd EeN_s(M_rows);
		IMEX3_SS(Ees_vec, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, EeN_s, M, m_rows, M_rows, EesTP_0, EesTP_1, EesTP_2, EesTP_3, EesTP_4, EesTP_5, Ee_vec, Eef_vec, EeN_f, tau, a0, a1, solver);

		space_i(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, F1_alphac, F2_alphac, h_x, h_y, t_2, Mode);
		get_L2_projection_2D(EisTP_0, EisTP_1, EisTP_2, EisTP_3, EisTP_4, EisTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Eis_vec(M_rows); Eigen::VectorXd EiN_s(M_rows);
		IMEX3_SS(Eis_vec, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, EiN_s, M, m_rows, M_rows, EisTP_0, EisTP_1, EisTP_2, EisTP_3, EisTP_4, EisTP_5, Ei_vec, Eif_vec, EiN_f, tau, a0, a1, solver);

		space_r(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, F1_alphac, F2_alphac, h_x, h_y, t_2, Mode);
		get_L2_projection_2D(ErsTP_0, ErsTP_1, ErsTP_2, ErsTP_3, ErsTP_4, ErsTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Ers_vec(M_rows); Eigen::VectorXd ErN_s(M_rows);
		IMEX3_SS(Ers_vec, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, ErN_s, M, m_rows, M_rows, ErsTP_0, ErsTP_1, ErsTP_2, ErsTP_3, ErsTP_4, ErsTP_5, Er_vec, Erf_vec, ErN_f, tau, a0, a1, solver);

		// update boundary
		update_boundary(rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5); update_boundary(m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5);
		update_boundary(m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5); update_boundary(Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5);
		update_boundary(Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5); update_boundary(Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(sigma_0, sigma_1, sigma_2, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, h_x, h_y);

			OE_procedure(rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, sigma_0, sigma_1, sigma_2, tau);

			Matrix_to_vector(Ee0_vec, Ees_0, N_y + 2, N_x + 2); Matrix_to_vector(Ee1_vec, Ees_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Ee2_vec, Ees_2, N_y + 2, N_x + 2); Matrix_to_vector(Ee3_vec, Ees_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Ee4_vec, Ees_4, N_y + 2, N_x + 2); Matrix_to_vector(Ee5_vec, Ees_5, N_y + 2, N_x + 2);

			Matrix_to_vector(Ei0_vec, Eis_0, N_y + 2, N_x + 2); Matrix_to_vector(Ei1_vec, Eis_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Ei2_vec, Eis_2, N_y + 2, N_x + 2); Matrix_to_vector(Ei3_vec, Eis_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Ei4_vec, Eis_4, N_y + 2, N_x + 2); Matrix_to_vector(Ei5_vec, Eis_5, N_y + 2, N_x + 2);

			Matrix_to_vector(Er0_vec, Ers_0, N_y + 2, N_x + 2); Matrix_to_vector(Er1_vec, Ers_1, N_y + 2, N_x + 2);
			Matrix_to_vector(Er2_vec, Ers_2, N_y + 2, N_x + 2); Matrix_to_vector(Er3_vec, Ers_3, N_y + 2, N_x + 2);
			Matrix_to_vector(Er4_vec, Ers_4, N_y + 2, N_x + 2); Matrix_to_vector(Er5_vec, Ers_5, N_y + 2, N_x + 2);

			for (int j = 0; j < m_rows; j++)
			{
				Ees_vec(j) = Ee0_vec[j]; Ees_vec(j + m_rows) = Ee1_vec[j]; Ees_vec(j + 2 * m_rows) = Ee2_vec[j];
				Ees_vec(j + 3 * m_rows) = Ee3_vec[j]; Ees_vec(j + 4 * m_rows) = Ee4_vec[j]; Ees_vec(j + 5 * m_rows) = Ee5_vec[j];

				Eis_vec(j) = Ei0_vec[j]; Eis_vec(j + m_rows) = Ei1_vec[j]; Eis_vec(j + 2 * m_rows) = Ei2_vec[j];
				Eis_vec(j + 3 * m_rows) = Ei3_vec[j]; Eis_vec(j + 4 * m_rows) = Ei4_vec[j]; Eis_vec(j + 5 * m_rows) = Ei5_vec[j];

				Ers_vec(j) = Er0_vec[j]; Ers_vec(j + m_rows) = Er1_vec[j]; Ers_vec(j + 2 * m_rows) = Er2_vec[j];
				Ers_vec(j + 3 * m_rows) = Er3_vec[j]; Ers_vec(j + 4 * m_rows) = Er4_vec[j]; Ers_vec(j + 5 * m_rows) = Er5_vec[j];
			}

			// update boundary
			update_boundary(rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5); update_boundary(m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5);
			update_boundary(m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5); update_boundary(Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5);
			update_boundary(Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5); update_boundary(Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		}

		delete[] Ee0_vec; delete[] Ee1_vec; delete[] Ee2_vec; delete[] Ee3_vec; delete[] Ee4_vec; delete[] Ee5_vec;
		delete[] Ei0_vec; delete[] Ei1_vec; delete[] Ei2_vec; delete[] Ei3_vec; delete[] Ei4_vec; delete[] Ei5_vec;
		delete[] Er0_vec; delete[] Er1_vec; delete[] Er2_vec; delete[] Er3_vec; delete[] Er4_vec; delete[] Er5_vec;


		// Third step of IMEX3
		double** rhot_0; createMatrix(mrows, mcols, rhot_0); double** rhot_1; createMatrix(mrows, mcols, rhot_1); double** rhot_2; createMatrix(mrows, mcols, rhot_2);
		double** rhot_3; createMatrix(mrows, mcols, rhot_3); double** rhot_4; createMatrix(mrows, mcols, rhot_4); double** rhot_5; createMatrix(mrows, mcols, rhot_5);
		double** m1t_0; createMatrix(mrows, mcols, m1t_0); double** m1t_1; createMatrix(mrows, mcols, m1t_1); double** m1t_2; createMatrix(mrows, mcols, m1t_2);
		double** m1t_3; createMatrix(mrows, mcols, m1t_3); double** m1t_4; createMatrix(mrows, mcols, m1t_4); double** m1t_5; createMatrix(mrows, mcols, m1t_5);
		double** m2t_0; createMatrix(mrows, mcols, m2t_0); double** m2t_1; createMatrix(mrows, mcols, m2t_1); double** m2t_2; createMatrix(mrows, mcols, m2t_2);
		double** m2t_3; createMatrix(mrows, mcols, m2t_3); double** m2t_4; createMatrix(mrows, mcols, m2t_4); double** m2t_5; createMatrix(mrows, mcols, m2t_5);
		double** Eet_0; createMatrix(mrows, mcols, Eet_0); double** Eet_1; createMatrix(mrows, mcols, Eet_1); double** Eet_2; createMatrix(mrows, mcols, Eet_2);
		double** Eet_3; createMatrix(mrows, mcols, Eet_3); double** Eet_4; createMatrix(mrows, mcols, Eet_4); double** Eet_5; createMatrix(mrows, mcols, Eet_5);
		double** Eit_0; createMatrix(mrows, mcols, Eit_0); double** Eit_1; createMatrix(mrows, mcols, Eit_1); double** Eit_2; createMatrix(mrows, mcols, Eit_2);
		double** Eit_3; createMatrix(mrows, mcols, Eit_3); double** Eit_4; createMatrix(mrows, mcols, Eit_4); double** Eit_5; createMatrix(mrows, mcols, Eit_5);
		double** Ert_0; createMatrix(mrows, mcols, Ert_0); double** Ert_1; createMatrix(mrows, mcols, Ert_1); double** Ert_2; createMatrix(mrows, mcols, Ert_2);
		double** Ert_3; createMatrix(mrows, mcols, Ert_3); double** Ert_4; createMatrix(mrows, mcols, Ert_4); double** Ert_5; createMatrix(mrows, mcols, Ert_5);

		compute_coefx(coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, h_x);
		compute_coefy(coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, h_y);

		compute_F1_alphaC(F1_alphac1, GQxi_1, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F1_alphaC(F1_alphac2, GQxi_2, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F1_alphaC(F1_alphac3, GQxi_3, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F1_alphaC(F1_alphac4, GQxi_4, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F1_alphaC(F1_alphac5, GQxi_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);

		compute_F2_alphaC(F2_alphac1, GQxi_1, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F2_alphaC(F2_alphac2, GQxi_2, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F2_alphaC(F2_alphac3, GQxi_3, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F2_alphaC(F2_alphac4, GQxi_4, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);
		compute_F2_alphaC(F2_alphac5, GQxi_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5);

		compute_F1_alphaC_max(F1_alphac, F1_alphac1, F1_alphac2, F1_alphac3, F1_alphac4, F1_alphac5);
		compute_F2_alphaC_max(F2_alphac, F2_alphac1, F2_alphac2, F2_alphac3, F2_alphac4, F2_alphac5);

		double** rhotTP_0; createMatrix(mrows, mcols, rhotTP_0); double** rhotTP_1; createMatrix(mrows, mcols, rhotTP_1); double** rhotTP_2; createMatrix(mrows, mcols, rhotTP_2);
		double** rhotTP_3; createMatrix(mrows, mcols, rhotTP_3); double** rhotTP_4; createMatrix(mrows, mcols, rhotTP_4); double** rhotTP_5; createMatrix(mrows, mcols, rhotTP_5);
		double** m1tTP_0; createMatrix(mrows, mcols, m1tTP_0); double** m1tTP_1; createMatrix(mrows, mcols, m1tTP_1); double** m1tTP_2; createMatrix(mrows, mcols, m1tTP_2);
		double** m1tTP_3; createMatrix(mrows, mcols, m1tTP_3); double** m1tTP_4; createMatrix(mrows, mcols, m1tTP_4); double** m1tTP_5; createMatrix(mrows, mcols, m1tTP_5);
		double** m2tTP_0; createMatrix(mrows, mcols, m2tTP_0); double** m2tTP_1; createMatrix(mrows, mcols, m2tTP_1); double** m2tTP_2; createMatrix(mrows, mcols, m2tTP_2);
		double** m2tTP_3; createMatrix(mrows, mcols, m2tTP_3); double** m2tTP_4; createMatrix(mrows, mcols, m2tTP_4); double** m2tTP_5; createMatrix(mrows, mcols, m2tTP_5);
		double** EetTP_0; createMatrix(mrows, mcols, EetTP_0); double** EetTP_1; createMatrix(mrows, mcols, EetTP_1); double** EetTP_2; createMatrix(mrows, mcols, EetTP_2);
		double** EetTP_3; createMatrix(mrows, mcols, EetTP_3); double** EetTP_4; createMatrix(mrows, mcols, EetTP_4); double** EetTP_5; createMatrix(mrows, mcols, EetTP_5);
		double** EitTP_0; createMatrix(mrows, mcols, EitTP_0); double** EitTP_1; createMatrix(mrows, mcols, EitTP_1); double** EitTP_2; createMatrix(mrows, mcols, EitTP_2);
		double** EitTP_3; createMatrix(mrows, mcols, EitTP_3); double** EitTP_4; createMatrix(mrows, mcols, EitTP_4); double** EitTP_5; createMatrix(mrows, mcols, EitTP_5);
		double** ErtTP_0; createMatrix(mrows, mcols, ErtTP_0); double** ErtTP_1; createMatrix(mrows, mcols, ErtTP_1); double** ErtTP_2; createMatrix(mrows, mcols, ErtTP_2);
		double** ErtTP_3; createMatrix(mrows, mcols, ErtTP_3); double** ErtTP_4; createMatrix(mrows, mcols, ErtTP_4); double** ErtTP_5; createMatrix(mrows, mcols, ErtTP_5);

		space_d(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, F1_alphac, F2_alphac, h_x, h_y, t_3);
		get_L2_projection_2D(rhotTP_0, rhotTP_1, rhotTP_2, rhotTP_3, rhotTP_4, rhotTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_t(rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, rhosTP_0, rhosTP_1, rhosTP_2, rhosTP_3, rhosTP_4, rhosTP_5, rhotTP_0, rhotTP_1, rhotTP_2, rhotTP_3, rhotTP_4, rhotTP_5, tau, a2);

		space_m1(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, F1_alphac, F2_alphac, h_x, h_y, t_3);
		get_L2_projection_2D(m1tTP_0, m1tTP_1, m1tTP_2, m1tTP_3, m1tTP_4, m1tTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_t(m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m1sTP_0, m1sTP_1, m1sTP_2, m1sTP_3, m1sTP_4, m1sTP_5, m1tTP_0, m1tTP_1, m1tTP_2, m1tTP_3, m1tTP_4, m1tTP_5, tau, a2);

		space_m2(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, F1_alphac, F2_alphac, h_x, h_y, t_3);
		get_L2_projection_2D(m2tTP_0, m2tTP_1, m2tTP_2, m2tTP_3, m2tTP_4, m2tTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_t(m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, m2sTP_0, m2sTP_1, m2sTP_2, m2sTP_3, m2sTP_4, m2sTP_5, m2tTP_0, m2tTP_1, m2tTP_2, m2tTP_3, m2tTP_4, m2tTP_5, tau, a2);

		space_e(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, F1_alphac, F2_alphac, h_x, h_y, t_3, Mode);
		get_L2_projection_2D(EetTP_0, EetTP_1, EetTP_2, EetTP_3, EetTP_4, EetTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Eet_vec(M_rows); Eigen::VectorXd EeN_t(M_rows);
		IMEX3_TS(Eet_vec, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, EeN_t, M, m_rows, M_rows, EetTP_0, EetTP_1, EetTP_2, EetTP_3, EetTP_4, EetTP_5, Ee_vec, Eef_vec, Ees_vec, EeN_s, tau, a0, b1, b2, a2, solver);

		space_i(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, F1_alphac, F2_alphac, h_x, h_y, t_3, Mode);
		get_L2_projection_2D(EitTP_0, EitTP_1, EitTP_2, EitTP_3, EitTP_4, EitTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Eit_vec(M_rows); Eigen::VectorXd EiN_t(M_rows);
		IMEX3_TS(Eit_vec, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, EiN_t, M, m_rows, M_rows, EitTP_0, EitTP_1, EitTP_2, EitTP_3, EitTP_4, EitTP_5, Ei_vec, Eif_vec, Eis_vec, EiN_s, tau, a0, b1, b2, a2, solver);

		space_r(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, F1_alphac, F2_alphac, h_x, h_y, t_3, Mode);
		get_L2_projection_2D(ErtTP_0, ErtTP_1, ErtTP_2, ErtTP_3, ErtTP_4, ErtTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		Eigen::VectorXd Ert_vec(M_rows); Eigen::VectorXd ErN_t(M_rows);
		IMEX3_TS(Ert_vec, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, ErN_t, M, m_rows, M_rows, ErtTP_0, ErtTP_1, ErtTP_2, ErtTP_3, ErtTP_4, ErtTP_5, Er_vec, Erf_vec, Ers_vec, ErN_s, tau, a0, b1, b2, a2, solver);

		// update boundary
		update_boundary(rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5); update_boundary(m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5);
		update_boundary(m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5); update_boundary(Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5);
		update_boundary(Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5); update_boundary(Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(sigma_0, sigma_1, sigma_2, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, h_x, h_y);

			OE_procedure(rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, sigma_0, sigma_1, sigma_2, tau);


			// update boundary
			update_boundary(rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5); update_boundary(m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5);
			update_boundary(m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5); update_boundary(Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5);
			update_boundary(Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5); update_boundary(Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		}

		// Last step of IMEX3
		double** rhol_0; createMatrix(mrows, mcols, rhol_0); double** rhol_1; createMatrix(mrows, mcols, rhol_1); double** rhol_2; createMatrix(mrows, mcols, rhol_2);
		double** rhol_3; createMatrix(mrows, mcols, rhol_3); double** rhol_4; createMatrix(mrows, mcols, rhol_4); double** rhol_5; createMatrix(mrows, mcols, rhol_5);
		double** m1l_0; createMatrix(mrows, mcols, m1l_0); double** m1l_1; createMatrix(mrows, mcols, m1l_1); double** m1l_2; createMatrix(mrows, mcols, m1l_2);
		double** m1l_3; createMatrix(mrows, mcols, m1l_3); double** m1l_4; createMatrix(mrows, mcols, m1l_4); double** m1l_5; createMatrix(mrows, mcols, m1l_5);
		double** m2l_0; createMatrix(mrows, mcols, m2l_0); double** m2l_1; createMatrix(mrows, mcols, m2l_1); double** m2l_2; createMatrix(mrows, mcols, m2l_2);
		double** m2l_3; createMatrix(mrows, mcols, m2l_3); double** m2l_4; createMatrix(mrows, mcols, m2l_4); double** m2l_5; createMatrix(mrows, mcols, m2l_5);
		double** Eel_0; createMatrix(mrows, mcols, Eel_0); double** Eel_1; createMatrix(mrows, mcols, Eel_1); double** Eel_2; createMatrix(mrows, mcols, Eel_2);
		double** Eel_3; createMatrix(mrows, mcols, Eel_3); double** Eel_4; createMatrix(mrows, mcols, Eel_4); double** Eel_5; createMatrix(mrows, mcols, Eel_5);
		double** Eil_0; createMatrix(mrows, mcols, Eil_0); double** Eil_1; createMatrix(mrows, mcols, Eil_1); double** Eil_2; createMatrix(mrows, mcols, Eil_2);
		double** Eil_3; createMatrix(mrows, mcols, Eil_3); double** Eil_4; createMatrix(mrows, mcols, Eil_4); double** Eil_5; createMatrix(mrows, mcols, Eil_5);
		double** Erl_0; createMatrix(mrows, mcols, Erl_0); double** Erl_1; createMatrix(mrows, mcols, Erl_1); double** Erl_2; createMatrix(mrows, mcols, Erl_2);
		double** Erl_3; createMatrix(mrows, mcols, Erl_3); double** Erl_4; createMatrix(mrows, mcols, Erl_4); double** Erl_5; createMatrix(mrows, mcols, Erl_5);

		compute_coefx(coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, h_x);
		compute_coefy(coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, h_y);

		compute_F1_alphaC(F1_alphac1, GQxi_1, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F1_alphaC(F1_alphac2, GQxi_2, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F1_alphaC(F1_alphac3, GQxi_3, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F1_alphaC(F1_alphac4, GQxi_4, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F1_alphaC(F1_alphac5, GQxi_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);

		compute_F2_alphaC(F2_alphac1, GQxi_1, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F2_alphaC(F2_alphac2, GQxi_2, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F2_alphaC(F2_alphac3, GQxi_3, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F2_alphaC(F2_alphac4, GQxi_4, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);
		compute_F2_alphaC(F2_alphac5, GQxi_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5);

		compute_F1_alphaC_max(F1_alphac, F1_alphac1, F1_alphac2, F1_alphac3, F1_alphac4, F1_alphac5);
		compute_F2_alphaC_max(F2_alphac, F2_alphac1, F2_alphac2, F2_alphac3, F2_alphac4, F2_alphac5);

		double** rholTP_0; createMatrix(mrows, mcols, rholTP_0); double** rholTP_1; createMatrix(mrows, mcols, rholTP_1); double** rholTP_2; createMatrix(mrows, mcols, rholTP_2);
		double** rholTP_3; createMatrix(mrows, mcols, rholTP_3); double** rholTP_4; createMatrix(mrows, mcols, rholTP_4); double** rholTP_5; createMatrix(mrows, mcols, rholTP_5);
		double** m1lTP_0; createMatrix(mrows, mcols, m1lTP_0); double** m1lTP_1; createMatrix(mrows, mcols, m1lTP_1); double** m1lTP_2; createMatrix(mrows, mcols, m1lTP_2);
		double** m1lTP_3; createMatrix(mrows, mcols, m1lTP_3); double** m1lTP_4; createMatrix(mrows, mcols, m1lTP_4); double** m1lTP_5; createMatrix(mrows, mcols, m1lTP_5);
		double** m2lTP_0; createMatrix(mrows, mcols, m2lTP_0); double** m2lTP_1; createMatrix(mrows, mcols, m2lTP_1); double** m2lTP_2; createMatrix(mrows, mcols, m2lTP_2);
		double** m2lTP_3; createMatrix(mrows, mcols, m2lTP_3); double** m2lTP_4; createMatrix(mrows, mcols, m2lTP_4); double** m2lTP_5; createMatrix(mrows, mcols, m2lTP_5);
		double** EelTP_0; createMatrix(mrows, mcols, EelTP_0); double** EelTP_1; createMatrix(mrows, mcols, EelTP_1); double** EelTP_2; createMatrix(mrows, mcols, EelTP_2);
		double** EelTP_3; createMatrix(mrows, mcols, EelTP_3); double** EelTP_4; createMatrix(mrows, mcols, EelTP_4); double** EelTP_5; createMatrix(mrows, mcols, EelTP_5);
		double** EilTP_0; createMatrix(mrows, mcols, EilTP_0); double** EilTP_1; createMatrix(mrows, mcols, EilTP_1); double** EilTP_2; createMatrix(mrows, mcols, EilTP_2);
		double** EilTP_3; createMatrix(mrows, mcols, EilTP_3); double** EilTP_4; createMatrix(mrows, mcols, EilTP_4); double** EilTP_5; createMatrix(mrows, mcols, EilTP_5);
		double** ErlTP_0; createMatrix(mrows, mcols, ErlTP_0); double** ErlTP_1; createMatrix(mrows, mcols, ErlTP_1); double** ErlTP_2; createMatrix(mrows, mcols, ErlTP_2);
		double** ErlTP_3; createMatrix(mrows, mcols, ErlTP_3); double** ErlTP_4; createMatrix(mrows, mcols, ErlTP_4); double** ErlTP_5; createMatrix(mrows, mcols, ErlTP_5);

		space_d(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, F1_alphac, F2_alphac, h_x, h_y, t_4);
		get_L2_projection_2D(rholTP_0, rholTP_1, rholTP_2, rholTP_3, rholTP_4, rholTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(rhol_0, rhol_1, rhol_2, rhol_3, rhol_4, rhol_5, rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, rhof_0, rhof_1, rhof_2, rhof_3, rhof_4, rhof_5, rhos_0, rhos_1, rhos_2, rhos_3, rhos_4, rhos_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, rhosTP_0, rhosTP_1, rhosTP_2, rhosTP_3, rhosTP_4, rhosTP_5, rhotTP_0, rhotTP_1, rhotTP_2, rhotTP_3, rhotTP_4, rhotTP_5, rholTP_0, rholTP_1, rholTP_2, rholTP_3, rholTP_4, rholTP_5, tau, b1, b2);

		space_m1(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, F1_alphac, F2_alphac, h_x, h_y, t_4);
		get_L2_projection_2D(m1lTP_0, m1lTP_1, m1lTP_2, m1lTP_3, m1lTP_4, m1lTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(m1l_0, m1l_1, m1l_2, m1l_3, m1l_4, m1l_5, m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, m1f_0, m1f_1, m1f_2, m1f_3, m1f_4, m1f_5, m1s_0, m1s_1, m1s_2, m1s_3, m1s_4, m1s_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m1sTP_0, m1sTP_1, m1sTP_2, m1sTP_3, m1sTP_4, m1sTP_5, m1tTP_0, m1tTP_1, m1tTP_2, m1tTP_3, m1tTP_4, m1tTP_5, m1lTP_0, m1lTP_1, m1lTP_2, m1lTP_3, m1lTP_4, m1lTP_5, tau, b1, b2);

		space_m2(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, F1_alphac, F2_alphac, h_x, h_y, t_4);
		get_L2_projection_2D(m2lTP_0, m2lTP_1, m2lTP_2, m2lTP_3, m2lTP_4, m2lTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(m2l_0, m2l_1, m2l_2, m2l_3, m2l_4, m2l_5, m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, m2f_0, m2f_1, m2f_2, m2f_3, m2f_4, m2f_5, m2s_0, m2s_1, m2s_2, m2s_3, m2s_4, m2s_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, m2sTP_0, m2sTP_1, m2sTP_2, m2sTP_3, m2sTP_4, m2sTP_5, m2tTP_0, m2tTP_1, m2tTP_2, m2tTP_3, m2tTP_4, m2tTP_5, m2lTP_0, m2lTP_1, m2lTP_2, m2lTP_3, m2lTP_4, m2lTP_5, tau, b1, b2);

		space_e(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, F1_alphac, F2_alphac, h_x, h_y, t_4, Mode);
		get_L2_projection_2D(EelTP_0, EelTP_1, EelTP_2, EelTP_3, EelTP_4, EelTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(Eel_0, Eel_1, Eel_2, Eel_3, Eel_4, Eel_5, Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, Eef_0, Eef_1, Eef_2, Eef_3, Eef_4, Eef_5, Ees_0, Ees_1, Ees_2, Ees_3, Ees_4, Ees_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, EesTP_0, EesTP_1, EesTP_2, EesTP_3, EesTP_4, EesTP_5, EetTP_0, EetTP_1, EetTP_2, EetTP_3, EetTP_4, EetTP_5, EelTP_0, EelTP_1, EelTP_2, EelTP_3, EelTP_4, EelTP_5, tau, b1, b2);

		space_i(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, F1_alphac, F2_alphac, h_x, h_y, t_4, Mode);
		get_L2_projection_2D(EilTP_0, EilTP_1, EilTP_2, EilTP_3, EilTP_4, EilTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(Eil_0, Eil_1, Eil_2, Eil_3, Eil_4, Eil_5, Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, Eif_0, Eif_1, Eif_2, Eif_3, Eif_4, Eif_5, Eis_0, Eis_1, Eis_2, Eis_3, Eis_4, Eis_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, EisTP_0, EisTP_1, EisTP_2, EisTP_3, EisTP_4, EisTP_5, EitTP_0, EitTP_1, EitTP_2, EitTP_3, EitTP_4, EitTP_5, EilTP_0, EilTP_1, EilTP_2, EilTP_3, EilTP_4, EilTP_5, tau, b1, b2);

		space_r(S_0, S_1, S_2, S_3, S_4, S_5, X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5, coef_L1, coef_R1, coef_L2, coef_R2, coef_L3, coef_R3, coef_L4, coef_R4, coef_L5, coef_R5, coef_B1, coef_T1, coef_B2, coef_T2, coef_B3, coef_T3, coef_B4, coef_T4, coef_B5, coef_T5, rhot_0, rhot_1, rhot_2, rhot_3, rhot_4, rhot_5, m1t_0, m1t_1, m1t_2, m1t_3, m1t_4, m1t_5, m2t_0, m2t_1, m2t_2, m2t_3, m2t_4, m2t_5, Eet_0, Eet_1, Eet_2, Eet_3, Eet_4, Eet_5, Eit_0, Eit_1, Eit_2, Eit_3, Eit_4, Eit_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, F1_alphac, F2_alphac, h_x, h_y, t_4, Mode);
		get_L2_projection_2D(ErlTP_0, ErlTP_1, ErlTP_2, ErlTP_3, ErlTP_4, ErlTP_5, S_0, S_1, S_2, S_3, S_4, S_5, h_x, h_y);
		IMEX3_l(Erl_0, Erl_1, Erl_2, Erl_3, Erl_4, Erl_5, Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, Erf_0, Erf_1, Erf_2, Erf_3, Erf_4, Erf_5, Ers_0, Ers_1, Ers_2, Ers_3, Ers_4, Ers_5, Ert_0, Ert_1, Ert_2, Ert_3, Ert_4, Ert_5, ErsTP_0, ErsTP_1, ErsTP_2, ErsTP_3, ErsTP_4, ErsTP_5, ErtTP_0, ErtTP_1, ErtTP_2, ErtTP_3, ErtTP_4, ErtTP_5, ErlTP_0, ErlTP_1, ErlTP_2, ErlTP_3, ErlTP_4, ErlTP_5, tau, b1, b2);

		// update boundary
		update_boundary(rhol_0, rhol_1, rhol_2, rhol_3, rhol_4, rhol_5); update_boundary(m1l_0, m1l_1, m1l_2, m1l_3, m1l_4, m1l_5);
		update_boundary(m2l_0, m2l_1, m2l_2, m2l_3, m2l_4, m2l_5); update_boundary(Eel_0, Eel_1, Eel_2, Eel_3, Eel_4, Eel_5);
		update_boundary(Eil_0, Eil_1, Eil_2, Eil_3, Eil_4, Eil_5); update_boundary(Erl_0, Erl_1, Erl_2, Erl_3, Erl_4, Erl_5);

		if (OE_index == 1)
		{
			OE_compute_system_sigma(sigma_0, sigma_1, sigma_2, rhol_0, rhol_1, rhol_2, rhol_3, rhol_4, rhol_5, m1l_0, m1l_1, m1l_2, m1l_3, m1l_4, m1l_5, m2l_0, m2l_1, m2l_2, m2l_3, m2l_4, m2l_5, Eel_0, Eel_1, Eel_2, Eel_3, Eel_4, Eel_5, Eil_0, Eil_1, Eil_2, Eil_3, Eil_4, Eil_5, Erl_0, Erl_1, Erl_2, Erl_3, Erl_4, Erl_5, h_x, h_y);

			OE_procedure(rhol_1, rhol_2, rhol_3, rhol_4, rhol_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m1l_1, m1l_2, m1l_3, m1l_4, m1l_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(m2l_1, m2l_2, m2l_3, m2l_4, m2l_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eel_1, Eel_2, Eel_3, Eel_4, Eel_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Eil_1, Eil_2, Eil_3, Eil_4, Eil_5, sigma_0, sigma_1, sigma_2, tau);
			OE_procedure(Erl_1, Erl_2, Erl_3, Erl_4, Erl_5, sigma_0, sigma_1, sigma_2, tau);

			// update boundary
			update_boundary(rhol_0, rhol_1, rhol_2, rhol_3, rhol_4, rhol_5); update_boundary(m1l_0, m1l_1, m1l_2, m1l_3, m1l_4, m1l_5);
			update_boundary(m2l_0, m2l_1, m2l_2, m2l_3, m2l_4, m2l_5); update_boundary(Eel_0, Eel_1, Eel_2, Eel_3, Eel_4, Eel_5);
			update_boundary(Eil_0, Eil_1, Eil_2, Eil_3, Eil_4, Eil_5); update_boundary(Erl_0, Erl_1, Erl_2, Erl_3, Erl_4, Erl_5);
		}

		deleteMatrix(N_y, sigma_0); deleteMatrix(N_y, sigma_1); deleteMatrix(N_y, sigma_2);

		deleteMatrix(N_y, coef_L1); deleteMatrix(N_y, coef_L2); deleteMatrix(N_y, coef_L3); deleteMatrix(N_y, coef_L4); deleteMatrix(N_y, coef_L5);
		deleteMatrix(N_y, coef_R1); deleteMatrix(N_y, coef_R2); deleteMatrix(N_y, coef_R3); deleteMatrix(N_y, coef_R4); deleteMatrix(N_y, coef_R5);
		deleteMatrix(N_y, coef_B1); deleteMatrix(N_y, coef_B2); deleteMatrix(N_y, coef_B3); deleteMatrix(N_y, coef_B4); deleteMatrix(N_y, coef_B5);
		deleteMatrix(N_y, coef_T1); deleteMatrix(N_y, coef_T2); deleteMatrix(N_y, coef_T3); deleteMatrix(N_y, coef_T4); deleteMatrix(N_y, coef_T5);

		deleteMatrix(N_y, F1_alphac1); deleteMatrix(N_y, F1_alphac2); deleteMatrix(N_y, F1_alphac3);
		deleteMatrix(N_y, F1_alphac4); deleteMatrix(N_y, F1_alphac5); deleteMatrix(N_y, F1_alphac);
		deleteMatrix(N_y + 1, F2_alphac1); deleteMatrix(N_y + 1, F2_alphac2); deleteMatrix(N_y + 1, F2_alphac3);
		deleteMatrix(N_y + 1, F2_alphac4); deleteMatrix(N_y + 1, F2_alphac5); deleteMatrix(N_y + 1, F2_alphac);
		deleteMatrix(mrows, S_0); deleteMatrix(mrows, S_1); deleteMatrix(mrows, S_2);
		deleteMatrix(mrows, S_3); deleteMatrix(mrows, S_4); deleteMatrix(mrows, S_5);

		deleteMatrix(mrows, rhofTP_0); deleteMatrix(mrows, rhofTP_1); deleteMatrix(mrows, rhofTP_2); deleteMatrix(mrows, rhofTP_3); deleteMatrix(mrows, rhofTP_4); deleteMatrix(mrows, rhofTP_5);
		deleteMatrix(mrows, rhosTP_0); deleteMatrix(mrows, rhosTP_1); deleteMatrix(mrows, rhosTP_2); deleteMatrix(mrows, rhosTP_3); deleteMatrix(mrows, rhosTP_4); deleteMatrix(mrows, rhosTP_5);
		deleteMatrix(mrows, rhotTP_0); deleteMatrix(mrows, rhotTP_1); deleteMatrix(mrows, rhotTP_2); deleteMatrix(mrows, rhotTP_3); deleteMatrix(mrows, rhotTP_4); deleteMatrix(mrows, rhotTP_5);
		deleteMatrix(mrows, rholTP_0); deleteMatrix(mrows, rholTP_1); deleteMatrix(mrows, rholTP_2); deleteMatrix(mrows, rholTP_3); deleteMatrix(mrows, rholTP_4); deleteMatrix(mrows, rholTP_5);

		deleteMatrix(mrows, m1fTP_0); deleteMatrix(mrows, m1fTP_1); deleteMatrix(mrows, m1fTP_2); deleteMatrix(mrows, m1fTP_3); deleteMatrix(mrows, m1fTP_4); deleteMatrix(mrows, m1fTP_5);
		deleteMatrix(mrows, m1sTP_0); deleteMatrix(mrows, m1sTP_1); deleteMatrix(mrows, m1sTP_2); deleteMatrix(mrows, m1sTP_3); deleteMatrix(mrows, m1sTP_4); deleteMatrix(mrows, m1sTP_5);
		deleteMatrix(mrows, m1tTP_0); deleteMatrix(mrows, m1tTP_1); deleteMatrix(mrows, m1tTP_2); deleteMatrix(mrows, m1tTP_3); deleteMatrix(mrows, m1tTP_4); deleteMatrix(mrows, m1tTP_5);
		deleteMatrix(mrows, m1lTP_0); deleteMatrix(mrows, m1lTP_1); deleteMatrix(mrows, m1lTP_2); deleteMatrix(mrows, m1lTP_3); deleteMatrix(mrows, m1lTP_4); deleteMatrix(mrows, m1lTP_5);

		deleteMatrix(mrows, m2fTP_0); deleteMatrix(mrows, m2fTP_1); deleteMatrix(mrows, m2fTP_2); deleteMatrix(mrows, m2fTP_3); deleteMatrix(mrows, m2fTP_4); deleteMatrix(mrows, m2fTP_5);
		deleteMatrix(mrows, m2sTP_0); deleteMatrix(mrows, m2sTP_1); deleteMatrix(mrows, m2sTP_2); deleteMatrix(mrows, m2sTP_3); deleteMatrix(mrows, m2sTP_4); deleteMatrix(mrows, m2sTP_5);
		deleteMatrix(mrows, m2tTP_0); deleteMatrix(mrows, m2tTP_1); deleteMatrix(mrows, m2tTP_2); deleteMatrix(mrows, m2tTP_3); deleteMatrix(mrows, m2tTP_4); deleteMatrix(mrows, m2tTP_5);
		deleteMatrix(mrows, m2lTP_0); deleteMatrix(mrows, m2lTP_1); deleteMatrix(mrows, m2lTP_2); deleteMatrix(mrows, m2lTP_3); deleteMatrix(mrows, m2lTP_4); deleteMatrix(mrows, m2lTP_5);

		deleteMatrix(mrows, EefTP_0); deleteMatrix(mrows, EefTP_1); deleteMatrix(mrows, EefTP_2); deleteMatrix(mrows, EefTP_3); deleteMatrix(mrows, EefTP_4); deleteMatrix(mrows, EefTP_5);
		deleteMatrix(mrows, EesTP_0); deleteMatrix(mrows, EesTP_1); deleteMatrix(mrows, EesTP_2); deleteMatrix(mrows, EesTP_3); deleteMatrix(mrows, EesTP_4); deleteMatrix(mrows, EesTP_5);
		deleteMatrix(mrows, EetTP_0); deleteMatrix(mrows, EetTP_1); deleteMatrix(mrows, EetTP_2); deleteMatrix(mrows, EetTP_3); deleteMatrix(mrows, EetTP_4); deleteMatrix(mrows, EetTP_5);
		deleteMatrix(mrows, EelTP_0); deleteMatrix(mrows, EelTP_1); deleteMatrix(mrows, EelTP_2); deleteMatrix(mrows, EelTP_3); deleteMatrix(mrows, EelTP_4); deleteMatrix(mrows, EelTP_5);

		deleteMatrix(mrows, EifTP_0); deleteMatrix(mrows, EifTP_1); deleteMatrix(mrows, EifTP_2); deleteMatrix(mrows, EifTP_3); deleteMatrix(mrows, EifTP_4); deleteMatrix(mrows, EifTP_5);
		deleteMatrix(mrows, EisTP_0); deleteMatrix(mrows, EisTP_1); deleteMatrix(mrows, EisTP_2); deleteMatrix(mrows, EisTP_3); deleteMatrix(mrows, EisTP_4); deleteMatrix(mrows, EisTP_5);
		deleteMatrix(mrows, EitTP_0); deleteMatrix(mrows, EitTP_1); deleteMatrix(mrows, EitTP_2); deleteMatrix(mrows, EitTP_3); deleteMatrix(mrows, EitTP_4); deleteMatrix(mrows, EitTP_5);
		deleteMatrix(mrows, EilTP_0); deleteMatrix(mrows, EilTP_1); deleteMatrix(mrows, EilTP_2); deleteMatrix(mrows, EilTP_3); deleteMatrix(mrows, EilTP_4); deleteMatrix(mrows, EilTP_5);

		deleteMatrix(mrows, ErfTP_0); deleteMatrix(mrows, ErfTP_1); deleteMatrix(mrows, ErfTP_2); deleteMatrix(mrows, ErfTP_3); deleteMatrix(mrows, ErfTP_4); deleteMatrix(mrows, ErfTP_5);
		deleteMatrix(mrows, ErsTP_0); deleteMatrix(mrows, ErsTP_1); deleteMatrix(mrows, ErsTP_2); deleteMatrix(mrows, ErsTP_3); deleteMatrix(mrows, ErsTP_4); deleteMatrix(mrows, ErsTP_5);
		deleteMatrix(mrows, ErtTP_0); deleteMatrix(mrows, ErtTP_1); deleteMatrix(mrows, ErtTP_2); deleteMatrix(mrows, ErtTP_3); deleteMatrix(mrows, ErtTP_4); deleteMatrix(mrows, ErtTP_5);
		deleteMatrix(mrows, ErlTP_0); deleteMatrix(mrows, ErlTP_1); deleteMatrix(mrows, ErlTP_2); deleteMatrix(mrows, ErlTP_3); deleteMatrix(mrows, ErlTP_4); deleteMatrix(mrows, ErlTP_5);

		// compute conservation
		double rho_whole = 0.0; double rhonew_whole = 0.0;
		for (int i = 1; i < N_y + 1; i++)
		{
			for (int j = 1; j < N_x + 1; j++)
			{
				rho_whole = rho_whole + rho_0[i][j];
				rhonew_whole = rhonew_whole + rhol_0[i][j];
			}

		}
		if (abs(rho_whole - rhonew_whole) > rho_c)
		{
			rho_c = abs(rho_whole - rhonew_whole);
		}

		double m1_whole = 0.0; double m1new_whole = 0.0;
		for (int i = 1; i < N_y + 1; i++)
		{
			for (int j = 1; j < N_x + 1; j++)
			{
				m1_whole = m1_whole + m1_0[i][j];
				m1new_whole = m1new_whole + m1l_0[i][j];
			}

		}
		if (abs(m1_whole - m1new_whole) > m1_c)
		{
			m1_c = abs(m1_whole - m1new_whole);
		}

		double m2_whole = 0.0; double m2new_whole = 0.0;
		for (int i = 1; i < N_y + 1; i++)
		{
			for (int j = 1; j < N_x + 1; j++)
			{
				m2_whole = m2_whole + m2_0[i][j];
				m2new_whole = m2new_whole + m2l_0[i][j];
			}

		}
		if (abs(m2_whole - m2new_whole) > m2_c)
		{
			m2_c = abs(m2_whole - m2new_whole);
		}

		for (int i = 0; i < N_y + 2; i++)
		{
			for (int j = 0; j < N_x + 2; j++)
			{
				rho_0[i][j] = rhol_0[i][j]; rho_1[i][j] = rhol_1[i][j]; rho_2[i][j] = rhol_2[i][j];
				rho_3[i][j] = rhol_3[i][j]; rho_4[i][j] = rhol_4[i][j]; rho_5[i][j] = rhol_5[i][j];

				m1_0[i][j] = m1l_0[i][j]; m1_1[i][j] = m1l_1[i][j]; m1_2[i][j] = m1l_2[i][j];
				m1_3[i][j] = m1l_3[i][j]; m1_4[i][j] = m1l_4[i][j]; m1_5[i][j] = m1l_5[i][j];

				m2_0[i][j] = m2l_0[i][j]; m2_1[i][j] = m2l_1[i][j]; m2_2[i][j] = m2l_2[i][j];
				m2_3[i][j] = m2l_3[i][j]; m2_4[i][j] = m2l_4[i][j]; m2_5[i][j] = m2l_5[i][j];

				Ee_0[i][j] = Eel_0[i][j]; Ee_1[i][j] = Eel_1[i][j]; Ee_2[i][j] = Eel_2[i][j];
				Ee_3[i][j] = Eel_3[i][j]; Ee_4[i][j] = Eel_4[i][j]; Ee_5[i][j] = Eel_5[i][j];

				Ei_0[i][j] = Eil_0[i][j]; Ei_1[i][j] = Eil_1[i][j]; Ei_2[i][j] = Eil_2[i][j];
				Ei_3[i][j] = Eil_3[i][j]; Ei_4[i][j] = Eil_4[i][j]; Ei_5[i][j] = Eil_5[i][j];

				Er_0[i][j] = Erl_0[i][j]; Er_1[i][j] = Erl_1[i][j]; Er_2[i][j] = Erl_2[i][j];
				Er_3[i][j] = Erl_3[i][j]; Er_4[i][j] = Erl_4[i][j]; Er_5[i][j] = Erl_5[i][j];

				E_new[i][j] = Ee_0[i][j] + Ei_0[i][j] + Er_0[i][j];
			}
		}

		double E_whole = 0.0; double Enew_whole = 0.0;
		for (int i = 1; i < N_y + 1; i++)
		{
			for (int j = 1; j < N_x + 1; j++)
			{
				E_whole = E_whole + E[i][j];
				Enew_whole = Enew_whole + E_new[i][j];
			}

		}
		if (abs(E_whole - Enew_whole) > E_c)
		{
			E_c = abs(E_whole - Enew_whole);
		}

		double rho_MAX = compute_MAX(rhol_0);
		double m1_MAX = compute_MAX(m1l_0);
		double m2_MAX = compute_MAX(m2l_0);
		double Ee_MAX = compute_MAX(Eel_0);
		double Ei_MAX = compute_MAX(Eil_0);
		double Er_MAX = compute_MAX(Erl_0);

		if (rho_MAX == 0.0 || m1_MAX == 0.0 || m2_MAX == 0.0 || Ee_MAX == 0.0 || Ei_MAX == 0.0 || Er_MAX == 0.0)
		{
			cout << "Blow up at t = " << t << endl;
			break;
		}

		t = t + tau;

		deleteMatrix(mrows, rhof_0); deleteMatrix(mrows, rhof_1); deleteMatrix(mrows, rhof_2); deleteMatrix(mrows, rhof_3); deleteMatrix(mrows, rhof_4); deleteMatrix(mrows, rhof_5);
		deleteMatrix(mrows, rhos_0); deleteMatrix(mrows, rhos_1); deleteMatrix(mrows, rhos_2); deleteMatrix(mrows, rhos_3); deleteMatrix(mrows, rhos_4); deleteMatrix(mrows, rhos_5);
		deleteMatrix(mrows, rhot_0); deleteMatrix(mrows, rhot_1); deleteMatrix(mrows, rhot_2); deleteMatrix(mrows, rhot_3); deleteMatrix(mrows, rhot_4); deleteMatrix(mrows, rhot_5);
		deleteMatrix(mrows, rhol_0); deleteMatrix(mrows, rhol_1); deleteMatrix(mrows, rhol_2); deleteMatrix(mrows, rhol_3); deleteMatrix(mrows, rhol_4); deleteMatrix(mrows, rhol_5);

		deleteMatrix(mrows, m1f_0); deleteMatrix(mrows, m1f_1); deleteMatrix(mrows, m1f_2); deleteMatrix(mrows, m1f_3); deleteMatrix(mrows, m1f_4); deleteMatrix(mrows, m1f_5);
		deleteMatrix(mrows, m1s_0); deleteMatrix(mrows, m1s_1); deleteMatrix(mrows, m1s_2); deleteMatrix(mrows, m1s_3); deleteMatrix(mrows, m1s_4); deleteMatrix(mrows, m1s_5);
		deleteMatrix(mrows, m1t_0); deleteMatrix(mrows, m1t_1); deleteMatrix(mrows, m1t_2); deleteMatrix(mrows, m1t_3); deleteMatrix(mrows, m1t_4); deleteMatrix(mrows, m1t_5);
		deleteMatrix(mrows, m1l_0); deleteMatrix(mrows, m1l_1); deleteMatrix(mrows, m1l_2); deleteMatrix(mrows, m1l_3); deleteMatrix(mrows, m1l_4); deleteMatrix(mrows, m1l_5);

		deleteMatrix(mrows, m2f_0); deleteMatrix(mrows, m2f_1); deleteMatrix(mrows, m2f_2); deleteMatrix(mrows, m2f_3); deleteMatrix(mrows, m2f_4); deleteMatrix(mrows, m2f_5);
		deleteMatrix(mrows, m2s_0); deleteMatrix(mrows, m2s_1); deleteMatrix(mrows, m2s_2); deleteMatrix(mrows, m2s_3); deleteMatrix(mrows, m2s_4); deleteMatrix(mrows, m2s_5);
		deleteMatrix(mrows, m2t_0); deleteMatrix(mrows, m2t_1); deleteMatrix(mrows, m2t_2); deleteMatrix(mrows, m2t_3); deleteMatrix(mrows, m2t_4); deleteMatrix(mrows, m2t_5);
		deleteMatrix(mrows, m2l_0); deleteMatrix(mrows, m2l_1); deleteMatrix(mrows, m2l_2); deleteMatrix(mrows, m2l_3); deleteMatrix(mrows, m2l_4); deleteMatrix(mrows, m2l_5);

		deleteMatrix(mrows, Eef_0); deleteMatrix(mrows, Eef_1); deleteMatrix(mrows, Eef_2); deleteMatrix(mrows, Eef_3); deleteMatrix(mrows, Eef_4); deleteMatrix(mrows, Eef_5);
		deleteMatrix(mrows, Ees_0); deleteMatrix(mrows, Ees_1); deleteMatrix(mrows, Ees_2); deleteMatrix(mrows, Ees_3); deleteMatrix(mrows, Ees_4); deleteMatrix(mrows, Ees_5);
		deleteMatrix(mrows, Eet_0); deleteMatrix(mrows, Eet_1); deleteMatrix(mrows, Eet_2); deleteMatrix(mrows, Eet_3); deleteMatrix(mrows, Eet_4); deleteMatrix(mrows, Eet_5);
		deleteMatrix(mrows, Eel_0); deleteMatrix(mrows, Eel_1); deleteMatrix(mrows, Eel_2); deleteMatrix(mrows, Eel_3); deleteMatrix(mrows, Eel_4); deleteMatrix(mrows, Eel_5);

		deleteMatrix(mrows, Eif_0); deleteMatrix(mrows, Eif_1); deleteMatrix(mrows, Eif_2); deleteMatrix(mrows, Eif_3); deleteMatrix(mrows, Eif_4); deleteMatrix(mrows, Eif_5);
		deleteMatrix(mrows, Eis_0); deleteMatrix(mrows, Eis_1); deleteMatrix(mrows, Eis_2); deleteMatrix(mrows, Eis_3); deleteMatrix(mrows, Eis_4); deleteMatrix(mrows, Eis_5);
		deleteMatrix(mrows, Eit_0); deleteMatrix(mrows, Eit_1); deleteMatrix(mrows, Eit_2); deleteMatrix(mrows, Eit_3); deleteMatrix(mrows, Eit_4); deleteMatrix(mrows, Eit_5);
		deleteMatrix(mrows, Eil_0); deleteMatrix(mrows, Eil_1); deleteMatrix(mrows, Eil_2); deleteMatrix(mrows, Eil_3); deleteMatrix(mrows, Eil_4); deleteMatrix(mrows, Eil_5);

		deleteMatrix(mrows, Erf_0); deleteMatrix(mrows, Erf_1); deleteMatrix(mrows, Erf_2); deleteMatrix(mrows, Erf_3); deleteMatrix(mrows, Erf_4); deleteMatrix(mrows, Erf_5);
		deleteMatrix(mrows, Ers_0); deleteMatrix(mrows, Ers_1); deleteMatrix(mrows, Ers_2); deleteMatrix(mrows, Ers_3); deleteMatrix(mrows, Ers_4); deleteMatrix(mrows, Ers_5);
		deleteMatrix(mrows, Ert_0); deleteMatrix(mrows, Ert_1); deleteMatrix(mrows, Ert_2); deleteMatrix(mrows, Ert_3); deleteMatrix(mrows, Ert_4); deleteMatrix(mrows, Ert_5);
		deleteMatrix(mrows, Erl_0); deleteMatrix(mrows, Erl_1); deleteMatrix(mrows, Erl_2); deleteMatrix(mrows, Erl_3); deleteMatrix(mrows, Erl_4); deleteMatrix(mrows, Erl_5);

	}

	cout << tau << endl;
    
	delete[] X_1; delete[] X_2; delete[] X_3; delete[] X_4; delete[] X_5;
	delete[] Y_1; delete[] Y_2; delete[] Y_3; delete[] Y_4; delete[] Y_5; deleteMatrix(mrows, E); deleteMatrix(mrows, E_new);
	deleteMatrix(N_y, s_j); deleteMatrix(N_y, beta_ij_x); deleteMatrix(N_y, beta_ij_y);

	if (t == t_fin)
	{
		cout << "The computing of numerical solution is now finished ! \n" << endl;
	}

	// whether conservative
	cout << "the case of rho: " << rho_c << endl;
	cout << "the case of m: " << m1_c << endl;
	cout << "the case of m: " << m2_c << endl;
	cout << "the case of E: " << E_c << endl;

	// compute error
	cout << "Start computing the error ..." << endl;
	//error of rho
	double* rho_ev = rho_err_2D(rho_0, rho_1, rho_2, rho_3, rho_4, rho_5, x, y, h_x, h_y, t);
	cout << "L1 error of rho under mesh " << N_x << "¡Á" << N_y << " = " << rho_ev[0] << endl;
	cout << "L2 error of rho under mesh " << N_x << "¡Á" << N_y << " = " << rho_ev[1] << endl;
	cout << "Li error of rho under mesh " << N_x << "¡Á" << N_y << " = " << rho_ev[2] << "\n" << endl;
	deleteMatrix(mrows, rho_0); deleteMatrix(mrows, rho_1); deleteMatrix(mrows, rho_2);
	deleteMatrix(mrows, rho_3); deleteMatrix(mrows, rho_4); deleteMatrix(mrows, rho_5); delete[] rho_ev;

	//error of m1
	double* m1_ev = m1_err_2D(m1_0, m1_1, m1_2, m1_3, m1_4, m1_5, x, y, h_x, h_y, t);
	cout << "L1 error of m1 under mesh " << N_x << "¡Á" << N_y << " = " << m1_ev[0] << endl;
	cout << "L2 error of m1 under mesh " << N_x << "¡Á" << N_y << " = " << m1_ev[1] << endl;
	cout << "Li error of m1 under mesh " << N_x << "¡Á" << N_y << " = " << m1_ev[2] << "\n" << endl;
	deleteMatrix(mrows, m1_0); deleteMatrix(mrows, m1_1); deleteMatrix(mrows, m1_2);
	deleteMatrix(mrows, m1_3); deleteMatrix(mrows, m1_4); deleteMatrix(mrows, m1_5); delete[] m1_ev;

	//error of m2
	double* m2_ev = m2_err_2D(m2_0, m2_1, m2_2, m2_3, m2_4, m2_5, x, y, h_x, h_y, t);
	cout << "L1 error of m2 under mesh " << N_x << "¡Á" << N_y << " = " << m2_ev[0] << endl;
	cout << "L2 error of m2 under mesh " << N_x << "¡Á" << N_y << " = " << m2_ev[1] << endl;
	cout << "Li error of m2 under mesh " << N_x << "¡Á" << N_y << " = " << m2_ev[2] << "\n" << endl;
	deleteMatrix(mrows, m2_0); deleteMatrix(mrows, m2_1); deleteMatrix(mrows, m2_2);
	deleteMatrix(mrows, m2_3); deleteMatrix(mrows, m2_4); deleteMatrix(mrows, m2_5); delete[] m2_ev;

	//error of Ee
	double* Ee_ev = Ee_err_2D(Ee_0, Ee_1, Ee_2, Ee_3, Ee_4, Ee_5, x, y, h_x, h_y, t);
	cout << "L1 error of Ee under mesh " << N_x << "¡Á" << N_y << " = " << Ee_ev[0] << endl;
	cout << "L2 error of Ee under mesh " << N_x << "¡Á" << N_y << " = " << Ee_ev[1] << endl;
	cout << "Li error of Ee under mesh " << N_x << "¡Á" << N_y << " = " << Ee_ev[2] << "\n" << endl;
	deleteMatrix(mrows, Ee_0); deleteMatrix(mrows, Ee_1); deleteMatrix(mrows, Ee_2);
	deleteMatrix(mrows, Ee_3); deleteMatrix(mrows, Ee_4); deleteMatrix(mrows, Ee_5); delete[] Ee_ev;

	//error of Ei
	double* Ei_ev = Ei_err_2D(Ei_0, Ei_1, Ei_2, Ei_3, Ei_4, Ei_5, x, y, h_x, h_y, t);
	cout << "L1 error of Ei under mesh " << N_x << "¡Á" << N_y << " = " << Ei_ev[0] << endl;
	cout << "L2 error of Ei under mesh " << N_x << "¡Á" << N_y << " = " << Ei_ev[1] << endl;
	cout << "Li error of Ei under mesh " << N_x << "¡Á" << N_y << " = " << Ei_ev[2] << "\n" << endl;
	deleteMatrix(mrows, Ei_0); deleteMatrix(mrows, Ei_1); deleteMatrix(mrows, Ei_2);
	deleteMatrix(mrows, Ei_3); deleteMatrix(mrows, Ei_4); deleteMatrix(mrows, Ei_5); delete[] Ei_ev;

	//error of Er
	double* Er_ev = Er_err_2D(Er_0, Er_1, Er_2, Er_3, Er_4, Er_5, x, y, h_x, h_y, t);
	cout << "L1 error of Er under mesh " << N_x << "¡Á" << N_y << " = " << Er_ev[0] << endl;
	cout << "L2 error of Er under mesh " << N_x << "¡Á" << N_y << " = " << Er_ev[1] << endl;
	cout << "Li error of Er under mesh " << N_x << "¡Á" << N_y << " = " << Er_ev[2] << "\n" << endl;
	deleteMatrix(mrows, Er_0); deleteMatrix(mrows, Er_1); deleteMatrix(mrows, Er_2);
	deleteMatrix(mrows, Er_3); deleteMatrix(mrows, Er_4); deleteMatrix(mrows, Er_5); delete[] Er_ev;

	delete[] x; delete[] y;
	cout << "The computing of error is now finished ! \n" << endl;
	
	// timing ended
	end = std::time(NULL);
	cout << "Time has passed " << end - start << " seconds." << endl;
	return 0;
}
