#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG.h"

// Modal DG general framework
void get_Gaussian_points(double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* x, int N_index, double h_x)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		X_1[i] = x[i] + h_x / 2.0 * GQxi_1; X_2[i] = x[i] + h_x / 2.0 * GQxi_2; X_3[i] = x[i] + h_x / 2.0 * GQxi_3;
		X_4[i] = x[i] + h_x / 2.0 * GQxi_4; X_5[i] = x[i] + h_x / 2.0 * GQxi_5;
	}
}

void get_L2_projection(double* U_0, double* U_1, double* U_2, double* u_0, double* u_1, double* u_2, int N_index, double h_x)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = u_0[i] / h_x;
		U_1[i] = 3.0 * u_1[i] / h_x;
		U_2[i] = 5.0 * u_2[i] / h_x;
	}
}

void get_Gauss_quadrature(double* U_0, double* U_1, double* U_2, double* u_x1, double* u_x2, double* u_x3, double* u_x4, double* u_x5, int N_index, double h_x)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = h_x / 2.0 * (weight_1 * u_x1[i] + weight_2 * u_x2[i] + weight_3 * u_x3[i] + weight_4 * u_x4[i] + weight_5 * u_x5[i]);
		U_1[i] = h_x / 2.0 * (weight_1 * u_x1[i] * GQxi_1 + weight_2 * u_x2[i] * GQxi_2 + weight_3 * u_x3[i] * GQxi_3 + weight_4 * u_x4[i] * GQxi_4 + weight_5 * u_x5[i] * GQxi_5);
		U_2[i] = h_x / 2.0 * (weight_1 * u_x1[i] * (3.0 * pow(GQxi_1, 2) - 1.0) / 2.0 + weight_2 * u_x2[i] * (3.0 * pow(GQxi_2, 2) - 1.0) / 2.0 + weight_3 * u_x3[i] * (3.0 * pow(GQxi_3, 2) - 1.0) / 2.0 + weight_4 * u_x4[i] * (3.0 * pow(GQxi_4, 2) - 1.0) / 2.0 + weight_5 * u_x5[i] * (3.0 * pow(GQxi_5, 2) - 1.0) / 2.0);
	}
}

void get_Gauss_quadrature_prime(double* U_0, double* U_1, double* U_2, double* u_0, double* u_1, int N_index, double h_x)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = 0;
		U_1[i] = 2.0 * u_0[i] / h_x;
		U_2[i] = 6.0 * u_1[i] / h_x;
	}
}

void get_time_prime(double* U_0, double* U_1, double* U_2, double* S_0, double* S_1, double* S_2, int N_index, double h_x)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		U_0[i] = S_0[i] / h_x;
		U_1[i] = 3.0 * S_1[i] / h_x;
		U_2[i] = 5.0 * S_2[i] / h_x;
	}
}

void Mode_select(double* u_0, double* u_1, double* u_2, int index_mode, int N_index)
{
	int N = N_index;
	if (index_mode == 1)
	{
		for (int i = 0; i < N; i++)
		{
			u_2[i] = 0;
		}
	}
	if (index_mode == 0)
	{
		for (int i = 0; i < N; i++)
		{
			u_1[i] = 0;
			u_2[i] = 0;
		}
	}
}

void update_boundary(double* u_0, double* u_1, double* u_2)
{
	// periodic
	// left boundary
	u_0[0] = u_0[N_x];
	u_1[0] = u_1[N_x];
	u_2[0] = u_2[N_x];
	// right boundary
	u_0[N_x + 1] = u_0[1];
	u_1[N_x + 1] = u_1[1];
	u_2[N_x + 1] = u_2[1];
}


// accuracy solution
double rho_ex(double x, double t)
{
	return 1.0 + 0.1 * sin(x - t) * sin(x - t);
}

double u_ex(double x, double t)
{
	return 1.0;
}

double rhoe_e_ex(double x, double t)
{
	return 3.0;
}

double rhoe_i_ex(double x, double t)
{
	return 3.0;
}

double rhoe_r_ex(double x, double t)
{
	return 2.0;
}

double rho_prime_ex(double x, double t)
{
	return 0.2 * sin(x - t) * cos(x - t);
}

double rho_prime2_ex(double x, double t)
{
	return 0.2 * (cos(x - t) * cos(x - t) - sin(x - t) * sin(x - t));
}


// model solution
double m_ex(double x, double t)
{
	return rho_ex(x, t) * u_ex(x, t);
}

double Ee_ex(double x, double t)
{
	return rhoe_e_ex(x, t) + rho_ex(x, t) * u_ex(x, t) * u_ex(x, t) / 6.0;
}

double Ei_ex(double x, double t)
{
	return rhoe_i_ex(x, t) + rho_ex(x, t) * u_ex(x, t) * u_ex(x, t) / 6.0;
}

double Er_ex(double x, double t)
{
	return rhoe_r_ex(x, t) + rho_ex(x, t) * u_ex(x, t) * u_ex(x, t) / 6.0;
}

double S3_ex(double x, double t)
{
	double Tr4 = rhoe_r_ex(x, t) / a;
	double Te = rhoe_e_ex(x, t) / (C_ve * rho_ex(x, t));
	double Ti = rhoe_i_ex(x, t) / (C_vi * rho_ex(x, t));
	double Source_i = omega_ei * (Te - Ti);
	double Source_r = omega_er * (pow(Te, 4) - Tr4);
	double Source_e = Source_i + Source_r;
	double D1 = 2.0 * pow(rho_prime_ex(x, t), 2) / pow(rho_ex(x, t), 3) - rho_prime2_ex(x, t) / pow(rho_ex(x, t), 2);
	double Diff_e = kappa_e / C_ve * D1 * rhoe_e_ex(x, t);
	double S3 = Source_e - Diff_e;
	return S3;
}

double S4_ex(double x, double t)
{
	double Te = rhoe_e_ex(x, t) / (C_ve * rho_ex(x, t));
	double Ti = rhoe_i_ex(x, t) / (C_vi * rho_ex(x, t));
	double Source_i = omega_ei * (Te - Ti);
	double D1 = 2.0 * pow(rho_prime_ex(x, t), 2) / pow(rho_ex(x, t), 3) - rho_prime2_ex(x, t) / pow(rho_ex(x, t), 2);
	double Diff_i = kappa_i / C_vi * D1 * rhoe_i_ex(x, t);
	double S4 = -1.0 * Source_i - Diff_i;
	return S4;
}

double S5_ex(double x, double t)
{
	double Tr4 = rhoe_r_ex(x, t) / a;
	double Te = rhoe_e_ex(x, t) / (C_ve * rho_ex(x, t));
	double Source_r = omega_er * (pow(Te, 4) - Tr4);
	double S5 = -1.0 * Source_r;
	return S5;
}


// Compute error and order
double* rho_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		int n = 9; // the number of points selected in each cell
		double err = 0.0;
		for (int k = 0; k < n; k++)
		{
			switch (Mode)
			{
			case 0:
				err = abs(rho_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i]));
				break;
			case 1:
				err = abs(rho_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0)));
				break;
			case 2:
				err = abs(rho_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0) + u_2[i] * (3 * pow((2.0 / (n - 1.0) * k - 1.0), 2) - 1.0) / 2.0));
				break;
			default:
				cout << "please choose the correct mode !" << endl;
				break;
			}
			L1_error = L1_error + err / n;
			L2_error = L2_error + pow(err, 2) / n;
			if (err > Linf_error)
			{
				Linf_error = err;
			}
		}
	}
	L1_error = L1_error / N_x;
	L2_error = sqrt(L2_error / N_x);
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* m_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		int n = 9; // the number of points selected in each cell
		double err = 0.0;
		for (int k = 0; k < n; k++)
		{
			switch (Mode)
			{
			case 0:
				err = abs(m_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i]));
				break;
			case 1:
				err = abs(m_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0)));
				break;
			case 2:
				err = abs(m_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0) + u_2[i] * (3 * pow((2.0 / (n - 1.0) * k - 1.0), 2) - 1.0) / 2.0));
				break;
			default:
				cout << "please choose the correct mode !" << endl;
				break;
			}
			L1_error = L1_error + err / n;
			L2_error = L2_error + pow(err, 2) / n;
			if (err > Linf_error)
			{
				Linf_error = err;
			}
		}
	}
	L1_error = L1_error / N_x;
	L2_error = sqrt(L2_error / N_x);
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Ee_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		int n = 9; // the number of points selected in each cell
		double err = 0.0;
		for (int k = 0; k < n; k++)
		{
			switch (Mode)
			{
			case 0:
				err = abs(Ee_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i]));
				break;
			case 1:
				err = abs(Ee_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0)));
				break;
			case 2:
				err = abs(Ee_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0) + u_2[i] * (3 * pow((2.0 / (n - 1.0) * k - 1.0), 2) - 1.0) / 2.0));
				break;
			default:
				cout << "please choose the correct mode !" << endl;
				break;
			}
			L1_error = L1_error + err / n;
			L2_error = L2_error + pow(err, 2) / n;
			if (err > Linf_error)
			{
				Linf_error = err;
			}
		}
	}
	L1_error = L1_error / N_x;
	L2_error = sqrt(L2_error / N_x);
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Ei_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		int n = 9; // the number of points selected in each cell
		double err = 0.0;
		for (int k = 0; k < n; k++)
		{
			switch (Mode)
			{
			case 0:
				err = abs(Ei_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i]));
				break;
			case 1:
				err = abs(Ei_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0)));
				break;
			case 2:
				err = abs(Ei_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0) + u_2[i] * (3 * pow((2.0 / (n - 1.0) * k - 1.0), 2) - 1.0) / 2.0));
				break;
			default:
				cout << "please choose the correct mode !" << endl;
				break;
			}
			L1_error = L1_error + err / n;
			L2_error = L2_error + pow(err, 2) / n;
			if (err > Linf_error)
			{
				Linf_error = err;
			}
		}
	}
	L1_error = L1_error / N_x;
	L2_error = sqrt(L2_error / N_x);
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}

double* Er_err(double* u_0, double* u_1, double* u_2, double* x, double h_x, double t)
{
	double L1_error = 0.0, L2_error = 0.0, Linf_error = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		int n = 9; // the number of points selected in each cell
		double err = 0.0;
		for (int k = 0; k < n; k++)
		{
			switch (Mode)
			{
			case 0:
				err = abs(Er_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i]));
				break;
			case 1:
				err = abs(Er_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0)));
				break;
			case 2:
				err = abs(Er_ex(x[i] + h_x / 2.0 * (2.0 / (n - 1.0) * k - 1.0), t) - (u_0[i] + u_1[i] * (2.0 / (n - 1.0) * k - 1.0) + u_2[i] * (3 * pow((2.0 / (n - 1.0) * k - 1.0), 2) - 1.0) / 2.0));
				break;
			default:
				cout << "please choose the correct mode !" << endl;
				break;
			}
			L1_error = L1_error + err / n;
			L2_error = L2_error + pow(err, 2) / n;
			if (err > Linf_error)
			{
				Linf_error = err;
			}
		}
	}
	L1_error = L1_error / N_x;
	L2_error = sqrt(L2_error / N_x);
	double* err_vector = new double[3];
	err_vector[0] = L1_error; err_vector[1] = L2_error; err_vector[2] = Linf_error;
	return err_vector;
}
