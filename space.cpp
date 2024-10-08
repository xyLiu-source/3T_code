#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG.h"
#include "space.h"

void F_d(double* u, double* u_0, double* u_1, double* u_2, double xi, int N_index)
{
	int N = N_index;
	for (int i = 0; i < N; i++)
	{
		u[i] = u_0[i] + u_1[i] * xi + u_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
	}
}

void C_s(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ee = Ee_0[i] + Ee_1[i] * xi + Ee_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_Ei = Ei_0[i] + Ei_1[i] * xi + Ei_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_Er = Er_0[i] + Er_1[i] * xi + Er_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double e_e = (temp_Ee - temp_rhou2 / 6.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double e_i = (temp_Ei - temp_rhou2 / 6.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double e_r = (temp_Er - temp_rhou2 / 6.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = sqrt(gamma_e * (gamma_e - 1.0) * e_e + gamma_i * (gamma_i - 1.0) * e_i + gamma_r * (gamma_r - 1.0) * e_r);
	}
}

void compute_alphaC(double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* m_plus = new double[N_x + 2]; double* m_minus = new double[N_x + 2];
	double* rho_plus = new double[N_x + 2]; double* rho_minus = new double[N_x + 2];
	F_d(m_plus, m_0, m_1, m_2, -1.0, N_x + 2); F_d(m_minus, m_0, m_1, m_2, 1.0, N_x + 2);
	F_d(rho_plus, rho_0, rho_1, rho_2, -1.0, N_x + 2); F_d(rho_minus, rho_0, rho_1, rho_2, 1.0, N_x + 2);
	double* u_plus = new double[N_x + 2]; double* u_minus = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		u_plus[i] = m_plus[i] / rho_plus[i];
		u_minus[i] = m_minus[i] / rho_minus[i];
	}
	double* Cs_plus = new double[N_x + 2]; double* Cs_minus = new double[N_x + 2];
	C_s(Cs_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, -1.0);
	C_s(Cs_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, 1.0);
	for (int i = 1; i < N_x + 2; i++)
	{
		double temp_max = abs(u_plus[i] - Cs_plus[i]);
		if (abs(u_plus[i]) > temp_max)
		{
			temp_max = abs(u_plus[i]);
		}
		if (abs(u_plus[i] + Cs_plus[i]) > temp_max)
		{
			temp_max = abs(u_plus[i] + Cs_plus[i]);
		}
		if (abs(u_minus[i - 1] - Cs_minus[i - 1]) > temp_max)
		{
			temp_max = abs(u_minus[i - 1] - Cs_minus[i - 1]);
		}
		if (abs(u_minus[i - 1]) > temp_max)
		{
			temp_max = abs(u_minus[i - 1]);
		}
		if (abs(u_minus[i - 1] + Cs_minus[i - 1]) > temp_max)
		{
			temp_max = abs(u_minus[i - 1] + Cs_minus[i - 1]);
		}
		alpha_c[i - 1] = temp_max;
	}

	delete[] m_plus; delete[] m_minus;
	delete[] rho_plus; delete[] rho_minus;
	delete[] u_plus; delete[] u_minus;
	delete[] Cs_plus; delete[] Cs_minus;
}

void compute_uhalf(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x, int index_mode)
{
	double* m1 = new double[N_x + 2]; double* m2 = new double[N_x + 2]; double* m3 = new double[N_x + 2];
	double* m4 = new double[N_x + 2]; double* m5 = new double[N_x + 2];
	F_d(m1, m_0, m_1, m_2, GQxi_1, N_x + 2); F_d(m2, m_0, m_1, m_2, GQxi_2, N_x + 2); F_d(m3, m_0, m_1, m_2, GQxi_3, N_x + 2);
	F_d(m4, m_0, m_1, m_2, GQxi_4, N_x + 2); F_d(m5, m_0, m_1, m_2, GQxi_5, N_x + 2);
	double* rho1 = new double[N_x + 2]; double* rho2 = new double[N_x + 2]; double* rho3 = new double[N_x + 2];
	double* rho4 = new double[N_x + 2]; double* rho5 = new double[N_x + 2];
	F_d(rho1, rho_0, rho_1, rho_2, GQxi_1, N_x + 2); F_d(rho2, rho_0, rho_1, rho_2, GQxi_2, N_x + 2); F_d(rho3, rho_0, rho_1, rho_2, GQxi_3, N_x + 2);
	F_d(rho4, rho_0, rho_1, rho_2, GQxi_4, N_x + 2); F_d(rho5, rho_0, rho_1, rho_2, GQxi_5, N_x + 2);
	double* u_1 = new double[N_x + 2]; double* u_2 = new double[N_x + 2]; double* u_3 = new double[N_x + 2]; double* u_4 = new double[N_x + 2]; double* u_5 = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		u_1[i] = m1[i] / rho1[i]; u_2[i] = m2[i] / rho2[i]; u_3[i] = m3[i] / rho3[i];
		u_4[i] = m4[i] / rho4[i]; u_5[i] = m5[i] / rho5[i];
	}
	double* u_GQ_0 = new double[N_x + 2]; double* u_GQ_1 = new double[N_x + 2]; double* u_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(u_GQ_0, u_GQ_1, u_GQ_2, u_1, u_2, u_3, u_4, u_5, N_x + 2, h_x);
	double* omega_0 = new double[N_x + 2]; double* omega_1 = new double[N_x + 2]; double* omega_2 = new double[N_x + 2];
	get_L2_projection(omega_0, omega_1, omega_2, u_GQ_0, u_GQ_1, u_GQ_2, N_x + 2, h_x);
	if (index_mode == 0)
	{
		for (int i = 0; i < N_x + 1; i++)
		{
			u[i] = (omega_0[i] + omega_0[i + 1]) / 2.0;
		}
	}
	if (index_mode == 1)
	{
		for (int i = 0; i < N_x + 1; i++)
		{
			u[i] = omega_0[i] + omega_1[i] + (omega_1[i + 1] - omega_1[i]) / 6.0 - (omega_0[i] + omega_1[i] - omega_0[i + 1] + omega_1[i + 1]) / 2.0;
		}
	}
	if (index_mode == 2)
	{
		double* w_3 = new double[N_x + 1]; double* w_4 = new double[N_x + 1]; double* w_5 = new double[N_x + 1];
		for (int i = 0; i < N_x + 1; i++)
		{
			w_3[i] = 7.0 * (omega_1[i + 1] - 3 * omega_2[i + 1] - omega_1[i] - 3 * omega_2[i]) / 64.0 + 5.0 * (omega_0[i + 1] - omega_1[i + 1] + omega_2[i + 1] - omega_0[i] - omega_1[i] - omega_2[i]) / 6.0 - (omega_2[i + 1] - omega_2[i]) / 15.0;
			w_4[i] = -3.0 * (omega_0[i + 1] - omega_1[i + 1] + omega_2[i + 1] - omega_0[i] - omega_1[i] - omega_2[i]) / 8.0 + 3.0 * (omega_2[i + 1] - omega_2[i]) / 40.0 - (omega_1[i + 1] - 3 * omega_2[i + 1] - omega_1[i] - 3 * omega_2[i]) / 64.0;
			w_5[i] = (omega_0[i + 1] - omega_1[i + 1] + omega_2[i + 1] - omega_0[i] - omega_1[i] - omega_2[i] - omega_2[i + 1] / 5.0 + omega_2[i] / 5.0) / 24.0;
		}

		for (int i = 0; i < N_x + 1; i++)
		{
			u[i] = omega_0[i] + omega_1[i] + omega_2[i] + w_3[i] + w_4[i] + w_5[i];
		}
		delete[] w_3; delete[] w_4; delete[] w_5;
	}

	delete[] m1; delete[] m2; delete[] m3; delete[] m4; delete[] m5;
	delete[] rho1; delete[] rho2; delete[] rho3; delete[] rho4; delete[] rho5;
	delete[] u_1; delete[] u_2; delete[] u_3; delete[] u_4; delete[] u_5;
	delete[] u_GQ_0; delete[] u_GQ_1; delete[] u_GQ_2;
	delete[] omega_0; delete[] omega_1; delete[] omega_2;
}

void compute_coef(double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x, int index_mode)
{
	double* u_half = new double[N_x + 1];
	compute_uhalf(u_half, rho_0, rho_1, rho_2, m_0, m_1, m_2, h_x, index_mode);
	for (int i = 0; i < N_x; i++)
	{
		double temp_L = -u_half[i] / 3.0;
		if (temp_L < 0)
		{
			temp_L = 0;

		}
		coef_L[i] = temp_L;

		double temp_R = -u_half[i + 1] / 3.0;
		if (temp_R > 0)
		{
			temp_R = 0;
		}
		coef_R[i] = temp_R;
	}
	delete[] u_half;
}

void compute_jump(double* u_jump, double* u_plus, double* u_minus, int N_index)
{
	int N = N_index;
	for (int i = 1; i < N - 1; i++)
	{
		u_jump[i] = u_plus[i] - u_minus[i - 1];
	}
	u_jump[0] = 0.0; u_jump[N - 1] = 0.0;
}

void compute_ave(double* u_ave, double* u_plus, double* u_minus, int N_index)
{
	int N = N_index;
	for (int i = 1; i < N - 1; i++)
	{
		u_ave[i] = (u_plus[i] + u_minus[i - 1]) / 2.0;
	}
	u_ave[0] = 0.0; u_ave[N - 1] = 0.0;
}

void compute_uxi(double* u, double* u_1, double* u_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		u[i] = u_1[i] + 3.0 * u_2[i] * xi;
	}
}

double Compute_alpha_max(double* alpha_c)
{
	double temp_max = 0.0;
	for (int i = 0; i < N_x + 1; i++)
	{
		if (alpha_c[i] > temp_max)
		{
			temp_max = alpha_c[i];
		}
	}
	return temp_max;
}

double Compute_alpha_min(double* alpha_c)
{
	double temp_min = 10.0;
	for (int i = 0; i < N_x + 1; i++)
	{
		if (alpha_c[i] < temp_min)
		{
			temp_min = alpha_c[i];
		}
	}
	return temp_min;
}

// Convection part
void Fd_hat(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* alpha_c)
{
	double* Fd_plus = new double[N_x + 2]; double* Fd_minus = new double[N_x + 2];
	double* d_plus = new double[N_x + 2]; double* d_minus = new double[N_x + 2];
	F_d(Fd_plus, m_0, m_1, m_2, -1.0, N_x + 2); F_d(Fd_minus, m_0, m_1, m_2, 1.0, N_x + 2);
	F_d(d_plus, rho_0, rho_1, rho_2, -1.0, N_x + 2); F_d(d_minus, rho_0, rho_1, rho_2, 1.0, N_x + 2);
	for (int i = 0; i < N_x + 1; i++)
	{
		u[i] = 0.5 * (Fd_plus[i + 1] + Fd_minus[i] - 1.0 * alpha_c[i] * (d_plus[i + 1] - d_minus[i]));
	}
	delete[] Fd_plus; delete[] Fd_minus; delete[] d_plus; delete[] d_minus;
}

void convection(double* U_0, double* U_1, double* U_2, double* F_hat, double* GP_0, double* GP_1, double* GP_2)
{
	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = -F_hat[i] + F_hat[i - 1] + GP_0[i];
		U_1[i] = -F_hat[i] - F_hat[i - 1] + GP_1[i];
		U_2[i] = -F_hat[i] + F_hat[i - 1] + GP_2[i];
	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;
}

void convection_d(double* con_d_0, double* con_d_1, double* con_d_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* alpha_c, double h_x)
{
	double* fd_hat = new double[N_x + 1];
	Fd_hat(fd_hat, rho_0, rho_1, rho_2, m_0, m_1, m_2, alpha_c);
	double* Fd_1 = new double[N_x + 2]; double* Fd_2 = new double[N_x + 2]; double* Fd_3 = new double[N_x + 2];
	double* Fd_4 = new double[N_x + 2]; double* Fd_5 = new double[N_x + 2];
	F_d(Fd_1, m_0, m_1, m_2, GQxi_1, N_x + 2); F_d(Fd_2, m_0, m_1, m_2, GQxi_2, N_x + 2); F_d(Fd_3, m_0, m_1, m_2, GQxi_3, N_x + 2);
	F_d(Fd_4, m_0, m_1, m_2, GQxi_4, N_x + 2); F_d(Fd_5, m_0, m_1, m_2, GQxi_5, N_x + 2);
	double* Fd_GQ_0 = new double[N_x + 2]; double* Fd_GQ_1 = new double[N_x + 2]; double* Fd_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Fd_GQ_0, Fd_GQ_1, Fd_GQ_2, Fd_1, Fd_2, Fd_3, Fd_4, Fd_5, N_x + 2, h_x);
	double* Fd_GQP_0 = new double[N_x + 2]; double* Fd_GQP_1 = new double[N_x + 2]; double* Fd_GQP_2 = new double[N_x + 2];
	get_Gauss_quadrature_prime(Fd_GQP_0, Fd_GQP_1, Fd_GQP_2, Fd_GQ_0, Fd_GQ_1, N_x + 2, h_x);
	convection(con_d_0, con_d_1, con_d_2, fd_hat, Fd_GQP_0, Fd_GQP_1, Fd_GQP_2);

	delete[] fd_hat;
	delete[] Fd_1; delete[] Fd_2; delete[] Fd_3; delete[] Fd_4; delete[] Fd_5;
	delete[] Fd_GQ_0; delete[] Fd_GQ_1; delete[] Fd_GQ_2;
	delete[] Fd_GQP_0; delete[] Fd_GQP_1; delete[] Fd_GQP_2;
}

void F_m(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double g_t = gamma_e + gamma_i + gamma_r - 3.0;
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ee = (gamma_e - 1.0) * (Ee_0[i] + Ee_1[i] * xi + Ee_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ei = (gamma_i - 1.0) * (Ei_0[i] + Ei_1[i] * xi + Ei_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Er = (gamma_r - 1.0) * (Er_0[i] + Er_1[i] * xi + Er_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = temp_rhou2 + temp_Ee + temp_Ei + temp_Er - g_t * temp_rhou2 / 6.0;
	}
}

void Fm_hat(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double* alpha_c)
{
	double* Fm_plus = new double[N_x + 2]; double* Fm_minus = new double[N_x + 2];
	double* m_plus = new double[N_x + 2]; double* m_minus = new double[N_x + 2];
	F_m(Fm_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, -1.0);
	F_m(Fm_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, 1.0);
	F_d(m_plus, m_0, m_1, m_2, -1.0, N_x + 2); F_d(m_minus, m_0, m_1, m_2, 1.0, N_x + 2);
	for (int i = 0; i < N_x + 1; i++)
	{
		u[i] = 0.5 * (Fm_plus[i + 1] + Fm_minus[i] - 1.0 * alpha_c[i] * (m_plus[i + 1] - m_minus[i]));
	}
	delete[] Fm_plus; delete[] Fm_minus; delete[] m_plus; delete[] m_minus;
}

void convection_m(double* con_m_0, double* con_m_1, double* con_m_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double* alpha_c, double h_x)
{
	double* fm_hat = new double[N_x + 1];
	Fm_hat(fm_hat, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, alpha_c);
	double* Fm_1 = new double[N_x + 2]; double* Fm_2 = new double[N_x + 2]; double* Fm_3 = new double[N_x + 2];
	double* Fm_4 = new double[N_x + 2]; double* Fm_5 = new double[N_x + 2];
	F_m(Fm_1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, GQxi_1);
	F_m(Fm_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, GQxi_2);
	F_m(Fm_3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, GQxi_3);
	F_m(Fm_4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, GQxi_4);
	F_m(Fm_5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, GQxi_5);
	double* Fm_GQ_0 = new double[N_x + 2]; double* Fm_GQ_1 = new double[N_x + 2]; double* Fm_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Fm_GQ_0, Fm_GQ_1, Fm_GQ_2, Fm_1, Fm_2, Fm_3, Fm_4, Fm_5, N_x + 2, h_x);
	double* Fm_GQP_0 = new double[N_x + 2]; double* Fm_GQP_1 = new double[N_x + 2]; double* Fm_GQP_2 = new double[N_x + 2];
	get_Gauss_quadrature_prime(Fm_GQP_0, Fm_GQP_1, Fm_GQP_2, Fm_GQ_0, Fm_GQ_1, N_x + 2, h_x);
	convection(con_m_0, con_m_1, con_m_2, fm_hat, Fm_GQP_0, Fm_GQP_1, Fm_GQP_2);

	delete[] fm_hat;
	delete[] Fm_1; delete[] Fm_2; delete[] Fm_3; delete[] Fm_4; delete[] Fm_5;
	delete[] Fm_GQ_0; delete[] Fm_GQ_1; delete[] Fm_GQ_2;
	delete[] Fm_GQP_0; delete[] Fm_GQP_1; delete[] Fm_GQP_2;
}

void F_alpha(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_u = (m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = gamma_alpha * (Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = (temp_Ealpha - (gamma_alpha - 1.0) * temp_rhou2 / 6.0) * temp_u;
	}
}

void Falpha_hat(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double* alpha_c)
{
	double* Falpha_plus = new double[N_x + 2]; double* Falpha_minus = new double[N_x + 2];
	double* Ealpha_plus = new double[N_x + 2]; double* Ealpha_minus = new double[N_x + 2];
	F_alpha(Falpha_plus, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, -1.0);
	F_alpha(Falpha_minus, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, 1.0);
	F_d(Ealpha_plus, Ealpha_0, Ealpha_1, Ealpha_2, -1.0, N_x + 2); F_d(Ealpha_minus, Ealpha_0, Ealpha_1, Ealpha_2, 1.0, N_x + 2);
	for (int i = 0; i < N_x + 1; i++)
	{
		u[i] = 0.5 * (Falpha_plus[i + 1] + Falpha_minus[i] - 1.0 * alpha_c[i] * (Ealpha_plus[i + 1] - Ealpha_minus[i]));
	}
	delete[] Falpha_plus; delete[] Falpha_minus; delete[] Ealpha_plus; delete[] Ealpha_minus;
}

void convection_alpha(double* con_alpha_0, double* con_alpha_1, double* con_alpha_2, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double* alpha_c, double h_x)
{
	double* falpha_hat = new double[N_x + 1];
	Falpha_hat(falpha_hat, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, alpha_c);
	double* Falpha_1 = new double[N_x + 2]; double* Falpha_2 = new double[N_x + 2]; double* Falpha_3 = new double[N_x + 2];
	double* Falpha_4 = new double[N_x + 2]; double* Falpha_5 = new double[N_x + 2];
	F_alpha(Falpha_1, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, GQxi_1);
	F_alpha(Falpha_2, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, GQxi_2);
	F_alpha(Falpha_3, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, GQxi_3);
	F_alpha(Falpha_4, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, GQxi_4);
	F_alpha(Falpha_5, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, GQxi_5);
	double* Falpha_GQ_0 = new double[N_x + 2]; double* Falpha_GQ_1 = new double[N_x + 2]; double* Falpha_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Falpha_GQ_0, Falpha_GQ_1, Falpha_GQ_2, Falpha_1, Falpha_2, Falpha_3, Falpha_4, Falpha_5, N_x + 2, h_x);
	double* Falpha_GQP_0 = new double[N_x + 2]; double* Falpha_GQP_1 = new double[N_x + 2]; double* Falpha_GQP_2 = new double[N_x + 2];
	get_Gauss_quadrature_prime(Falpha_GQP_0, Falpha_GQP_1, Falpha_GQP_2, Falpha_GQ_0, Falpha_GQ_1, N_x + 2, h_x);
	convection(con_alpha_0, con_alpha_1, con_alpha_2, falpha_hat, Falpha_GQP_0, Falpha_GQP_1, Falpha_GQP_2);

	delete[] falpha_hat;
	delete[] Falpha_1; delete[] Falpha_2; delete[] Falpha_3; delete[] Falpha_4; delete[] Falpha_5;
	delete[] Falpha_GQ_0; delete[] Falpha_GQ_1; delete[] Falpha_GQ_2;
	delete[] Falpha_GQP_0; delete[] Falpha_GQP_1; delete[] Falpha_GQP_2;
}

// Nonconservetive part
void p_alpha(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = (gamma_alpha - 1.0) * (temp_Ealpha - temp_rhou2 / 6.0);
	}
}

void p_alphaxi(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_1, double* Ealpha_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_u = (m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealphaxi = Ealpha_1[i] + 3.0 * Ealpha_2[i] * xi;
		double temp_rhoxi = rho_1[i] + 3.0 * rho_2[i] * xi;
		double temp_mxi = m_1[i] + 3.0 * m_2[i] * xi;
		u[i] = (gamma_alpha - 1.0) * (temp_Ealphaxi + temp_rhoxi * pow(temp_u, 2) / 6.0 - temp_mxi * temp_u / 3.0);
	}
}

void F_Ne(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi)
{
	double* pe_xi = new double[N_x + 2]; double* pi_xi = new double[N_x + 2]; double* pr_xi = new double[N_x + 2];
	p_alphaxi(pe_xi, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, xi);
	p_alphaxi(pi_xi, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_1, Ei_2, xi);
	p_alphaxi(pr_xi, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_1, Er_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_u = (m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = temp_u / 3.0 * (2 * pe_xi[i] - pi_xi[i] - pr_xi[i]);
	}

	delete[] pe_xi; delete[] pi_xi; delete[] pr_xi;
}

void F_Ni(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi)
{
	double* pe_xi = new double[N_x + 2]; double* pi_xi = new double[N_x + 2]; double* pr_xi = new double[N_x + 2];
	p_alphaxi(pe_xi, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, xi);
	p_alphaxi(pi_xi, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_1, Ei_2, xi);
	p_alphaxi(pr_xi, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_1, Er_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_u = (m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = temp_u / 3.0 * (2 * pi_xi[i] - pe_xi[i] - pr_xi[i]);
	}

	delete[] pe_xi; delete[] pi_xi; delete[] pr_xi;
}

void F_Nr(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_1, double* Ee_2, double* Ei_1, double* Ei_2, double* Er_1, double* Er_2, double xi)
{
	double* pe_xi = new double[N_x + 2]; double* pi_xi = new double[N_x + 2]; double* pr_xi = new double[N_x + 2];
	p_alphaxi(pe_xi, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, xi);
	p_alphaxi(pi_xi, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_1, Ei_2, xi);
	p_alphaxi(pr_xi, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_1, Er_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_u = (m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		u[i] = temp_u / 3.0 * (2 * pr_xi[i] - pe_xi[i] - pi_xi[i]);
	}

	delete[] pe_xi; delete[] pi_xi; delete[] pr_xi;
}

void palpha_jump(double* u, double gamma_alpha, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2)
{
	double* palpha_plus = new double[N_x + 2]; double* palpha_minus = new double[N_x + 2];
	p_alpha(palpha_plus, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, -1.0);
	p_alpha(palpha_minus, gamma_alpha, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, 1.0);
	compute_jump(u, palpha_plus, palpha_minus, N_x + 3);

	delete[] palpha_plus; delete[] palpha_minus;
}

void Nonconservetive_Ne(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode)
{
	double* F_Ne1 = new double[N_x + 2]; double* F_Ne2 = new double[N_x + 2]; double* F_Ne3 = new double[N_x + 2];
	double* F_Ne4 = new double[N_x + 2]; double* F_Ne5 = new double[N_x + 2];
	F_Ne(F_Ne1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_1);
	F_Ne(F_Ne2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_2);
	F_Ne(F_Ne3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_3);
	F_Ne(F_Ne4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_4);
	F_Ne(F_Ne5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_5);

	double* Ne_GQ_0 = new double[N_x + 2]; double* Ne_GQ_1 = new double[N_x + 2]; double* Ne_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Ne_GQ_0, Ne_GQ_1, Ne_GQ_2, F_Ne1, F_Ne2, F_Ne3, F_Ne4, F_Ne5, N_x + 2, h_x);

	double* pe_jump = new double[N_x + 3]; double* pi_jump = new double[N_x + 3]; double* pr_jump = new double[N_x + 3];
	palpha_jump(pe_jump, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2);
	palpha_jump(pi_jump, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2);
	palpha_jump(pr_jump, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);

	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = 2.0 * Ne_GQ_0[i] / h_x - coef_L[i - 1] * (2.0 * pe_jump[i] - pi_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pe_jump[i + 1] - pi_jump[i + 1] - pr_jump[i + 1]);
		U_1[i] = 2.0 * Ne_GQ_1[i] / h_x + coef_L[i - 1] * (2.0 * pe_jump[i] - pi_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pe_jump[i + 1] - pi_jump[i + 1] - pr_jump[i + 1]);
		U_2[i] = 2.0 * Ne_GQ_2[i] / h_x - coef_L[i - 1] * (2.0 * pe_jump[i] - pi_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pe_jump[i + 1] - pi_jump[i + 1] - pr_jump[i + 1]);

	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;

	delete[] F_Ne1; delete[] F_Ne2; delete[] F_Ne3; delete[] F_Ne4; delete[] F_Ne5;
	delete[] Ne_GQ_0; delete[] Ne_GQ_1; delete[] Ne_GQ_2;
	delete[] pe_jump; delete[] pi_jump; delete[] pr_jump;

}

void Nonconservetive_Ni(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode)
{
	double* F_Ni1 = new double[N_x + 2]; double* F_Ni2 = new double[N_x + 2]; double* F_Ni3 = new double[N_x + 2];
	double* F_Ni4 = new double[N_x + 2]; double* F_Ni5 = new double[N_x + 2];
	F_Ni(F_Ni1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_1);
	F_Ni(F_Ni2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_2);
	F_Ni(F_Ni3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_3);
	F_Ni(F_Ni4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_4);
	F_Ni(F_Ni5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_5);

	double* Ni_GQ_0 = new double[N_x + 2]; double* Ni_GQ_1 = new double[N_x + 2]; double* Ni_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Ni_GQ_0, Ni_GQ_1, Ni_GQ_2, F_Ni1, F_Ni2, F_Ni3, F_Ni4, F_Ni5, N_x + 2, h_x);

	double* pe_jump = new double[N_x + 3]; double* pi_jump = new double[N_x + 3]; double* pr_jump = new double[N_x + 3];
	palpha_jump(pe_jump, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2);
	palpha_jump(pi_jump, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2);
	palpha_jump(pr_jump, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);

	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = 2.0 * Ni_GQ_0[i] / h_x - coef_L[i - 1] * (2.0 * pi_jump[i] - pe_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pi_jump[i + 1] - pe_jump[i + 1] - pr_jump[i + 1]);
		U_1[i] = 2.0 * Ni_GQ_1[i] / h_x + coef_L[i - 1] * (2.0 * pi_jump[i] - pe_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pi_jump[i + 1] - pe_jump[i + 1] - pr_jump[i + 1]);
		U_2[i] = 2.0 * Ni_GQ_2[i] / h_x - coef_L[i - 1] * (2.0 * pi_jump[i] - pe_jump[i] - pr_jump[i]) - coef_R[i - 1] * (2.0 * pi_jump[i + 1] - pe_jump[i + 1] - pr_jump[i + 1]);
	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;

	delete[] F_Ni1; delete[] F_Ni2; delete[] F_Ni3; delete[] F_Ni4; delete[] F_Ni5;
	delete[] Ni_GQ_0; delete[] Ni_GQ_1; delete[] Ni_GQ_2;
	delete[] pe_jump; delete[] pi_jump; delete[] pr_jump;
}

void Nonconservetive_Nr(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, int index_mode)
{
	double* F_Nr1 = new double[N_x + 2]; double* F_Nr2 = new double[N_x + 2]; double* F_Nr3 = new double[N_x + 2];
	double* F_Nr4 = new double[N_x + 2]; double* F_Nr5 = new double[N_x + 2];
	F_Nr(F_Nr1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_1);
	F_Nr(F_Nr2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_2);
	F_Nr(F_Nr3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_3);
	F_Nr(F_Nr4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_4);
	F_Nr(F_Nr5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_1, Ee_2, Ei_1, Ei_2, Er_1, Er_2, GQxi_5);

	double* Nr_GQ_0 = new double[N_x + 2]; double* Nr_GQ_1 = new double[N_x + 2]; double* Nr_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Nr_GQ_0, Nr_GQ_1, Nr_GQ_2, F_Nr1, F_Nr2, F_Nr3, F_Nr4, F_Nr5, N_x + 2, h_x);

	double* pe_jump = new double[N_x + 3]; double* pi_jump = new double[N_x + 3]; double* pr_jump = new double[N_x + 3];
	palpha_jump(pe_jump, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2);
	palpha_jump(pi_jump, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2);
	palpha_jump(pr_jump, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);

	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = 2.0 * Nr_GQ_0[i] / h_x - coef_L[i - 1] * (2.0 * pr_jump[i] - pe_jump[i] - pi_jump[i]) - coef_R[i - 1] * (2.0 * pr_jump[i + 1] - pe_jump[i + 1] - pi_jump[i + 1]);
		U_1[i] = 2.0 * Nr_GQ_1[i] / h_x + coef_L[i - 1] * (2.0 * pr_jump[i] - pe_jump[i] - pi_jump[i]) - coef_R[i - 1] * (2.0 * pr_jump[i + 1] - pe_jump[i + 1] - pi_jump[i + 1]);
		U_2[i] = 2.0 * Nr_GQ_2[i] / h_x - coef_L[i - 1] * (2.0 * pr_jump[i] - pe_jump[i] - pi_jump[i]) - coef_R[i - 1] * (2.0 * pr_jump[i + 1] - pe_jump[i + 1] - pi_jump[i + 1]);
	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;

	delete[] F_Nr1; delete[] F_Nr2; delete[] F_Nr3; delete[] F_Nr4; delete[] F_Nr5;
	delete[] Nr_GQ_0; delete[] Nr_GQ_1; delete[] Nr_GQ_2;
	delete[] pe_jump; delete[] pi_jump; delete[] pr_jump;
}

// Source part
void compute_Tr4(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Er = Er_0[i] + Er_1[i] * xi + Er_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = (temp_Er - temp_rhou2 / 6.0) / a;
	}
}

void compute_Talpha(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_rho = rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = (temp_Ealpha - temp_rhou2 / 6.0) / (C_valpha * temp_rho);
	}
}

void F_Si(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double xi)
{
	double* Te = new double[N_x + 2]; double* Ti = new double[N_x + 2];
	compute_Talpha(Te, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, C_ve, xi);
	compute_Talpha(Ti, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2, C_vi, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		u[i] = omega_ei * (Te[i] - Ti[i]);
	}

	delete[] Te; delete[] Ti;
}

void F_Sr(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	double* Te = new double[N_x + 2]; double* Tr4 = new double[N_x + 2];
	compute_Talpha(Te, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, C_ve, xi);
	compute_Tr4(Tr4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		u[i] = omega_er * (pow(Te[i], 4) - Tr4[i]);
	}
	delete[] Te; delete[] Tr4;
}

void Source_i(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double h_x)
{
	double* F_Si1 = new double[N_x + 2]; double* F_Si2 = new double[N_x + 2]; double* F_Si3 = new double[N_x + 2];
	double* F_Si4 = new double[N_x + 2]; double* F_Si5 = new double[N_x + 2];
	F_Si(F_Si1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, GQxi_1);
	F_Si(F_Si2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, GQxi_2);
	F_Si(F_Si3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, GQxi_3);
	F_Si(F_Si4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, GQxi_4);
	F_Si(F_Si5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, GQxi_5);

	double* Si_GQ_0 = new double[N_x + 2]; double* Si_GQ_1 = new double[N_x + 2]; double* Si_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Si_GQ_0, Si_GQ_1, Si_GQ_2, F_Si1, F_Si2, F_Si3, F_Si4, F_Si5, N_x + 2, h_x);

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = Si_GQ_0[i];
		U_1[i] = Si_GQ_1[i];
		U_2[i] = Si_GQ_2[i];
	}
	delete[] F_Si1; delete[] F_Si2; delete[] F_Si3; delete[] F_Si4; delete[] F_Si5;
	delete[] Si_GQ_0; delete[] Si_GQ_1; delete[] Si_GQ_2;
}

void Source_r(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Er_0, double* Er_1, double* Er_2, double h_x)
{
	double* F_Sr1 = new double[N_x + 2]; double* F_Sr2 = new double[N_x + 2]; double* F_Sr3 = new double[N_x + 2];
	double* F_Sr4 = new double[N_x + 2]; double* F_Sr5 = new double[N_x + 2];
	F_Sr(F_Sr1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, GQxi_1);
	F_Sr(F_Sr2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, GQxi_2);
	F_Sr(F_Sr3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, GQxi_3);
	F_Sr(F_Sr4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, GQxi_4);
	F_Sr(F_Sr5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, GQxi_5);

	double* Sr_GQ_0 = new double[N_x + 2]; double* Sr_GQ_1 = new double[N_x + 2]; double* Sr_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Sr_GQ_0, Sr_GQ_1, Sr_GQ_2, F_Sr1, F_Sr2, F_Sr3, F_Sr4, F_Sr5, N_x + 2, h_x);

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = Sr_GQ_0[i];
		U_1[i] = Sr_GQ_1[i];
		U_2[i] = Sr_GQ_2[i];
	}

	delete[] F_Sr1; delete[] F_Sr2; delete[] F_Sr3; delete[] F_Sr4; delete[] F_Sr5;
	delete[] Sr_GQ_0; delete[] Sr_GQ_1; delete[] Sr_GQ_2;
}

void Source_e(double* U_0, double* U_1, double* U_2, double* Si_0, double* Si_1, double* Si_2, double* Sr_0, double* Sr_1, double* Sr_2)
{
	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = -1.0 * (Si_0[i] + Sr_0[i]);
		U_1[i] = -1.0 * (Si_1[i] + Sr_1[i]);
		U_2[i] = -1.0 * (Si_2[i] + Sr_2[i]);
	}
}


// Diffusion part
void compute_Talphaxi(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi)
{
	double* rhoxi = new double[N_x + 2]; double* Ealphaxi = new double[N_x + 2]; double* mxi = new double[N_x + 2];
	compute_uxi(rhoxi, rho_1, rho_2, xi); compute_uxi(Ealphaxi, Ealpha_1, Ealpha_2, xi);
	compute_uxi(mxi, m_1, m_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_rho = rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_m = m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = -1.0 * rhoxi[i] * (temp_Ealpha - temp_rhou2 / 6.0) / (C_valpha * pow(temp_rho, 2)) + (Ealphaxi[i] + (pow(temp_m, 2) * rhoxi[i]) / (6.0 * pow(temp_rho, 2)) - (temp_m * mxi[i]) / (3.0 * temp_rho)) / (C_valpha * temp_rho);
	}

	delete[] rhoxi; delete[] Ealphaxi; delete[] mxi;
}

void compute_Tr4xi(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	double* rhoxi = new double[N_x + 2]; double* Erxi = new double[N_x + 2]; double* mxi = new double[N_x + 2];
	compute_uxi(rhoxi, rho_1, rho_2, xi); compute_uxi(Erxi, Er_1, Er_2, xi); compute_uxi(mxi, m_1, m_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rho = rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_m = m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = (Erxi[i] + (pow(temp_m, 2) * rhoxi[i]) / (6.0 * pow(temp_rho, 2)) - (temp_m * mxi[i]) / (3.0 * temp_rho)) / a;
	}

	delete[] rhoxi; delete[] Erxi; delete[] mxi;
}

void compute_Talphaxi2(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double xi)
{
	double* rhoxi = new double[N_x + 2]; double* Ealphaxi = new double[N_x + 2]; double* mxi = new double[N_x + 2];
	compute_uxi(rhoxi, rho_1, rho_2, xi); compute_uxi(Ealphaxi, Ealpha_1, Ealpha_2, xi);
	compute_uxi(mxi, m_1, m_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhoxi2 = 3.0 * rho_2[i]; double temp_Ealphaxi2 = 3.0 * Ealpha_2[i]; double temp_mxi2 = 3.0 * m_2[i];
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_rho = rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_m = m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp1 = (2.0 * pow(rhoxi[i], 2) / pow(temp_rho, 3) - temp_rhoxi2 / pow(temp_rho, 2)) * (temp_Ealpha - temp_rhou2 / 6.0) / C_valpha;
		double temp2 = 2.0 * rhoxi[i] * (Ealphaxi[i] + (pow(temp_m, 2) * rhoxi[i]) / (6.0 * pow(temp_rho, 2)) - (temp_m * mxi[i]) / (3.0 * temp_rho)) / (C_valpha * pow(temp_rho, 2));
		double temp3 = (temp_Ealphaxi2 - (pow(temp_m, 2) * pow(rhoxi[i], 2)) / (3.0 * pow(temp_rho, 3)) + (2.0 * temp_m * mxi[i] * rhoxi[i]) / (3.0 * pow(temp_rho, 2)) + (pow(temp_m, 2) * temp_rhoxi2) / (6.0 * pow(temp_rho, 2)) - (pow(mxi[i], 2) + temp_m * temp_mxi2) / (3.0 * temp_rho)) / (C_valpha * temp_rho);
		u[i] = temp1 - temp2 + temp3;
	}

	delete[] rhoxi; delete[] Ealphaxi; delete[] mxi;
}

void compute_Tr4xi2(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double xi)
{
	double* rhoxi = new double[N_x + 2]; double* Erxi = new double[N_x + 2]; double* mxi = new double[N_x + 2];
	compute_uxi(rhoxi, rho_1, rho_2, xi); compute_uxi(Erxi, Er_1, Er_2, xi); compute_uxi(mxi, m_1, m_2, xi);
	for (int i = 0; i < N_x + 2; i++)
	{
		double temp_rhoxi2 = 3.0 * rho_2[i]; double temp_Erxi2 = 3.0 * Er_2[i]; double temp_mxi2 = 3.0 * m_2[i];
		double temp_rho = rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		double temp_m = m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		u[i] = (temp_Erxi2 - (pow(temp_m, 2) * pow(rhoxi[i], 2)) / (3.0 * pow(temp_rho, 3)) + (2.0 * temp_m * mxi[i] * rhoxi[i]) / (3.0 * pow(temp_rho, 2)) + (pow(temp_m, 2) * temp_rhoxi2) / (6.0 * pow(temp_rho, 2)) - (pow(mxi[i], 2) + temp_m * temp_mxi2) / (3.0 * temp_rho)) / a;
	}

	delete[] rhoxi; delete[] Erxi; delete[] mxi;
}

void Talpha_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha)
{
	double* Talpha_plus = new double[N_x + 2]; double* Talpha_minus = new double[N_x + 2];
	compute_Talpha(Talpha_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, -1.0);
	compute_Talpha(Talpha_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, 1.0);
	compute_jump(u, Talpha_plus, Talpha_minus, N_x + 3);

	delete[] Talpha_plus; delete[] Talpha_minus;
}

void Talphaxi_average(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha)
{
	double* Talphaxi_plus = new double[N_x + 2]; double* Talphaxi_minus = new double[N_x + 2];
	compute_Talphaxi(Talphaxi_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, -1.0);
	compute_Talphaxi(Talphaxi_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, 1.0);
	compute_ave(u, Talphaxi_plus, Talphaxi_minus, N_x + 3);

	delete[] Talphaxi_plus; delete[] Talphaxi_minus;
}

void Talphaxi2_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha)
{
	double* Talphaxi2_plus = new double[N_x + 2]; double* Talphaxi2_minus = new double[N_x + 2];
	compute_Talphaxi2(Talphaxi2_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, -1.0);
	compute_Talphaxi2(Talphaxi2_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, 1.0);
	compute_jump(u, Talphaxi2_plus, Talphaxi2_minus, N_x + 3);

	delete[] Talphaxi2_plus; delete[] Talphaxi2_minus;
}

void Tr4_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* Tr4_plus = new double[N_x + 2]; double* Tr4_minus = new double[N_x + 2];
	compute_Tr4(Tr4_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, -1.0);
	compute_Tr4(Tr4_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, 1.0);
	compute_jump(u, Tr4_plus, Tr4_minus, N_x + 3);

	delete[] Tr4_plus; delete[] Tr4_minus;
}

void Tr4xi_average(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* Tr4xi_plus = new double[N_x + 2]; double* Tr4xi_minus = new double[N_x + 2];
	compute_Tr4xi(Tr4xi_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, -1.0);
	compute_Tr4xi(Tr4xi_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, 1.0);
	compute_ave(u, Tr4xi_plus, Tr4xi_minus, N_x + 3);

	delete[] Tr4xi_plus; delete[] Tr4xi_minus;
}

void Tr4xi2_jump(double* u, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* Tr4xi2_plus = new double[N_x + 2]; double* Tr4xi2_minus = new double[N_x + 2];
	compute_Tr4xi2(Tr4xi2_plus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, -1.0);
	compute_Tr4xi2(Tr4xi2_minus, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, 1.0);
	compute_jump(u, Tr4xi2_plus, Tr4xi2_minus, N_x + 3);

	delete[] Tr4xi2_plus; delete[] Tr4xi2_minus;
}

void Diffusion_alpha(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2, double C_valpha, double kappa_alpha, double diff_alpha, double h_x)
{
	double* Talpha_j = new double[N_x + 3]; double* Talphaxi_a = new double[N_x + 3]; double* Talphaxi2_j = new double[N_x + 3];
	Talpha_jump(Talpha_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha);
	Talphaxi_average(Talphaxi_a, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha);
	Talphaxi2_jump(Talphaxi2_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha);

	double* Talphaxi_1 = new double[N_x + 2]; double* Talphaxi_2 = new double[N_x + 2]; double* Talphaxi_3 = new double[N_x + 2];
	double* Talphaxi_4 = new double[N_x + 2]; double* Talphaxi_5 = new double[N_x + 2];
	compute_Talphaxi(Talphaxi_1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, GQxi_1);
	compute_Talphaxi(Talphaxi_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, GQxi_2);
	compute_Talphaxi(Talphaxi_3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, GQxi_3);
	compute_Talphaxi(Talphaxi_4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, GQxi_4);
	compute_Talphaxi(Talphaxi_5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2, C_valpha, GQxi_5);

	double* Talphaxi_GQ_0 = new double[N_x + 2]; double* Talphaxi_GQ_1 = new double[N_x + 2]; double* Talphaxi_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Talphaxi_GQ_0, Talphaxi_GQ_1, Talphaxi_GQ_2, Talphaxi_1, Talphaxi_2, Talphaxi_3, Talphaxi_4, Talphaxi_5, N_x + 2, h_x);
	double* Talphaxi_GQP_0 = new double[N_x + 2]; double* Talphaxi_GQP_1 = new double[N_x + 2]; double* Talphaxi_GQP_2 = new double[N_x + 2];
	get_Gauss_quadrature_prime(Talphaxi_GQP_0, Talphaxi_GQP_1, Talphaxi_GQP_2, Talphaxi_GQ_0, Talphaxi_GQ_1, N_x + 2, h_x);

	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = -2.0 * kappa_alpha / h_x * Talphaxi_GQP_0[i] + 2.0 * kappa_alpha / h_x * (Talphaxi_a[i + 1] - Talphaxi_a[i]) + diff_alpha / h_x * (Talpha_j[i + 1] - Talpha_j[i]) + 4.0 * diff_beta * kappa_alpha / h_x * (Talphaxi2_j[i + 1] - Talphaxi2_j[i]);
		U_1[i] = -2.0 * kappa_alpha / h_x * Talphaxi_GQP_1[i] + 2.0 * kappa_alpha / h_x * (Talphaxi_a[i + 1] + Talphaxi_a[i]) + (diff_alpha - kappa_alpha) / h_x * (Talpha_j[i + 1] + Talpha_j[i]) + 4.0 * diff_beta * kappa_alpha / h_x * (Talphaxi2_j[i + 1] + Talphaxi2_j[i]);
		U_2[i] = -2.0 * kappa_alpha / h_x * Talphaxi_GQP_2[i] + 2.0 * kappa_alpha / h_x * (Talphaxi_a[i + 1] - Talphaxi_a[i]) + (diff_alpha - 3.0 * kappa_alpha) / h_x * (Talpha_j[i + 1] - Talpha_j[i]) + 4.0 * diff_beta * kappa_alpha / h_x * (Talphaxi2_j[i + 1] - Talphaxi2_j[i]);
	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;

	delete[] Talpha_j; delete[] Talphaxi_a; delete[] Talphaxi2_j;
	delete[] Talphaxi_1; delete[] Talphaxi_2; delete[] Talphaxi_3; delete[] Talphaxi_4; delete[] Talphaxi_5;
	delete[] Talphaxi_GQ_0; delete[] Talphaxi_GQ_1; delete[] Talphaxi_GQ_2;
	delete[] Talphaxi_GQP_0; delete[] Talphaxi_GQP_1; delete[] Talphaxi_GQP_2;
}

void Diffusion_r(double* U_0, double* U_1, double* U_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Er_0, double* Er_1, double* Er_2, double diff_alpha, double h_x)
{
	double* Tr4_j = new double[N_x + 3]; double* Tr4xi_a = new double[N_x + 3]; double* Tr4xi2_j = new double[N_x + 3];
	Tr4_jump(Tr4_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);
	Tr4xi_average(Tr4xi_a, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);
	Tr4xi2_jump(Tr4xi2_j, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);

	double* Tr4xi_1 = new double[N_x + 2]; double* Tr4xi_2 = new double[N_x + 2]; double* Tr4xi_3 = new double[N_x + 2];
	double* Tr4xi_4 = new double[N_x + 2]; double* Tr4xi_5 = new double[N_x + 2];
	compute_Tr4xi(Tr4xi_1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, GQxi_1);
	compute_Tr4xi(Tr4xi_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, GQxi_2);
	compute_Tr4xi(Tr4xi_3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, GQxi_3);
	compute_Tr4xi(Tr4xi_4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, GQxi_4);
	compute_Tr4xi(Tr4xi_5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, GQxi_5);

	double* Tr4xi_GQ_0 = new double[N_x + 2]; double* Tr4xi_GQ_1 = new double[N_x + 2]; double* Tr4xi_GQ_2 = new double[N_x + 2];
	get_Gauss_quadrature(Tr4xi_GQ_0, Tr4xi_GQ_1, Tr4xi_GQ_2, Tr4xi_1, Tr4xi_2, Tr4xi_3, Tr4xi_4, Tr4xi_5, N_x + 2, h_x);
	double* Tr4xi_GQP_0 = new double[N_x + 2]; double* Tr4xi_GQP_1 = new double[N_x + 2]; double* Tr4xi_GQP_2 = new double[N_x + 2];
	get_Gauss_quadrature_prime(Tr4xi_GQP_0, Tr4xi_GQP_1, Tr4xi_GQP_2, Tr4xi_GQ_0, Tr4xi_GQ_1, N_x + 2, h_x);

	for (int i = 1; i < N_x + 1; i++)
	{
		U_0[i] = -2.0 * kappa_r / h_x * Tr4xi_GQP_0[i] + 2.0 * kappa_r / h_x * (Tr4xi_a[i + 1] - Tr4xi_a[i]) + diff_alpha / h_x * (Tr4_j[i + 1] - Tr4_j[i]) + 4.0 * diff_beta * kappa_r / h_x * (Tr4xi2_j[i + 1] - Tr4xi2_j[i]);
		U_1[i] = -2.0 * kappa_r / h_x * Tr4xi_GQP_1[i] + 2.0 * kappa_r / h_x * (Tr4xi_a[i + 1] + Tr4xi_a[i]) + (diff_alpha - kappa_r) / h_x * (Tr4_j[i + 1] + Tr4_j[i]) + 4.0 * diff_beta * kappa_r / h_x * (Tr4xi2_j[i + 1] + Tr4xi2_j[i]);
		U_2[i] = -2.0 * kappa_r / h_x * Tr4xi_GQP_2[i] + 2.0 * kappa_r / h_x * (Tr4xi_a[i + 1] - Tr4xi_a[i]) + (diff_alpha - 3.0 * kappa_r) / h_x * (Tr4_j[i + 1] - Tr4_j[i]) + 4.0 * diff_beta * kappa_r / h_x * (Tr4xi2_j[i + 1] - Tr4xi2_j[i]);
	}
	U_0[0] = 0.0; U_1[0] = 0.0; U_2[0] = 0.0;
	U_0[N_x + 1] = 0.0; U_1[N_x + 1] = 0.0; U_2[N_x + 1] = 0.0;

	delete[] Tr4_j; delete[] Tr4xi_a; delete[] Tr4xi2_j;
	delete[] Tr4xi_1; delete[] Tr4xi_2; delete[] Tr4xi_3; delete[] Tr4xi_4; delete[] Tr4xi_5;
	delete[] Tr4xi_GQ_0; delete[] Tr4xi_GQ_1; delete[] Tr4xi_GQ_2;
	delete[] Tr4xi_GQP_0; delete[] Tr4xi_GQP_1; delete[] Tr4xi_GQP_2;
}

void space_d(double* U_0, double* U_1, double* U_2, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double h_x)
{
	double* con_d_0 = new double[N_x + 2]; double* con_d_1 = new double[N_x + 2]; double* con_d_2 = new double[N_x + 2];

	convection_d(con_d_0, con_d_1, con_d_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, alpha_c, h_x);

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = con_d_0[i];
		U_1[i] = con_d_1[i];
		U_2[i] = con_d_2[i];
	}

	delete[] con_d_0; delete[] con_d_1; delete[] con_d_2;
}

void space_m(double* U_0, double* U_1, double* U_2, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x)
{
	//Convection part
	double* con_m_0 = new double[N_x + 2]; double* con_m_1 = new double[N_x + 2]; double* con_m_2 = new double[N_x + 2];

	convection_m(con_m_0, con_m_1, con_m_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, alpha_c, h_x);

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = con_m_0[i];
		U_1[i] = con_m_1[i];
		U_2[i] = con_m_2[i];
	}

	delete[] con_m_0; delete[] con_m_1; delete[] con_m_2;
}

void space_e(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode)
{
	double* con_e_0 = new double[N_x + 2]; double* con_e_1 = new double[N_x + 2]; double* con_e_2 = new double[N_x + 2];
	double* Non_e_0 = new double[N_x + 2]; double* Non_e_1 = new double[N_x + 2]; double* Non_e_2 = new double[N_x + 2];
	double* S_e_0 = new double[N_x + 2]; double* S_e_1 = new double[N_x + 2]; double* S_e_2 = new double[N_x + 2];
	double* D_e_0 = new double[N_x + 2]; double* D_e_1 = new double[N_x + 2]; double* D_e_2 = new double[N_x + 2];
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

	convection_alpha(con_e_0, con_e_1, con_e_2, gamma_e, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, alpha_c, h_x);

	Nonconservetive_Ne(Non_e_0, Non_e_1, Non_e_2, coef_L, coef_R, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, index_mode);

	double* S_i_0 = new double[N_x + 2]; double* S_i_1 = new double[N_x + 2]; double* S_i_2 = new double[N_x + 2];
	double* S_r_0 = new double[N_x + 2]; double* S_r_1 = new double[N_x + 2]; double* S_r_2 = new double[N_x + 2];
	Source_i(S_i_0, S_i_1, S_i_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, h_x);
	Source_r(S_r_0, S_r_1, S_r_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, h_x);
	Source_e(S_e_0, S_e_1, S_e_2, S_i_0, S_i_1, S_i_2, S_r_0, S_r_1, S_r_2);
	delete[] S_i_0; delete[] S_i_1; delete[] S_i_2; delete[] S_r_0; delete[] S_r_1; delete[] S_r_2;

	Diffusion_alpha(D_e_0, D_e_1, D_e_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, C_ve, kappa_e, diff_alpha_e, h_x);


	// Accuracy source part
	double* S3_GQ_0 = new double[N_x + 2]; double* S3_GQ_1 = new double[N_x + 2]; double* S3_GQ_2 = new double[N_x + 2];
	double* S3_x1 = new double[N_x + 2]; double* S3_x2 = new double[N_x + 2]; double* S3_x3 = new double[N_x + 2];
	double* S3_x4 = new double[N_x + 2]; double* S3_x5 = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		S3_x1[i] = S3_ex(X_1[i], t); S3_x2[i] = S3_ex(X_2[i], t); S3_x3[i] = S3_ex(X_3[i], t);
		S3_x4[i] = S3_ex(X_4[i], t); S3_x5[i] = S3_ex(X_5[i], t);
	}
	get_Gauss_quadrature(S3_GQ_0, S3_GQ_1, S3_GQ_2, S3_x1, S3_x2, S3_x3, S3_x4, S3_x5, N_x + 2, h_x);
	delete[] S3_x1; delete[] S3_x2; delete[] S3_x3; delete[] S3_x4; delete[] S3_x5;

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = con_e_0[i] + Non_e_0[i] + S_e_0[i] + D_e_0[i] + S3_GQ_0[i];
		U_1[i] = con_e_1[i] + Non_e_1[i] + S_e_1[i] + D_e_1[i] + S3_GQ_1[i];
		U_2[i] = con_e_2[i] + Non_e_2[i] + S_e_2[i] + D_e_2[i] + S3_GQ_2[i];
	}


	delete[] con_e_0; delete[] con_e_1; delete[] con_e_2;
	delete[] Non_e_0; delete[] Non_e_1; delete[] Non_e_2;
	delete[] S_e_0; delete[] S_e_1; delete[] S_e_2;
	delete[] D_e_0; delete[] D_e_1; delete[] D_e_2;
	delete[] S3_GQ_0; delete[] S3_GQ_1; delete[] S3_GQ_2;
}

void space_i(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode)
{
	double* con_i_0 = new double[N_x + 2]; double* con_i_1 = new double[N_x + 2]; double* con_i_2 = new double[N_x + 2];
	double* Non_i_0 = new double[N_x + 2]; double* Non_i_1 = new double[N_x + 2]; double* Non_i_2 = new double[N_x + 2];
	double* S_i_0 = new double[N_x + 2]; double* S_i_1 = new double[N_x + 2]; double* S_i_2 = new double[N_x + 2];
	double* D_i_0 = new double[N_x + 2]; double* D_i_1 = new double[N_x + 2]; double* D_i_2 = new double[N_x + 2];
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

	convection_alpha(con_i_0, con_i_1, con_i_2, gamma_i, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2, alpha_c, h_x);

	Nonconservetive_Ni(Non_i_0, Non_i_1, Non_i_2, coef_L, coef_R, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, index_mode);

	Source_i(S_i_0, S_i_1, S_i_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, h_x);

	Diffusion_alpha(D_i_0, D_i_1, D_i_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2, C_vi, kappa_i, diff_alpha_i, h_x);

	// Accuracy source part
	double* S4_GQ_0 = new double[N_x + 2]; double* S4_GQ_1 = new double[N_x + 2]; double* S4_GQ_2 = new double[N_x + 2];
	double* S4_x1 = new double[N_x + 2]; double* S4_x2 = new double[N_x + 2]; double* S4_x3 = new double[N_x + 2];
	double* S4_x4 = new double[N_x + 2]; double* S4_x5 = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		S4_x1[i] = S4_ex(X_1[i], t); S4_x2[i] = S4_ex(X_2[i], t); S4_x3[i] = S4_ex(X_3[i], t);
		S4_x4[i] = S4_ex(X_4[i], t); S4_x5[i] = S4_ex(X_5[i], t);
	}
	get_Gauss_quadrature(S4_GQ_0, S4_GQ_1, S4_GQ_2, S4_x1, S4_x2, S4_x3, S4_x4, S4_x5, N_x + 2, h_x);
	delete[] S4_x1; delete[] S4_x2; delete[] S4_x3; delete[] S4_x4; delete[] S4_x5;


	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = con_i_0[i] + Non_i_0[i] + S_i_0[i] + D_i_0[i] + S4_GQ_0[i];
		U_1[i] = con_i_1[i] + Non_i_1[i] + S_i_1[i] + D_i_1[i] + S4_GQ_1[i];
		U_2[i] = con_i_2[i] + Non_i_2[i] + S_i_2[i] + D_i_2[i] + S4_GQ_2[i];
	}

	delete[] con_i_0; delete[] con_i_1; delete[] con_i_2;
	delete[] Non_i_0; delete[] Non_i_1; delete[] Non_i_2;
	delete[] S_i_0; delete[] S_i_1; delete[] S_i_2;
	delete[] D_i_0; delete[] D_i_1; delete[] D_i_2;
	delete[] S4_GQ_0; delete[] S4_GQ_1; delete[] S4_GQ_2;
}

void space_r(double* U_0, double* U_1, double* U_2, double* coef_L, double* coef_R, double* X_1, double* X_2, double* X_3, double* X_4, double* X_5, double* alpha_c, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2, double h_x, double t, int index_mode)
{
	double* con_r_0 = new double[N_x + 2]; double* con_r_1 = new double[N_x + 2]; double* con_r_2 = new double[N_x + 2];
	double* Non_r_0 = new double[N_x + 2]; double* Non_r_1 = new double[N_x + 2]; double* Non_r_2 = new double[N_x + 2];
	double* S_r_0 = new double[N_x + 2]; double* S_r_1 = new double[N_x + 2]; double* S_r_2 = new double[N_x + 2];
	double* D_r_0 = new double[N_x + 2]; double* D_r_1 = new double[N_x + 2]; double* D_r_2 = new double[N_x + 2];
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

	convection_alpha(con_r_0, con_r_1, con_r_2, gamma_r, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, alpha_c, h_x);

	Nonconservetive_Nr(Non_r_0, Non_r_1, Non_r_2, coef_L, coef_R, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2, h_x, index_mode);

	Source_r(S_r_0, S_r_1, S_r_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Er_0, Er_1, Er_2, h_x);

	Diffusion_r(D_r_0, D_r_1, D_r_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2, diff_alpha_r, h_x);

	// Accuracy source part
	double* S5_GQ_0 = new double[N_x + 2]; double* S5_GQ_1 = new double[N_x + 2]; double* S5_GQ_2 = new double[N_x + 2];
	double* S5_x1 = new double[N_x + 2]; double* S5_x2 = new double[N_x + 2]; double* S5_x3 = new double[N_x + 2];
	double* S5_x4 = new double[N_x + 2]; double* S5_x5 = new double[N_x + 2];
	for (int i = 0; i < N_x + 2; i++)
	{
		S5_x1[i] = S5_ex(X_1[i], t); S5_x2[i] = S5_ex(X_2[i], t); S5_x3[i] = S5_ex(X_3[i], t);
		S5_x4[i] = S5_ex(X_4[i], t); S5_x5[i] = S5_ex(X_5[i], t);
	}
	get_Gauss_quadrature(S5_GQ_0, S5_GQ_1, S5_GQ_2, S5_x1, S5_x2, S5_x3, S5_x4, S5_x5, N_x + 2, h_x);
	delete[] S5_x1; delete[] S5_x2; delete[] S5_x3; delete[] S5_x4; delete[] S5_x5;

	for (int i = 0; i < N_x + 2; i++)
	{
		U_0[i] = con_r_0[i] + Non_r_0[i] + S_r_0[i] + D_r_0[i] + S5_GQ_0[i];
		U_1[i] = con_r_1[i] + Non_r_1[i] + S_r_1[i] + D_r_1[i] + S5_GQ_1[i];
		U_2[i] = con_r_2[i] + Non_r_2[i] + S_r_2[i] + D_r_2[i] + S5_GQ_2[i];
	}

	delete[] con_r_0; delete[] con_r_1; delete[] con_r_2;
	delete[] Non_r_0; delete[] Non_r_1; delete[] Non_r_2;
	delete[] S_r_0; delete[] S_r_1; delete[] S_r_2;
	delete[] D_r_0; delete[] D_r_1; delete[] D_r_2;
	delete[] S5_GQ_0; delete[] S5_GQ_1; delete[] S5_GQ_2;
}

