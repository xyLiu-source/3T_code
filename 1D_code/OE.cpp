#define _CRT_SECURE_NO_WARNINGS 1
#include "modalDG.h"
#include "space.h"
#include "OE.h"

void OE_compute_u(double* u, double xi, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2)
{
	double* m = new double[N_x + 2]; double* rho = new double[N_x + 2];
	F_d(m, m_0, m_1, m_2, xi, N_x + 2); F_d(rho, rho_0, rho_1, rho_2, xi, N_x + 2);

	for (int i = 1; i < N_x + 1; i++)
	{
		u[i - 1] = m[i] / rho[i];
	}
	delete[] m; delete[] rho;
}

void OE_compute_u_ave(double* u_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2)
{
	double* u_GQ1 = new double[N_x]; double* u_GQ2 = new double[N_x]; double* u_GQ3 = new double[N_x];
	double* u_GQ4 = new double[N_x]; double* u_GQ5 = new double[N_x];

	OE_compute_u(u_GQ1, GQxi_1, rho_0, rho_1, rho_2, m_0, m_1, m_2);
	OE_compute_u(u_GQ2, GQxi_2, rho_0, rho_1, rho_2, m_0, m_1, m_2);
	OE_compute_u(u_GQ3, GQxi_3, rho_0, rho_1, rho_2, m_0, m_1, m_2);
	OE_compute_u(u_GQ4, GQxi_4, rho_0, rho_1, rho_2, m_0, m_1, m_2);
	OE_compute_u(u_GQ5, GQxi_5, rho_0, rho_1, rho_2, m_0, m_1, m_2);

	for (int i = 0; i < N_x; i++)
	{
		u_ave[i] = 0.5 * (weight_1 * u_GQ1[i] + weight_2 * u_GQ2[i] + weight_3 * u_GQ3[i] + weight_4 * u_GQ4[i] + weight_5 * u_GQ5[i]);
	}

	delete[] u_GQ1; delete[] u_GQ2; delete[] u_GQ3; delete[] u_GQ4; delete[] u_GQ5;
}

void OE_compute_e_alpha(double* e_alpha, double xi, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2)
{
	for (int i = 1; i < N_x + 1; i++)
	{
		double temp_rhou2 = pow((m_0[i] + m_1[i] * xi + m_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0), 2) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
		double temp_Ealpha = Ealpha_0[i] + Ealpha_1[i] * xi + Ealpha_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0;
		e_alpha[i - 1] = (temp_Ealpha - temp_rhou2 / 6.0) / (rho_0[i] + rho_1[i] * xi + rho_2[i] * (3.0 * pow(xi, 2) - 1.0) / 2.0);
	}
}

void OE_compute_e_alpha_ave(double* e_alpha_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ealpha_0, double* Ealpha_1, double* Ealpha_2)
{
	double* e_alpha_GQ1 = new double[N_x]; double* e_alpha_GQ2 = new double[N_x]; double* e_alpha_GQ3 = new double[N_x];
	double* e_alpha_GQ4 = new double[N_x]; double* e_alpha_GQ5 = new double[N_x];

	OE_compute_e_alpha(e_alpha_GQ1, GQxi_1, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2);
	OE_compute_e_alpha(e_alpha_GQ2, GQxi_2, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2);
	OE_compute_e_alpha(e_alpha_GQ3, GQxi_3, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2);
	OE_compute_e_alpha(e_alpha_GQ4, GQxi_4, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2);
	OE_compute_e_alpha(e_alpha_GQ5, GQxi_5, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ealpha_0, Ealpha_1, Ealpha_2);

	for (int i = 0; i < N_x; i++)
	{
		e_alpha_ave[i] = 0.5 * (weight_1 * e_alpha_GQ1[i] + weight_2 * e_alpha_GQ2[i] + weight_3 * e_alpha_GQ3[i] + weight_4 * e_alpha_GQ4[i] + weight_5 * e_alpha_GQ5[i]);
	}

	delete[] e_alpha_GQ1; delete[] e_alpha_GQ2; delete[] e_alpha_GQ3; delete[] e_alpha_GQ4; delete[] e_alpha_GQ5;
}

void OE_compute_Cs_ave(double* Cs_ave, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* e_e_ave = new double[N_x]; double* e_i_ave = new double[N_x]; double* e_r_ave = new double[N_x];

	OE_compute_e_alpha_ave(e_e_ave, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2);
	OE_compute_e_alpha_ave(e_i_ave, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ei_0, Ei_1, Ei_2);
	OE_compute_e_alpha_ave(e_r_ave, rho_0, rho_1, rho_2, m_0, m_1, m_2, Er_0, Er_1, Er_2);

	for (int i = 0; i < N_x; i++)
	{
		Cs_ave[i] = sqrt(gamma_e * (gamma_e - 1.0) * e_e_ave[i] + gamma_i * (gamma_i - 1.0) * e_i_ave[i] + gamma_r * (gamma_r - 1.0) * e_r_ave[i]);
	}
	delete[] e_e_ave; delete[] e_i_ave; delete[] e_r_ave;
}

void OE_compute_beta_j(double* beta_j, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* u_ave = new double[N_x]; double* Cs_ave = new double[N_x];

	OE_compute_u_ave(u_ave, rho_0, rho_1, rho_2, m_0, m_1, m_2);
	OE_compute_Cs_ave(Cs_ave, rho_0, rho_1, rho_2, m_0, m_1, m_2, Ee_0, Ee_1, Ee_2, Ei_0, Ei_1, Ei_2, Er_0, Er_1, Er_2);

	for (int i = 0; i < N_x; i++)
	{
		double temp_max = abs(u_ave[i] - Cs_ave[i]);
		if (abs(u_ave[i]) > temp_max)
		{
			temp_max = abs(u_ave[i]);
		}
		if (abs(u_ave[i] + Cs_ave[i]) > temp_max)
		{
			temp_max = abs(u_ave[i] + Cs_ave[i]);
		}
		beta_j[i] = temp_max;
	}

	delete[] u_ave; delete[] Cs_ave;
}

void OE_compute_scalar_sigma(double* Linf, double* sigma_0, double* sigma_1, double* sigma_2, double* u_0, double* u_1, double* u_2)
{
	double u_Gave = 0.0;
	for (int i = 1; i < N_x + 1; i++)
	{
		u_Gave = u_Gave + u_0[i];
	}
	u_Gave = u_Gave / N_x;
	double* U_L = new double[N_x + 2]; double* U_R = new double[N_x + 2];
	double* U_GQ1 = new double[N_x + 2]; double* U_GQ2 = new double[N_x + 2]; double* U_GQ3 = new double[N_x + 2];
	double* U_GQ4 = new double[N_x + 2]; double* U_GQ5 = new double[N_x + 2];

	double* u_L = new double[N_x]; double* u_R = new double[N_x];
	double* u_GQ1 = new double[N_x]; double* u_GQ2 = new double[N_x]; double* u_GQ3 = new double[N_x];
	double* u_GQ4 = new double[N_x]; double* u_GQ5 = new double[N_x];

	F_d(U_L, u_0, u_1, u_2, -1.0, N_x + 2); F_d(U_R, u_0, u_1, u_2, 1.0, N_x + 2);
	F_d(U_GQ1, u_0, u_1, u_2, GQxi_1, N_x + 2); F_d(U_GQ2, u_0, u_1, u_2, GQxi_2, N_x + 2); F_d(U_GQ3, u_0, u_1, u_2, GQxi_3, N_x + 2);
	F_d(U_GQ4, u_0, u_1, u_2, GQxi_4, N_x + 2); F_d(U_GQ5, u_0, u_1, u_2, GQxi_5, N_x + 2);

	for (int i = 1; i < N_x + 1; i++)
	{
		u_L[i - 1] = U_L[i]; u_R[i - 1] = U_R[i];
		u_GQ1[i - 1] = U_GQ1[i]; u_GQ2[i - 1] = U_GQ2[i]; u_GQ3[i - 1] = U_GQ3[i];
		u_GQ4[i - 1] = U_GQ4[i]; u_GQ5[i - 1] = U_GQ5[i];
	}

	double u_GLinf = 0.0;
	for (int i = 0; i < N_x; i++)
	{
		if (abs(u_L[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_L[i] - u_Gave);
		}
		if (abs(u_GQ1[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_GQ1[i] - u_Gave);
		}
		if (abs(u_GQ2[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_GQ2[i] - u_Gave);
		}
		if (abs(u_GQ3[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_GQ3[i] - u_Gave);
		}
		if (abs(u_GQ4[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_GQ4[i] - u_Gave);
		}
		if (abs(u_GQ5[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_GQ5[i] - u_Gave);
		}
		if (abs(u_R[i] - u_Gave) > u_GLinf)
		{
			u_GLinf = abs(u_R[i] - u_Gave);
		}
	}
	Linf[0] = u_GLinf;
	//cout << u_Gave << "  " << u_GLinf << endl;

	if (u_GLinf < pow(10, -12))
	{
		for (int i = 0; i < N_x; i++)
		{
			sigma_0[i] = 0.0; sigma_1[i] = 0.0; sigma_2[i] = 0.0;
		}
	}
	else
	{
		double* U_plus = new double[N_x + 2]; double* U_minus = new double[N_x + 2]; double* U_jump = new double[N_x + 3];
		F_d(U_plus, u_0, u_1, u_2, -1.0, N_x + 2); F_d(U_minus, u_0, u_1, u_2, 1.0, N_x + 2); compute_jump(U_jump, U_plus, U_minus, N_x + 3);

		double* uxi_plus = new double[N_x + 2]; double* uxi_minus = new double[N_x + 2]; double* uxi_jump = new double[N_x + 3];

		compute_uxi(uxi_plus, u_1, u_2, -1.0); compute_uxi(uxi_minus, u_1, u_2, 1.0); compute_jump(uxi_jump, uxi_plus, uxi_minus, N_x + 3);

		double* uxi2_plus = new double[N_x + 2]; double* uxi2_minus = new double[N_x + 2]; double* uxi2_jump = new double[N_x + 3];
		for (int i = 0; i < N_x + 2; i++)
		{
			uxi2_plus[i] = 3.0 * u_2[i];
			uxi2_minus[i] = uxi2_plus[i];
		}
		compute_jump(uxi2_jump, uxi2_plus, uxi2_minus, N_x + 3);

		u_GLinf = 1.0;

		for (int i = 0; i < N_x; i++)
		{
			sigma_0[i] = Omega * 1.0 / (6.0 * u_GLinf) * (abs(U_jump[i + 1]) + abs(U_jump[i + 2]));
			sigma_1[i] = Omega * 1.0 / u_GLinf * (abs(uxi_jump[i + 1]) + abs(uxi_jump[i + 2]));
			sigma_2[i] = Omega * 5.0 / (3.0 * u_GLinf) * (abs(uxi2_jump[i + 1]) + abs(uxi2_jump[i + 2]));
		}

		delete[] U_plus; delete[] U_minus; delete[] U_jump;
		delete[] uxi_plus; delete[] uxi_minus; delete[] uxi_jump;
		delete[] uxi2_plus; delete[] uxi2_minus; delete[] uxi2_jump;
	}

	delete[] U_GQ1; delete[] U_GQ2; delete[] U_GQ3; delete[] U_GQ4; delete[] U_GQ5; delete[] U_L; delete[] U_R;
	delete[] u_GQ1; delete[] u_GQ2; delete[] u_GQ3; delete[] u_GQ4; delete[] u_GQ5; delete[] u_L; delete[] u_R;
}

void OE_compute_system_sigma(double* sigma_0, double* sigma_1, double* sigma_2, double* rho_0, double* rho_1, double* rho_2, double* m_0, double* m_1, double* m_2, double* Ee_0, double* Ee_1, double* Ee_2, double* Ei_0, double* Ei_1, double* Ei_2, double* Er_0, double* Er_1, double* Er_2)
{
	double* rho_sigma_0 = new double[N_x]; double* rho_sigma_1 = new double[N_x]; double* rho_sigma_2 = new double[N_x];
	double* m_sigma_0 = new double[N_x]; double* m_sigma_1 = new double[N_x]; double* m_sigma_2 = new double[N_x];
	double* Ee_sigma_0 = new double[N_x]; double* Ee_sigma_1 = new double[N_x]; double* Ee_sigma_2 = new double[N_x];
	double* Ei_sigma_0 = new double[N_x]; double* Ei_sigma_1 = new double[N_x]; double* Ei_sigma_2 = new double[N_x];
	double* Er_sigma_0 = new double[N_x]; double* Er_sigma_1 = new double[N_x]; double* Er_sigma_2 = new double[N_x];

	double* rhoL_inf = new double[1]; double* mL_inf = new double[1]; double* EeL_inf = new double[1];
	double* EiL_inf = new double[1]; double* ErL_inf = new double[1];
	OE_compute_scalar_sigma(rhoL_inf, rho_sigma_0, rho_sigma_1, rho_sigma_2, rho_0, rho_1, rho_2);
	OE_compute_scalar_sigma(mL_inf, m_sigma_0, m_sigma_1, m_sigma_2, m_0, m_1, m_2);
	OE_compute_scalar_sigma(EeL_inf, Ee_sigma_0, Ee_sigma_1, Ee_sigma_2, Ee_0, Ee_1, Ee_2);
	OE_compute_scalar_sigma(EiL_inf, Ei_sigma_0, Ei_sigma_1, Ei_sigma_2, Ei_0, Ei_1, Ei_2);
	OE_compute_scalar_sigma(ErL_inf, Er_sigma_0, Er_sigma_1, Er_sigma_2, Er_0, Er_1, Er_2);

	double L_inf = rhoL_inf[0];
	if (mL_inf[0] > L_inf)
	{
		L_inf = mL_inf[0];
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

	for (int i = 0; i < N_x; i++)
	{
		rho_sigma_0[i] = rho_sigma_0[i] / L_inf; rho_sigma_1[i] = rho_sigma_1[i] / L_inf; rho_sigma_2[i] = rho_sigma_2[i] / L_inf;
		m_sigma_0[i] = m_sigma_0[i] / L_inf; m_sigma_1[i] = m_sigma_1[i] / L_inf; m_sigma_2[i] = m_sigma_2[i] / L_inf;
		Ee_sigma_0[i] = Ee_sigma_0[i] / L_inf; Ee_sigma_1[i] = Ee_sigma_1[i] / L_inf; Ee_sigma_2[i] = Ee_sigma_2[i] / L_inf;
		Ei_sigma_0[i] = Ei_sigma_0[i] / L_inf; Ei_sigma_1[i] = Ei_sigma_1[i] / L_inf; Ei_sigma_2[i] = Ei_sigma_2[i] / L_inf;
		Er_sigma_0[i] = Er_sigma_0[i] / L_inf; Er_sigma_1[i] = Er_sigma_1[i] / L_inf; Er_sigma_2[i] = Er_sigma_2[i] / L_inf;

		double temp_max0 = 0.0;
		if (rho_sigma_0[i] > temp_max0)
		{
			temp_max0 = rho_sigma_0[i];
		}
		if (m_sigma_0[i] > temp_max0)
		{
			temp_max0 = m_sigma_0[i];
		}
		if (Ee_sigma_0[i] > temp_max0)
		{
			temp_max0 = Ee_sigma_0[i];
		}
		if (Ei_sigma_0[i] > temp_max0)
		{
			temp_max0 = Ei_sigma_0[i];
		}
		if (Er_sigma_0[i] > temp_max0)
		{
			temp_max0 = Er_sigma_0[i];
		}
		sigma_0[i] = temp_max0;

		double temp_max1 = 0.0;
		if (rho_sigma_1[i] > temp_max1)
		{
			temp_max1 = rho_sigma_1[i];
		}
		if (m_sigma_1[i] > temp_max1)
		{
			temp_max1 = m_sigma_1[i];
		}
		if (Ee_sigma_1[i] > temp_max1)
		{
			temp_max1 = Ee_sigma_1[i];
		}
		if (Ei_sigma_1[i] > temp_max1)
		{
			temp_max1 = Ei_sigma_1[i];
		}
		if (Er_sigma_1[i] > temp_max1)
		{
			temp_max1 = Er_sigma_1[i];
		}
		sigma_1[i] = temp_max1;

		double temp_max2 = 0.0;
		if (rho_sigma_2[i] > temp_max2)
		{
			temp_max2 = rho_sigma_2[i];
		}
		if (m_sigma_2[i] > temp_max2)
		{
			temp_max2 = m_sigma_2[i];
		}
		if (Ee_sigma_2[i] > temp_max2)
		{
			temp_max2 = Ee_sigma_2[i];
		}
		if (Ei_sigma_2[i] > temp_max2)
		{
			temp_max2 = Ei_sigma_2[i];
		}
		if (Er_sigma_2[i] > temp_max2)
		{
			temp_max2 = Er_sigma_2[i];
		}
		sigma_2[i] = temp_max2;
	}

	delete[] rho_sigma_0; delete[] rho_sigma_1; delete[] rho_sigma_2;
	delete[] m_sigma_0; delete[] m_sigma_1; delete[] m_sigma_2;
	delete[] Ee_sigma_0; delete[] Ee_sigma_1; delete[] Ee_sigma_2;
	delete[] Ei_sigma_0; delete[] Ei_sigma_1; delete[] Ei_sigma_2;
	delete[] Er_sigma_0; delete[] Er_sigma_1; delete[] Er_sigma_2;
	delete[] rhoL_inf; delete[] mL_inf; delete[] EeL_inf;
	delete[] EiL_inf; delete[] ErL_inf;
}

void OE_procedure(double* u_1, double* u_2, double* beta_j, double* sigma_0, double* sigma_1, double* sigma_2, double tau, double h_x)
{
	for (int i = 0; i < N_x; i++)
	{
		double C1 = exp(-1.0 * beta_j[i] * tau * (sigma_0[i] + sigma_1[i]) / h_x);
		double C2 = exp(-1.0 * beta_j[i] * tau * (sigma_0[i] + sigma_1[i] + sigma_2[i]) / h_x);
		u_1[i + 1] = C1 * u_1[i + 1];
		u_2[i + 1] = C2 * u_2[i + 1];
	}
}
