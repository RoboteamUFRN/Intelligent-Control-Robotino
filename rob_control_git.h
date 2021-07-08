Authorship: 	Gabriel S. Lima - gabriel.lima.095@ufrn.edu.br, 
		Victor R. F. Moreira - victorramon@ufrn.edu.br, 
		Wallace M. Bessa - wmobes@utu.fi

#define M_PI 3.14159265358979323846

arma::mat M(3, 3);

double l = 0.294;
double beta = M_PI / 3.0;
double r = 0.051;
double mr = 11.2; 
double Rr = 0.27;

//Derivative filter
double mi = 5.0; 
double gama = 1.0; 
double xi = 0.5; 
double estim_00 = 0.0;
double estim_01 = 0.0;
double estim_02 = 0.0;
double estim_10 = 0.0;
double estim_11 = 0.0;
double estim_12 = 0.0;

double csc(double x)
{
	double out;

	out = 1.0 / sin(x);

	return out;
}

double sec(double x)
{
	double out;

	out = 1.0 / cos(x);

	return out;
}

double sgn(double x)
{
	if (x > 0.0)
	{
		return 1.0;
	}
	if (x < 0.0)
	{
		return -1.0;
	}
	if (x == 0.0)
	{
		return 0.0;
	}
}

double sat(double x)
{
	if (x > 1.0)
	{
		return 1.0;
	}
	else if (x < -1.0)
	{
		return -1.0;
	}
	else
	{
		return x;
	}

}

double max(double a, double b)
{
	if (a > b)
	{
		return a;
	}
	else if (b > a)
	{
		return b;
	}
	else
	{
		return a;
	}
}

double min(double a, double b)
{
	if (a < b)
	{
		return a;
	}
	else if (b < a)
	{
		return b;
	}
	else
	{
		return a;
	}
}

/*Filter*/
void filtro(arma::vec tau, arma::vec omega, arma::vec& omega_p, double theta_tau, double theta_omega)
{
	omega_p = theta_tau * tau - theta_omega * omega;
}

void rk4_f(double t, double tf, arma::vec tau, arma::vec& omega, double theta_tau, double theta_omega)
{
	arma::vec k1_p(3), k2_p(3), k3_p(3), k4_p(3);
	arma::vec auxP(3), omega_p(3);
	double h;

	h = tf - t;

	filtro(tau, omega, omega_p, theta_tau, theta_omega);

	for (int i = 0; i < 3; i++)
	{
		k1_p(i) = h * omega_p(i);
	}

	auxP = 0.5 * k1_p + tau;
	filtro(auxP, omega, omega_p, theta_tau, theta_omega);

	for (int i = 0; i < 3; i++)
	{
		k2_p(i) = h * omega_p(i);
	}

	auxP = 0.5 * k2_p + tau;
	filtro(auxP, omega, omega_p, theta_tau, theta_omega);

	for (int i = 0; i < 3; i++)
	{
		k3_p(i) = h * omega_p(i);
	}

	auxP = k3_p + tau;
	filtro(auxP, omega, omega_p, theta_tau, theta_omega);

	for (int i = 0; i < 3; i++)
	{
		k4_p(i) = h * omega_p(i);
	}

	omega = omega + (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) / 6.0;
}

/*Derivative filter*/
double derivador(double x, double x_est)
{
	double xp_est;

	xp_est = -mi * sat((x_est - x) / xi);

	return xp_est;
}

double rk4_d(double t, double tf, double x, double* x_est)
{
	double k1_p, k2_p, k3_p, k4_p, X_est, xp_est;
	double h;
	h = tf - t;
	X_est = *x_est;

	xp_est = derivador(x, X_est);
	k1_p = h * xp_est;

	xp_est = derivador(x, X_est + k1_p / 2.0);
	k2_p = h * xp_est;

	xp_est = derivador(x, X_est + k2_p / 2.0);
	k3_p = h * xp_est;

	xp_est = derivador(x, X_est + k3_p);
	k4_p = h * xp_est;

	X_est = X_est + (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) / 6.0;

	xp_est = derivador(x, X_est);
	*x_est = X_est;

	return xp_est;
}

void matrizM()
{
	M(0, 0) = r * r * pow(sec(beta / 2.0), 4) * mr * ((2.0 * l * l * csc(beta / 2.0) * csc(beta / 2.0)) + (Rr * Rr)) / (32.0 * l * l);
	M(0, 1) = r * r * pow(sec(beta / 2.0), 4) * mr * ((cos(beta) * Rr * Rr) - (2.0 * l * l)) / (16.0 * l * l);
	M(0, 2) = r * r * pow(sec(beta / 2.0), 4) * mr * (((4.0 * cos(beta) * l * l) / (cos(beta) - 1.0)) + (Rr * Rr)) / (32.0 * l * l);
	M(1, 0) = M(0, 1);
	M(1, 1) = r * r * pow(sec(beta / 2.0), 4) * mr * ((cos(beta) * cos(beta) * Rr * Rr) + (2.0 * l * l)) / (8.0 * l * l);
	M(1, 2) = M(0, 1);
	M(2, 0) = M(0, 2);
	M(2, 1) = M(0, 1);
	M(2, 2) = M(0, 0);
}

void matrizTr(arma::mat& Tr, arma::vec Q)
{
	Tr(0, 0) = r * (sin(beta - Q(2)) + sin(Q(2))) / (cos(2.0 * beta) - 1.0);
	Tr(0, 1) = r * sin(Q(2)) / (cos(beta) + 1.0);
	Tr(0, 2) = r * (sin(Q(2)) - sin(beta + Q(2))) / (cos(2.0 * beta) - 1.0);
	Tr(1, 0) = r * (cos(beta - Q(2)) - cos(Q(2))) / (cos(2.0 * beta) - 1.0);
	Tr(1, 1) = -r * cos(Q(2)) / (cos(beta) + 1.0);
	Tr(1, 2) = r * (cos(beta + Q(2)) - cos(Q(2))) / (cos(2.0 * beta) - 1.0);
	Tr(2, 0) = r / (2.0 * l * (cos(beta) + 1.0));
	Tr(2, 1) = r * cos(beta) / (l * (cos(beta) + 1.0));
	Tr(2, 2) = Tr(2, 0);
}

void matrizTrp(arma::mat& Trp, arma::vec Q, arma::vec Qp)
{
	Trp(0, 0) = r * (-Qp(2) * cos(beta - Q(2)) + Qp(2) * cos(Q(2))) / (cos(2.0 * beta) - 1.0);
	Trp(0, 1) = r * Qp(2) * cos(Q(2)) / (cos(beta) + 1.0);
	Trp(0, 2) = r * (Qp(2) * cos(Q(2)) - Qp(2) * cos(beta + Q(2))) / (cos(2.0 * beta) - 1.0);
	Trp(1, 0) = r * (Qp(2) * sin(beta - Q(2)) + Qp(2) * sin(Q(2))) / (cos(2.0 * beta) - 1.0);
	Trp(1, 1) = r * Qp(2) * sin(Q(2)) / (cos(beta) + 1.0);
	Trp(1, 2) = r * (-Qp(2) * sin(beta + Q(2)) + Qp(2) * sin(Q(2))) / (cos(2.0 * beta) - 1.0);
	Trp(2, 0) = 0.0;
	Trp(2, 1) = 0.0;
	Trp(2, 2) = Trp(2, 0);
}

/*ANN*/
void ativ_gauss(arma::vec x, arma::mat center, arma::mat width, arma::vec cam, arma::mat& out)
{
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			out(i, j) = exp(-0.5 * (x(j) - center(i, j)) * (x(j) - center(i, j)) / (width(i, j) * width(i, j)));
		}
	}
}

void atual_gauss(double t, double tf, arma::vec input, arma::vec cam, arma::mat& weight, arma::mat center, arma::mat width, arma::vec learn, arma::vec& d_est, double limit_weight)
{
	arma::mat func_ativ(7, 3), wwT(7, 7);
	double h = tf - t;

	ativ_gauss(input, center, width, cam, func_ativ);

	for (int j = 0; j < 3; j++)
	{
		if (norm(weight.col(j)) < limit_weight || (norm(weight.col(j)) == limit_weight && learn(j) * input(j) * dot(weight.col(j), func_ativ.col(j)) <= 0.0))
		{
			weight.col(j) += learn(j) * input(j) * func_ativ.col(j) * h;
		}
		else
		{
			wwT = weight.col(j) * weight.col(j).t();
			weight.col(j) += (arma::eye(7, 7) - wwT / dot(weight.col(j), weight.col(j))) * learn(j) * input(j) * func_ativ.col(j) * h;
		}
		d_est(j) = dot(weight.col(j), func_ativ.col(j));
	}
}

/*Desired trajectories*/
void trajetoriaP(double t, double tf, double ampli, double freq, arma::vec& Pd)
{
	Pd(0) =   ampli * cos(freq*t);
	Pd(1) = - ampli * sin(freq * t);
	Pd(2) = - (freq * t) - M_PI/2.0;
}

void trajetoriaV(double t, double tf, double ampli, double freq, arma::vec& Pd, arma::vec& Vd)
{
	Vd(0) = - ampli * freq * sin(freq * t);
	Vd(1) = - ampli * freq * cos(freq * t);
	Vd(2) = - freq;
}

void trajetoriaA(double t, double tf, double ampli, double freq, arma::vec& Pd, arma::vec& Vd, arma::vec& Ad)
{
	Ad(0) = - ampli * pow(freq, 2) * cos(freq * t);
	Ad(1) =   ampli * pow(freq, 2) * sin(freq * t);
	Ad(2) = 0.0;
}

/*Controller*/
void FBL_ANN_gauss(double t, double tf, arma::vec Pd, arma::vec Vd, arma::vec Ad, arma::vec& Ep, arma::vec& Ev, arma::vec Q, arma::vec& Qp, arma::vec& Qfp, arma::vec& s, arma::vec cam, arma::mat center, arma::mat width, arma::mat& weight, arma::vec learn, double limit_weight, arma::vec& d_est, arma::vec& ang_Vel, arma::vec lambda, arma::vec& tau)
{
	const double theta_tau = 28.95;
	const double theta_omega = 1.0;
	const double alph1 = 1.0;
	const double alph2 = 1.0;

	arma::mat Tr(3, 3), Trp(3, 3), i_Tr(3, 3), matrixb(3, 3);
	arma::vec F(3);
	arma::vec tau_m(3), omega_m(3);

	matrizTr(Tr, Q);
	i_Tr = inv(Tr);
	matrixb = M * i_Tr;

	Qp(0) = rk4_d(t, tf, Q(0), &estim_00);
	Qp(1) = rk4_d(t, tf, Q(1), &estim_01);
	Qp(2) = rk4_d(t, tf, Q(2), &estim_02);

	rk4_f(t, tf, Qp, Qfp, alph1, alph2);

	Ep = Q - Pd;
	Ev = Qfp - Vd;

	for (int i = 0; i < 3; i++)
	{
		s(i) = Ev(i) + lambda(i) * Ep(i);
		ang_Vel(i) = Ad(i) - 2 * lambda(i) * Ev(i) - pow(lambda(i), 2) * Ep(i);
	}

	atual_gauss(t, tf, s, cam, weight, center, width, learn, d_est, limit_weight);

	matrizTrp(Trp, Q, Qfp);

	F = Trp * i_Tr * Qfp;

	ang_Vel = matrixb * (-F - d_est + ang_Vel);

	tau = ang_Vel;

	rk4_f(t, tf, ang_Vel, ang_Vel, theta_tau, theta_omega);

	ang_Vel = (ang_Vel * 60.0) / (2 * M_PI);

	for (int i = 0; i < 3; i++)
	{
		if (fabs(ang_Vel(i)) > 3600.0)
		{
			ang_Vel(i) = 3600.0 * sgn(ang_Vel(i));
		}
	}
}