/*
Authorship: 	Gabriel S. Lima - gabriel.lima.095@ufrn.edu.br, 
		Victor R. F. Moreira - victor.moreira.086@ufrn.edu.br, 
		Wallace M. Bessa - wmobes@utu.fi
*/

#define _USE_MATH_DEFINES // for C++
#define ARMA_USE_CXX11

#include <cstdlib>
#include <random>
#include <cmath>
#include <iostream>
#include <time.h>
#include <fstream>
#include <cstdio>
#include <chrono>
#include <armadillo>
#include <thread>
#include <iomanip>

#include "rec/robotino/com/all.h"
#include "rec/core_lt/utils.h"
#include "rec/core_lt/Timer.h"

FILE *out;

#include "rob_control.h"

class MyCom : public rec::robotino::com::Com
{
public:
	MyCom()
	{
	}

	void errorEvent(Error error, const char* errorString)
	{
		std::cerr << "Error: " << errorString << std::endl;
	}

	void connectedEvent()
	{
		std::cout << "Connected." << std::endl;
	}

	void connectionClosedEvent()
	{
		std::cout << "Connection closed." << std::endl;
	}

	void modeChangedEvent(bool isPassiveMode)
	{
		if (isPassiveMode)
		{
			std::cout << "Connected int passive mode." << std::endl;
		}
	}
};

MyCom com;
rec::core_lt::Timer timer;				
rec::robotino::com::Motor motor1;		
rec::robotino::com::Motor motor2;
rec::robotino::com::Motor motor3;
rec::robotino::com::Odometry odometry;
rec::robotino::com::Bumper bumper;

void init(const std::string& hostname)
{
	// Initialize the actors

	motor1.setComId(com.id());
	motor1.setMotorNumber(0);

	motor2.setComId(com.id());
	motor2.setMotorNumber(1);

	motor3.setComId(com.id());
	motor3.setMotorNumber(2);

	//PID - angular velocities control

	float kp = 0.95, ki = 0.1, kd = 0.01;

	motor1.setPID(kp, ki, kd);
	motor2.setPID(kp, ki, kd);
	motor3.setPID(kp, ki, kd);

	odometry.set(0.0, 0.0, 0.0);

	// Connect
	std::cout << "Connecting..." << std::endl;
	com.setAddress(hostname.c_str());

	com.connect();

	std::cout << std::endl << "Connected" << std::endl;

	//odometry.set(0.0f, 0.0f, 0.0f);
}

void drive()
{
	timer.start();			   

	double t_init = 1.0;
	double t = 0.001 * timer.msecsElapsed();		
	double tf = 59.99 + t_init;					
	double dt = 0.0;
	double t_old;
	int flag = 0;
				
	double atf = 0.0;
	double cont_ang = 0.0;
	double Q_old = 0.0;
	double deltaQ = 0.0;
	double Phi = 0.0;

	matrizM();

	arma::vec ang_Vel(3), tau(3);
	arma::vec Q(3), Qp(3), Qfp(3);
	arma::vec Ep(3), Ev(3);					

	double amp = 0.5;						
	double freq = 2.0 * M_PI / 30.0;		

	arma::vec Pd(3), Vd(3), Ad(3);			
	arma::vec lambda(3), cam(3), s(3);

	Qfp.fill(0.0);
	Qp.fill(0.0);

	Q(0) = odometry.x() * 0.001;	/* Q(0) -> x */
	Q(1) = odometry.y() * 0.001;	/* Q(1) -> y */
	Q(2) = 0.0;			/* Q(2) -> phi */
	Phi = (odometry.phi() * M_PI) / 180.0;

	lambda.fill(1.0);						
	cam = { 0.5, 0.5, 3.0 };

	/*ANN*/
	arma::mat weight(7, 3), center(7, 3), width(7, 3);
	arma::vec learn(3), d_est(3);
	weight.fill(0.0);
	learn = { 30.0, 30.0, 35.0 };

	d_est.fill(0.0);

	//Centers
	center.row(0) = -cam.t() / 2.0;
	center.row(1) = -cam.t() / 4.0;
	center.row(2) = -cam.t() / 8.0;
	center.row(3).fill(0.0);
	center.row(4) = cam.t() / 8.0;
	center.row(5) = cam.t() / 4.0;
	center.row(6) = cam.t() / 2.0;

	//Widths
	width.row(0) = cam.t() / 3.0;
	width.row(1) = cam.t() / 5.0;
	width.row(2) = cam.t() / 15.0;
	width.row(3) = cam.t() / 30.0;
	width.row(4) = cam.t() / 15.0;
	width.row(5) = cam.t() / 5.0;
	width.row(6) = cam.t() / 3.0;

	double limit_weight = 0.1;					

	fprintf(out, "[1]t\t    ");
	fprintf(out, "[2]Q(0)\t[3]Q(1)\t[4]Q(2)\t");						
	fprintf(out, "[5]Qp(0)\t[6]Qp(2)\t[7]Qp(2)\t");					
	fprintf(out, "[8]Pd(0)\t[9]Pd(1)\t[10]Pd(2)\t");						
	fprintf(out, "[11]Vd(0)\t[12]Vd(1)\t[13]Vd(2)\t");						
	fprintf(out, "[14]Ad(0)\t[15]Ad(1)\t[16]Ad(2)\t");						
	fprintf(out, "[17]Ep(0)\t[18]Ep(1)\t[19]Ep(2)\t");						
	fprintf(out, "[20]Ev(0)\t[21]Ev(1)\t[22]Ev(2)\t");						
	fprintf(out, "[23]s(0)\t[24]s(1)\t[25]s(2)\t");						
	fprintf(out, "[26]tau0\t[27]tau1\t[28]tau2\t");		
	fprintf(out, "[29]angV0\t[30]angV1\t[31]angV2\t");		
	fprintf(out, "[32]d(0)\t[33]d(1)\t[34]d(2)\t");			
	fprintf(out, "\n");


	/*Loop*/
	while (t < tf && false == bumper.value())
	{
		t_old = t;
		t = 0.001 * timer.msecsElapsed();
		dt += t - t_old;

		if (dt >= 0.2)	
		{
			t_ant = t_c;
			t_c = t;
			dt_c = (t_c - t_ant);

			Q(0) = odometry.x() * 0.001;
			Q(1) = odometry.y() * 0.001;

			Q_old = Phi;
			Phi = ((odometry.phi() * M_PI) / 180.0);

			if (flag == 0)
			{
				Q(2) = Phi;
				flag = 1;
			}

			deltaQ = Phi - Q_old;

			if (fabs(deltaQ) >= 2.0)
			{
				Q(2) += deltaQ - 2.0 * M_PI * sgn(deltaQ);
			}
			else
			{
				Q(2) += deltaQ;
			}


			//Circunference
			trajetoriaP(t - t_init, tf - t_init, amp, freq, Pd);
			trajetoriaV(t - t_init, tf - t_init, amp, freq, Pd, Vd);
			trajetoriaA(t - t_init, tf - t_init, amp, freq, Pd, Vd, Ad);

			/*Controller*/
			if (t < t_init)
			{
				ang_Vel.fill(0.0);
				tau.fill(0.0);
			}
			else
			{
				FBL_ANN_gauss(t - t_init, t - t_init + dt, Pd, Vd, Ad, Ep, Ev, Q, Qp, Qfp, s, cam, center, width, weight, learn, limit_weight, d_est, ang_Vel, lambda, tau);
			}

			/*Control signal*/
			motor1.setSpeedSetPoint(ang_Vel(0));		/* ang_Vel(0) -> Motor 1 [RPM]*/
			motor2.setSpeedSetPoint(ang_Vel(1));		/* ang_Vel(0) -> Motor 2 [RPM]*/
			motor3.setSpeedSetPoint(ang_Vel(2));		/* ang_Vel(0) -> Motor 3 [RPM]*/

			if (t > atf && t > t_init)
			{
				fprintf(out, "%.6f ", t - t_init);										
				fprintf(out, "%.6f %.6f %.6f ", Q(0), Q(1), Q(2));						
				fprintf(out, "%.6f %.6f %.6f ", Qfp(0), Qfp(1), Qfp(2));					
				fprintf(out, "%.6f %.6f %.6f ", Pd(0), Pd(1), Pd(2));						
				fprintf(out, "%.6f %.6f %.6f ", Vd(0), Vd(1), Vd(2));						
				fprintf(out, "%.6f %.6f %.6f ", Ad(0), Ad(1), Ad(2));						
				fprintf(out, "%.6f %.6f %.6f ", Ep(0), Ep(1), Ep(2));						
				fprintf(out, "%.6f %.6f %.6f ", Ev(0), Ev(1), Ev(2));						
				fprintf(out, "%.6f %.6f %.6f ", s(0), s(1), s(2));						
				fprintf(out, "%.6f %.6f %.6f ", tau(0), tau(1), tau(2));		
				fprintf(out, "%.6f %.6f %.6f ", ang_Vel(0), ang_Vel(1), ang_Vel(2));		
				fprintf(out, "%.6f %.6f %.6f ", d_est(0), d_est(1), d_est(2));			
				fprintf(out, "\n");
				atf = atf + ast;
			}
			dt = 0.0;
		}
	}
}

void destroy()
{
	com.disconnect();
}

void print()
{

}

int main(int argc, char** argv)
{
	std::string hostname = "172.26.201.2";	// Robotino Real

	errno_t arquivo;
		
	arquivo = fopen_s( &out,"results.txt", "w");

	if (!out)
	{
		perror("File opening failed");
		return EXIT_FAILURE;
	}

	if (argc > 1)
	{
		hostname = argv[1];
	}

	try
	{
		init(hostname);
		drive();
		print();
		destroy();
	}

	catch (const rec::robotino::com::ComException& e)
	{
		std::cerr << "Com Error: " << e.what() << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknow Error" << std::endl;
	}

	fclose(out);

	std::cout << "Press any key to exit..." << std::endl;

	std::cin.get();
}
