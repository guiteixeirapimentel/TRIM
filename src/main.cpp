#include <iostream>
#include <vector>
#include <fstream>
#include "ElevatorFunction.h"

#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"

#include "Math/Matriz.h"
#include "Math/MathUtils.h"

int main()
{
	const double S = 16.0;
	const double meanChord = 1.5;
	const double mass = 1250.0;
	const double momentInertia = 1825.0;
	const double rho = 1.0;
	const double gravity = 9.80665;

	DragCoefficient CD(0.02806, 0.121);
	LiftCoefficient CL(0.2977, 4.41, 1.7, 3.9, 0.43);
	MomentCoefficient CM(0.03917, -0.613, -7.27, -12.4, -1.122);

	XForce xf(CL, CD);
	ZForce zf(CL, CD);
	YMoment ym(CM);

	ThrustForce T(1100.0);

	f1 errorFunc1(xf, T, rho, gravity, S);
	f2 errorFunc2(zf, rho, gravity, S);
	f3 errorFunc3(ym, rho, gravity, S);
	f4 errorFunc4;

	State Xref(70.0, 0.0, 0.0, 0.0);
	StateDot Xdotref(0.0, 0.0, 0.0, 0.0);
	double deltaEref = 2 * M_PI / 180.0;

	std::vector<double> EE;
	EE.resize(4*4);

	EE[0] = errorFunc1.JacobianXdot1(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[1] = errorFunc1.JacobianXdot2(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[2] = errorFunc1.JacobianXdot3(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[3] = errorFunc1.JacobianXdot4(meanChord, Xref, Xdotref, mass, deltaEref);

	EE[0 + (4 * 1)] = errorFunc2.JacobianXdot1(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[1 + (4 * 1)] = errorFunc2.JacobianXdot2(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[2 + (4 * 1)] = errorFunc2.JacobianXdot3(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[3 + (4 * 1)] = errorFunc2.JacobianXdot4(meanChord, Xref, Xdotref, mass, deltaEref);

	EE[0 + (4 * 2)] = errorFunc3.JacobianXdot1(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	EE[1 + (4 * 2)] = errorFunc3.JacobianXdot2(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	EE[2 + (4 * 2)] = errorFunc3.JacobianXdot3(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	EE[3 + (4 * 2)] = errorFunc3.JacobianXdot4(meanChord, Xref, Xdotref, momentInertia, deltaEref);

	EE[0 + (4 * 3)] = errorFunc4.JacobianXdot1(Xref, Xdotref);
	EE[1 + (4 * 3)] = errorFunc4.JacobianXdot2(Xref, Xdotref);
	EE[2 + (4 * 3)] = errorFunc4.JacobianXdot3(Xref, Xdotref);
	EE[3 + (4 * 3)] = errorFunc4.JacobianXdot4(Xref, Xdotref);

	std::vector<double> AA;
	AA.resize(4*4);

	AA[0] = errorFunc1.JacobianX1(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[1] = errorFunc1.JacobianX2(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[2] = errorFunc1.JacobianX3(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[3] = errorFunc1.JacobianX4(meanChord, Xref, Xdotref, mass, deltaEref);

	AA[0 + (4 * 1)] = errorFunc2.JacobianX1(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[1 + (4 * 1)] = errorFunc2.JacobianX2(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[2 + (4 * 1)] = errorFunc2.JacobianX3(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[3 + (4 * 1)] = errorFunc2.JacobianX4(meanChord, Xref, Xdotref, mass, deltaEref);

	AA[0 + (4 * 2)] = errorFunc3.JacobianX1(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	AA[1 + (4 * 2)] = errorFunc3.JacobianX2(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	AA[2 + (4 * 2)] = errorFunc3.JacobianX3(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	AA[3 + (4 * 2)] = errorFunc3.JacobianX4(meanChord, Xref, Xdotref, momentInertia, deltaEref);

	AA[0 + (4 * 3)] = errorFunc4.JacobianX1(Xref, Xdotref);
	AA[1 + (4 * 3)] = errorFunc4.JacobianX2(Xref, Xdotref);
	AA[2 + (4 * 3)] = errorFunc4.JacobianX3(Xref, Xdotref);
	AA[3 + (4 * 3)] = errorFunc4.JacobianX4(Xref, Xdotref);

	std::vector<double> BB;
	BB.resize(4);

	BB[0] = errorFunc1.JacobianU1(meanChord, Xref, Xdotref, mass, deltaEref);
	
	BB[1] = errorFunc2.JacobianU1(meanChord, Xref, Xdotref, mass, deltaEref);
	
	BB[2] = errorFunc3.JacobianU1(meanChord, Xref, Xdotref, momentInertia, deltaEref);
	
	BB[3] = errorFunc4.JacobianU1(Xref, Xdotref);


	Matriz mAA(AA, 4, 4);
	Matriz mBB(BB, 4, 1);
	Matriz mEE(EE, 4, 4);


	Matriz mEEinv = CalcInvMatriz(mEE*-1);
	
	Matriz I = mEEinv*mEE;

	Matriz A = mEEinv * mAA;
	Matriz B = mEEinv * mBB;
	
	std::cout << "Matriz EE ";
	MostrarMatriz(mEE);
	std::cout << "\nMatriz AA ";
	MostrarMatriz(mAA);
	std::cout << "\nmatriz BB ";
	MostrarMatriz(mBB);

	std::cout << "\nMatriz A ";
	MostrarMatriz(A);
	std::cout << "\nMatriz B ";
	MostrarMatriz(B);


	std::vector<Matriz> estados;

	const double dt = 0.001;
	const double tempoSim = 60.0;
	const size_t numEstados = tempoSim/dt;

	ElevatorFunction posProf;

	estados.resize(numEstados);

	estados[0] = Matriz({0.0, 0.0, 0.0, 0.0}, 4, 1);
	
	for(size_t i = 1; i < numEstados; i++)
	{
		const double t = i * dt;
		Matriz estadosDot = (A*(estados[i-1])) + (B*(posProf(t)-deltaEref));

		estados[i] = estados[i-1] + (estadosDot * dt);
	} 


	std::ofstream arq("saida.csv");
	arq << "t;u;w;q [graus/s];theta [graus]\n";

	for(size_t i = 0; i < numEstados; i+=100)
	{
		arq << i*dt << ";";
		arq << (*estados[i].GetPtrMatriz())[0] + Xref.u;
		arq << ";" << (*estados[i].GetPtrMatriz())[1] + Xref.w;
		arq << ";" << ((*estados[i].GetPtrMatriz())[2] + Xref.q) * 180.0 / M_PI;
		arq << ";" << ((*estados[i].GetPtrMatriz())[3] + Xref.theta) * 180.0 / M_PI;
		arq <<"\n";
	}

	arq.close();

	Matriz mA = mEEinv * mAA;
	Matriz MB = mEEinv * mBB;


    return 0;
}