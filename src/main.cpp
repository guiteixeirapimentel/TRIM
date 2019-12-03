#include <iostream>
#include <vector>
#include <fstream>
#include "ElevatorFunction.h"

#include "ftrim1.h"
#include "ftrim2.h"
#include "ftrim3.h"
#include "ftrim4.h"


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

	DragCoefficient CD(0.02, 0.07);
	LiftCoefficient CL(0.2977, 4.41, 1.7, 3.9, 0.43);
	MomentCoefficient CM(0.03917, -0.613, -7.27, -12.4, -1.122);

	XForce xf(CL, CD);
	ZForce zf(CL, CD);
	YMoment ym(CM);

	ThrustForce T(135e3);

	ftrim1 errorFunc1(xf, T, rho, gravity, S, meanChord, mass);
	ftrim2 errorFunc2(zf, rho, gravity, S, meanChord, mass);
	ftrim3 errorFunc3(ym, rho, gravity, S, meanChord, momentInertia);
	ftrim4 errorFunc4;

	//State Xref(70.0, 0.0, 0.0, 0.0);
	//StateDot Xdotref(0.0, 0.0, 0.0, 0.0);
	//double deltaEref = 0.0;

	// 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
	std::vector<double> args = {33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 33.0, 0.0, 0.0};
	
	// 0->alfa; 1->deltaE; 2->deltaT; 3->q
	Matriz X({0.0, 0.0, 1.0, 0.0}, 4, 1);
	Matriz Xp(X);

	double norma = 1;
	size_t n = 0;
	while(norma > 1e-7)
	{
		Matriz J = MatrizI(4);

		args[11] = (*X.GetPtrMatriz())[0];
		args[8] = (*X.GetPtrMatriz())[1];
		args[9] = (*X.GetPtrMatriz())[2];
		args[2] = (*X.GetPtrMatriz())[3];

		args[3] = args[11];
		args[1] = tan(args[11])*args[10];
		args[0] = sqrt((args[10]*args[10])-((args[1]*args[1])));

		(*J.GetPtrMatriz())[0 + (4 * 0)] = errorFunc1.GetPartialDerivative(args, 11);
		(*J.GetPtrMatriz())[0 + (4 * 1)] = errorFunc2.GetPartialDerivative(args, 11);
		(*J.GetPtrMatriz())[0 + (4 * 2)] = errorFunc3.GetPartialDerivative(args, 11);
		(*J.GetPtrMatriz())[0 + (4 * 3)] = errorFunc4.GetPartialDerivative(args, 11);

		(*J.GetPtrMatriz())[1 + (4 * 0)] = errorFunc1.GetPartialDerivative(args, 8);
		(*J.GetPtrMatriz())[1 + (4 * 1)] = errorFunc2.GetPartialDerivative(args, 8);
		(*J.GetPtrMatriz())[1 + (4 * 2)] = errorFunc3.GetPartialDerivative(args, 8);
		(*J.GetPtrMatriz())[1 + (4 * 3)] = errorFunc4.GetPartialDerivative(args, 8);

		(*J.GetPtrMatriz())[2 + (4 * 0)] = errorFunc1.GetPartialDerivative(args, 9);
		(*J.GetPtrMatriz())[2 + (4 * 1)] = errorFunc2.GetPartialDerivative(args, 9);
		(*J.GetPtrMatriz())[2 + (4 * 2)] = errorFunc3.GetPartialDerivative(args, 9);
		(*J.GetPtrMatriz())[2 + (4 * 3)] = errorFunc4.GetPartialDerivative(args, 9);

		(*J.GetPtrMatriz())[3 + (4 * 0)] = errorFunc1.GetPartialDerivative(args, 2);
		(*J.GetPtrMatriz())[3 + (4 * 1)] = errorFunc2.GetPartialDerivative(args, 2);
		(*J.GetPtrMatriz())[3 + (4 * 2)] = errorFunc3.GetPartialDerivative(args, 2);
		(*J.GetPtrMatriz())[3 + (4 * 3)] = errorFunc4.GetPartialDerivative(args, 2);

		Matriz JInv = CalcInvMatriz(J);

		Matriz F({errorFunc1(args), errorFunc2(args), 
		errorFunc3(args), errorFunc4(args)}, 4, 1);

		norma = CalculaNormaVetor(F); 

		Xp = X - (JInv*F);

		X = Xp;

		std::cout << "-----------------------";	
		std::cout << n << std::endl;	
		std::cout << "Matriz X "; 
		MostrarMatriz(X);
		std::cout << std::endl;
		std::cout << "Matriz F ";
		MostrarMatriz(F);
		std::cout << std::endl;
		std::cout << "Matriz J ";
		MostrarMatriz(J);
		std::cout << std::endl;
		std::cout << "Norma matriz F: " << norma << std::endl;
		std::cout << "-----------------------" <<std::endl;
		n++;
	}

	//MostrarMatriz(Xp);
	//std::cout << "Norma matriz F: " << norma << std::endl;

    return 0;
}