#include <iostream>
#include <vector>
#include <fstream>
#include "ElevatorFunction.h"

#include "ftrim1.h"
#include "ftrim2.h"
#include "ftrim3.h"
#include "ftrim4.h"

#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"


#include "Math/Matriz.h"
#include "Math/MathUtils.h"

#define _DEBUG

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

	f1 errorFunc1(xf, T, rho, gravity, S, meanChord, mass);
	f2 errorFunc2(zf, rho, gravity, S, meanChord, mass);
	f3 errorFunc3(ym, rho, gravity, S, meanChord, momentInertia);
	f4 errorFunc4;

	ftrim1 errorFuncTrim1(xf, T, rho, gravity, S, meanChord, mass);
	ftrim2 errorFuncTrim2(zf, rho, gravity, S, meanChord, mass);
	ftrim3 errorFuncTrim3(ym, rho, gravity, S, meanChord, momentInertia);
	ftrim4 errorFuncTrim4;

	std::vector<std::pair<double, Matriz>> condTrim;
	std::ofstream arqSaid("condtrim.csv");

	arqSaid << "TRIMAGEM AERONAVE\n";
	arqSaid << "V;gama;alfa;deltaE;deltaT;q\n";

	for(size_t i = 0; i < 58; i++)
	{
		double V = 33.0 + static_cast<double>(i);

		// 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    	// 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
		std::vector<double> args = {V, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, V, 0.0, 0.0};

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

			(*J.GetPtrMatriz())[0 + (4 * 0)] = errorFuncTrim1.GetPartialDerivative(args, 11);
			(*J.GetPtrMatriz())[0 + (4 * 1)] = errorFuncTrim2.GetPartialDerivative(args, 11);
			(*J.GetPtrMatriz())[0 + (4 * 2)] = errorFuncTrim3.GetPartialDerivative(args, 11);
			(*J.GetPtrMatriz())[0 + (4 * 3)] = errorFuncTrim4.GetPartialDerivative(args, 11);

			(*J.GetPtrMatriz())[1 + (4 * 0)] = errorFuncTrim1.GetPartialDerivative(args, 8);
			(*J.GetPtrMatriz())[1 + (4 * 1)] = errorFuncTrim2.GetPartialDerivative(args, 8);
			(*J.GetPtrMatriz())[1 + (4 * 2)] = errorFuncTrim3.GetPartialDerivative(args, 8);
			(*J.GetPtrMatriz())[1 + (4 * 3)] = errorFuncTrim4.GetPartialDerivative(args, 8);

			(*J.GetPtrMatriz())[2 + (4 * 0)] = errorFuncTrim1.GetPartialDerivative(args, 9);
			(*J.GetPtrMatriz())[2 + (4 * 1)] = errorFuncTrim2.GetPartialDerivative(args, 9);
			(*J.GetPtrMatriz())[2 + (4 * 2)] = errorFuncTrim3.GetPartialDerivative(args, 9);
			(*J.GetPtrMatriz())[2 + (4 * 3)] = errorFuncTrim4.GetPartialDerivative(args, 9);

			(*J.GetPtrMatriz())[3 + (4 * 0)] = errorFuncTrim1.GetPartialDerivative(args, 2);
			(*J.GetPtrMatriz())[3 + (4 * 1)] = errorFuncTrim2.GetPartialDerivative(args, 2);
			(*J.GetPtrMatriz())[3 + (4 * 2)] = errorFuncTrim3.GetPartialDerivative(args, 2);
			(*J.GetPtrMatriz())[3 + (4 * 3)] = errorFuncTrim4.GetPartialDerivative(args, 2);

			Matriz JInv = CalcInvMatriz(J);

			Matriz F({errorFuncTrim1(args), errorFuncTrim2(args), 
			errorFuncTrim3(args), errorFuncTrim4(args)}, 4, 1);

			norma = CalculaNormaVetor(F); 

			Xp = X - (JInv*F);

			X = Xp;

#ifdef _DEBUG
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
#endif
			n++;
		}

		arqSaid << V <<";"<<0.0 <<";" << X.cMatriz[0]<<";"
		<<X.cMatriz[1]<<";"<<X.cMatriz[2]<<";"<<X.cMatriz[3] <<"\n";

		condTrim.push_back({V, X});		
	}
	
	arqSaid.close();

	std::ofstream arqSaidMatrizes("matrizesLin.csv");
	// 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
	std::vector<double> args = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
	
	for(size_t i = 0; i < condTrim.size(); i++)
	{	
		args[11] = condTrim[i].second.cMatriz[0];
		args[8] = condTrim[i].second.cMatriz[1];
		args[9] = condTrim[i].second.cMatriz[2];
		args[2] = condTrim[i].second.cMatriz[3];

		args[10] = condTrim[i].first;

		args[1] = tan(args[11])*args[10];
		args[0] = sqrt((args[10]*args[10])-((args[1]*args[1])));
		args[3] = args[11];

		Matriz mEE = MatrizZeros(4,4);
		Matriz mAA = MatrizZeros(4,4);
		Matriz mBB = MatrizZeros(4, 2);

		for(size_t n = 0; n < 4; n++)
		{
			mEE.cMatriz[n + (4*0)] = -errorFunc1.GetPartialDerivative(args, 4+n);
			mEE.cMatriz[n + (4*1)] = -errorFunc2.GetPartialDerivative(args, 4+n);
			mEE.cMatriz[n + (4*2)] = -errorFunc3.GetPartialDerivative(args, 4+n);
			mEE.cMatriz[n + (4*3)] = -errorFunc4.GetPartialDerivative(args, 4+n);

			mAA.cMatriz[n + (4*0)] = errorFunc1.GetPartialDerivative(args, n);
			mAA.cMatriz[n + (4*1)] = errorFunc2.GetPartialDerivative(args, n);
			mAA.cMatriz[n + (4*2)] = errorFunc3.GetPartialDerivative(args, n);
			mAA.cMatriz[n + (4*3)] = errorFunc4.GetPartialDerivative(args, n);
		}

		for(size_t n = 0; n < 2; n++)
		{
			mBB.cMatriz[n + (2*0)] = errorFunc1.GetPartialDerivative(args, n+8);
			mBB.cMatriz[n + (2*1)] = errorFunc2.GetPartialDerivative(args, n+8);
			mBB.cMatriz[n + (2*2)] = errorFunc3.GetPartialDerivative(args, n+8);
			mBB.cMatriz[n + (2*3)] = errorFunc4.GetPartialDerivative(args, n+8);
		}

		Matriz Einv = CalcInvMatriz(mEE);

		Matriz A = Einv * mAA;
		Matriz B = Einv * mBB;

#ifdef _DEBUG
		Matriz I = Einv * mEE;
		MostrarMatriz(I);

		MostrarMatriz(mAA);
		MostrarMatriz(mBB);
#endif
		for(size_t n = 0; n < 4; n++)
		{
			arqSaidMatrizes << A.cMatriz[0+(n*4)] << ";" << A.cMatriz[1+(n*4)] << ";" 
			<< A.cMatriz[2+(n*4)] << ";" << A.cMatriz[3+(n*4)] << ";";

			arqSaidMatrizes << B.cMatriz[0+(n*4)] << ";" << B.cMatriz[1+(n*4)] << ";" 
			<< ";";

			arqSaidMatrizes << "\n";
		}
	}

	arqSaidMatrizes.close();

    return 0;
}