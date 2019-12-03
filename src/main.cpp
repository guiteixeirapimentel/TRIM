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

	DragCoefficient CD(0.02, 0.07);
	LiftCoefficient CL(0.2977, 4.41, 1.7, 3.9, 0.43);
	MomentCoefficient CM(0.03917, -0.613, -7.27, -12.4, -1.122);

	XForce xf(CL, CD);
	ZForce zf(CL, CD);
	YMoment ym(CM);

	ThrustForce T(135e3);

	f1 errorFunc1(xf, T, rho, gravity, S);
	f2 errorFunc2(zf, rho, gravity, S);
	f3 errorFunc3(ym, rho, gravity, S);
	f4 errorFunc4;

	State Xref(70.0, 0.0, 0.0, 0.0);
	StateDot Xdotref(0.0, 0.0, 0.0, 0.0);
	double deltaEref = 0.0;

	


    return 0;
}