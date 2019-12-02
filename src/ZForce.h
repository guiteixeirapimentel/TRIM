#ifndef ZFORCE_H
#define ZFORCE_H
#include "AerodynamicCoefficients/LiftCoefficient.h"
#include "AerodynamicCoefficients/DragCoefficient.h"
#include <cmath>

class ZForce
{
public:
    ZForce(LiftCoefficient liftCoeff, DragCoefficient dragCoeff)
    :
    cCL(liftCoeff),cCD(dragCoeff){}
    ~ZForce(){}

    inline double operator() (double alfa, double alfadothat, double qhat, double deltaE,
     double dynamicPressure, double S) const
    {
        return dynamicPressure*S*(
            (-cCL(alfa, alfadothat, qhat, deltaE)*cos(alfa)) - (cCD(alfa)*sin(alfa)));
    }

private:
    const LiftCoefficient cCL;
    const DragCoefficient cCD;
};

#endif