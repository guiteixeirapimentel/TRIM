#ifndef XFORCE_H
#define XFORCE_H
#include "AerodynamicCoefficients/LiftCoefficient.h"
#include "AerodynamicCoefficients/DragCoefficient.h"
#include <cmath>

class XForce
{
public:
    XForce(LiftCoefficient liftCoeff, DragCoefficient dragCoeff)
    :
    cCL(liftCoeff),cCD(dragCoeff){}
    ~XForce(){}

    inline double operator() (double alfa, double alfadothat, double qhat, double deltaE,
     double dynamicPressure, double S) const
    {
        return dynamicPressure*S*(
            (cCL(alfa, alfadothat, qhat, deltaE)*sin(alfa)) - (cCD(alfa)*cos(alfa)));
    }

private:
    const LiftCoefficient cCL;
    const DragCoefficient cCD;
};

#endif