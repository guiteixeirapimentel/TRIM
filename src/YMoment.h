#ifndef YMOMENT_H
#define YMOMENT_H
#include "AerodynamicCoefficients/MomentCoefficient.h"
#include <cmath>

class YMoment
{
public:
    YMoment(MomentCoefficient momentCoeff)
    :
    cCM(momentCoeff){};
    ~YMoment(){}

    inline double operator() (double alfa, double alfadothat, double qhat, double deltaE,
     double dynamicPressure, double S, double meanchord) const
    {
        return dynamicPressure*S*meanchord*cCM(alfa, alfadothat, qhat, deltaE);
    }

private:
    const MomentCoefficient cCM;
};

#endif