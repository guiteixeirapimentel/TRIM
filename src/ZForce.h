#ifndef ZFORCE_H
#define ZFORCE_H
#include "AerodynamicCoefficients/LiftCoefficient.h"
#include "AerodynamicCoefficients/DragCoefficient.h"

class ZForce : public MultiVarFunction
{
public:
    ZForce(LiftCoefficient liftCoeff, DragCoefficient dragCoeff)
    :
    MultiVarFunction(6),cCL(liftCoeff),cCD(dragCoeff){}
    ~ZForce(){}

    // 0->alfa; 1->alfadothat; 2->qhat; 3->deltaE; 4->dynamicPressure; 5->S;
    inline double operator() (const std::vector<double>& args) const override
    {
        const double CL = cCL(args);
        const double CD = cCD({CL});

        return args[4]*args[5]*(
            (-CL*cos(args[0])) - (CD*sin(args[0])));
    }

private:
    const LiftCoefficient cCL;
    const DragCoefficient cCD;
};

#endif