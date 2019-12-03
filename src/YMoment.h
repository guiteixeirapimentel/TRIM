#ifndef YMOMENT_H
#define YMOMENT_H
#include "AerodynamicCoefficients/MomentCoefficient.h"


class YMoment : public MultiVarFunction
{
public:
    YMoment(MomentCoefficient momentCoeff)
    :
    MultiVarFunction(7),
    cCM(momentCoeff){};
    ~YMoment(){}
    // 0->alfa; 1->alfadothat; 2->qhat; 3->deltaE; 4->dynamicPressure; 5->S; 6->meanchord
    inline double operator() (const std::vector<double>& args) const override
    {
        return args[4]*args[5]*args[6]*cCM(args);
    }

private:
    const MomentCoefficient cCM;
};

#endif