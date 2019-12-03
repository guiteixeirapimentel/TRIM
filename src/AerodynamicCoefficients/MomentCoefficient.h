#ifndef MOMENTCOEFFICIENT_H
#define MOMENTCOEFFICIENT_H
#include "../MultiVarFunction.h"

class MomentCoefficient : public MultiVarFunction
{
public:
    MomentCoefficient(double CM0, double CMAlfa, double CMAlfadot, double CMq, double CMdeltaE)
    :
    MultiVarFunction(4),
    cCM0(CM0),
    cCMAlfa(CMAlfa),
    cCMAlfadot(CMAlfadot),
    cCMQ(CMq),
    cCMDeltaE(CMdeltaE){}

    ~MomentCoefficient(){}

    // 0 -> alfa; 1 -> alfdothat; 2 -> qhat; 3 -> deltae;
    inline double operator() (const std::vector<double>& args) const override
    {
        return cCM0 + (args[0]*cCMAlfa)+(cCMAlfadot*args[1]) + (cCMQ*args[2])+
        (cCMDeltaE*args[3]);
    }

private:
    const double cCM0;
    const double cCMAlfa;
    const double cCMAlfadot;
    const double cCMQ;
    const double cCMDeltaE;
};

#endif