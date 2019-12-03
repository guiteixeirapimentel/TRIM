#ifndef LIFTCOEFFICIENT_H
#define LIFTCOEFFICIENT_H
#include "../MultiVarFunction.h"

class LiftCoefficient : public MultiVarFunction
{
public:
    LiftCoefficient(double CL0, double CLAlfa, double CLAlfadot, double CLq, double CLdeltaE)
    :
    MultiVarFunction(4),
    cCL0(CL0),
    cCLAlfa(CLAlfa),
    cCLAlfadot(CLAlfadot),
    cCLQ(CLq),
    cCLDeltaE(CLdeltaE){}
    ~LiftCoefficient(){}
    // 0 -> alfa; 1 -> alfdothat; 2 -> qhat; 3 -> deltae;
    inline double operator() (const std::vector<double>& args) const override
    {
        return cCL0 + (args[0]*cCLAlfa)+(cCLAlfadot*args[1]) + (cCLQ*args[2])+(cCLDeltaE*args[3]);
    }

private:
    const double cCL0;
    const double cCLAlfa;
    const double cCLAlfadot;
    const double cCLQ;
    const double cCLDeltaE;
};

#endif