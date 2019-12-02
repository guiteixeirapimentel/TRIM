#ifndef LIFTCOEFFICIENT_H
#define LIFTCOEFFICIENT_H

class LiftCoefficient
{
public:
    LiftCoefficient(double CL0, double CLAlfa, double CLAlfadot, double CLq, double CLdeltaE)
    :
    cCL0(CL0),
    cCLAlfa(CLAlfa),
    cCLAlfadot(CLAlfadot),
    cCLQ(CLq),
    cCLDeltaE(CLdeltaE){}
    ~LiftCoefficient(){}

    inline double operator() (double alfa, double alfadothat, double qhat, double deltae) const
    {
        return cCL0 + (alfa*cCLAlfa)+(cCLAlfadot*alfadothat) + (cCLQ*qhat)+(cCLDeltaE*deltae);
    }

private:
    const double cCL0;
    const double cCLAlfa;
    const double cCLAlfadot;
    const double cCLQ;
    const double cCLDeltaE;
};

#endif