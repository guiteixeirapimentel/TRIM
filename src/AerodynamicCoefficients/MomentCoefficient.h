#ifndef MOMENTCOEFFICIENT_H
#define MOMENTCOEFFICIENT_H

class MomentCoefficient
{
public:
    MomentCoefficient(double CM0, double CMAlfa, double CMAlfadot, double CMq, double CMdeltaE)
    :
    cCM0(CM0),
    cCMAlfa(CMAlfa),
    cCMAlfadot(CMAlfadot),
    cCMQ(CMq),
    cCMDeltaE(CMdeltaE){}

    ~MomentCoefficient(){}

    inline double operator() (double alfa, double alfadothat, double qhat, double deltae) const
    {
        return cCM0 + (alfa*cCMAlfa)+(cCMAlfadot*alfadothat) + (cCMQ*qhat)+
        (cCMDeltaE*deltae);
    }

private:
    const double cCM0;
    const double cCMAlfa;
    const double cCMAlfadot;
    const double cCMQ;
    const double cCMDeltaE;
};

#endif