#ifndef THRUSTFORCE_H
#define THRUSTFORCE_H
#include "MultiVarFunction.h"

class ThrustForce : public MultiVarFunction
{
public:
    ThrustForce(double valueMax):MultiVarFunction(2),cValueMax(valueMax){};
    ~ThrustForce(){};
    
    // 0 -> velocity; 1 -> deltaT
    inline double operator() (const std::vector<double>& args) const override
    {
        return (cValueMax* args[1])/args[0];
    }

private:
    const double cValueMax;
};

#endif