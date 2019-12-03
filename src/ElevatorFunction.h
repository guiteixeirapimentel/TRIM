#ifndef ELEVATORFUNCTION_H
#define ELEVATORFUNCTION_H
#define _USE_MATH_DEFINES_
#define _USE_MATH_DEFINES
#include "MultiVarFunction.h"

class ElevatorFunction : public MultiVarFunction
{
public:
    ElevatorFunction():MultiVarFunction(1){}
    ~ElevatorFunction(){};

    double operator() (const std::vector<double>& t) const override
    {
        if(t[0] >= 1.0 && t[0] <= 2.0)
        {
            return UMGRAU;
        }
        else
        {
            return DOISGRAUS;
        }
        
    }

private:
    static constexpr double UMGRAU = 1.0 * M_PI / 180.0;
    static constexpr double DOISGRAUS = 2.0 * M_PI / 180.0;
};

#endif