#ifndef ELEVATORFUNCTION_H
#define ELEVATORFUNCTION_H
#define _USE_MATH_DEFINES_
#define _USE_MATH_DEFINES
#include <cmath>

class ElevatorFunction
{
public:
    ElevatorFunction(){}
    ~ElevatorFunction(){};

    double operator() (double t) const
    {
        if(t >= 1.0 && t <= 2.0)
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