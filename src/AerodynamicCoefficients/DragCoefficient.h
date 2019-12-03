#ifndef DRAGCOEFFICIENT_H
#define DRAGCOEFFICIENT_H
#include "../MultiVarFunction.h"

class DragCoefficient : public MultiVarFunction
{
public:
    DragCoefficient(double CDmin, double k)
    :
    MultiVarFunction(1),
    cCDMin(CDmin),
    cK(k){};
    ~DragCoefficient(){};

    inline double operator ()(const std::vector<double>& CL) const override
    {
        return cCDMin + (cK*CL[0]*CL[0]);
    }

private:
    const double cCDMin;
    const double cK;
};

#endif