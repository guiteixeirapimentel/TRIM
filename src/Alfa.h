#ifndef ALFA_H
#define ALFA_H
#include "MultiVarFunction.h"

class Alfa : public MultiVarFunction
{
public:
    Alfa():MultiVarFunction(2){};
    ~Alfa(){};

    // 0 -> u; 1 -> w;
    double operator() (const std::vector<double>& args) const override
    {
        return atan2(args[1], sqrt((args[1]*args[1])+(args[0]*args[0])));
    }
};

#endif