#ifndef ALFADOT_H
#define ALFADOT_H
#include "MultiVarFunction.h"

class AlfaDot : public MultiVarFunction
{
public:
    AlfaDot(): MultiVarFunction(4){};
    ~AlfaDot(){};

    // 0 -> u; 1 -> w; 2 -> udot; 3 -> wdot;
    double operator() (const std::vector<double>& args) const override
    {
        const double V = sqrt((args[0]*args[0])+(args[1]*args[1]));

        return ((args[3]/V) - ((args[1]*((args[0]*args[2]) + (args[1]*args[3]))/
        pow((args[0]*args[0])+(args[1]*args[1]), 1.5)))) / 
        (1 + ((args[1]*args[1])/((args[1]*args[1])+(args[0]*args[0]))));
    }
};

#endif