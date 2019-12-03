#ifndef FTRIM4_H
#define FTRIM4_H
#include "State.h"
#include "StateDot.h"
#include "MultiVarFunction.h"

class ftrim4 : public MultiVarFunction
{
public:
    ftrim4():MultiVarFunction(13){};
    ~ftrim4(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
    inline double operator() (const std::vector<double>& args) const override
    {
        return args[7] - args[2];
    }
};

#endif