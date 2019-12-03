#ifndef F4_H
#define F4_H
#include "MultiVarFunction.h"

class f4 : public MultiVarFunction
{
public:
    f4():MultiVarFunction(13){};
    ~f4(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT
    inline double operator() (const std::vector<double>& args) const override
    {
        return args[7] - args[2];
    }
};

#endif