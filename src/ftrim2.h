#ifndef FTRIM2_H
#define FTRIM2_H

#include "ZForce.h"


class ftrim2 : public MultiVarFunction
{
public:
    ftrim2(ZForce Z, double rho, double g, double S, double meanChord, double mass)
    :
    MultiVarFunction(13), cZ(Z), cRho(rho), cG(g), cS(S), cMeanChord(meanChord), 
    cMass(mass){};
    ~ftrim2(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
    inline double operator() (const std::vector<double>& args) const override
    {
        const double V = args[10];
        const double dynPress = 0.5 * cRho * (V*V);

        const double alfadothat = args[12]*cMeanChord/(2.0*V);
        const double qhat = args[2]*cMeanChord/(2.0*V);

        return args[5] - ((1/cMass)*
        (cZ({args[11], alfadothat,
        qhat, args[8], dynPress, cS })))
            - (cG * cos(args[3])) - (args[2]*args[0]);

    }
    
private:
    const ZForce cZ;
    const double cRho;
    const double cG;
    const double cS;
    const double cMeanChord;
    const double cMass;
};

#endif