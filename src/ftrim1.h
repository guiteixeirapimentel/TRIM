#ifndef F1_H
#define F1_H
#include "XForce.h"
#include "ThrustForce.h"
#include "State.h"
#include "StateDot.h"

#include "Alfa.h"
#include "AlfaDot.h"

class ftrim1 : public MultiVarFunction
{
public:
    ftrim1(XForce X, ThrustForce T, double rho, double g, double S, double meanChord, double mass)
    :
    MultiVarFunction(13),cX(X), cThrust(T), cRho(rho), cG(g), 
    cS(S),cMeanChord(meanChord), cMass(mass){};
    ~ftrim1(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
    inline double operator() (const std::vector<double>& args) const override
    {
        const double V = args[10];
        const double dynPress = 0.5 * cRho * (V*V);

        const double alfa = args[11];
        const double alfadothat = args[12]*cMeanChord/(2.0*V);
        const double qhat = args[2]*cMeanChord/(2.0*V);

        return args[4] - ((1/cMass)*
        (cX({ alfa, alfadothat, qhat, args[8], dynPress, cS }) + 
        cThrust({V, args[9]}))) + (cG * sin(args[3])) + (args[2]*args[1]);

    }

private:
    const XForce cX;
    const ThrustForce cThrust;
    const double cRho;
    const double cG;
    const double cS;
    const double cMeanChord;
    const double cMass;
};

#endif