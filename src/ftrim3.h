#ifndef F3_H
#define F3_H
#include "YMoment.h"
#include "State.h"
#include "StateDot.h"

#include "Alfa.h"
#include "AlfaDot.h"

class ftrim3 : public MultiVarFunction
{
public:
    ftrim3(YMoment M, double rho, double g, double S, double meanchord, double momentInertia)
    :
    MultiVarFunction(13),cM(M), cRho(rho), cG(g), cS(S), cMeanChord(meanchord), cMomentInertia(momentInertia){}
    ~ftrim3(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT; 10->vel; 11->alfa; 12->alfadot;
    inline double operator() (const std::vector<double>& args) const override
    {
        const double V = args[10];
        const double dynPress = 0.5 * cRho * (V*V);

        const double alfa = args[11];
        const double alfadothat = args[12]*cMeanChord/(2.0*V);
        const double qhat = args[2]*cMeanChord/(2.0*V);

        return args[6] - 
        (cM({alfa, alfadothat, qhat, args[8], dynPress, cS, cMeanChord})/cMomentInertia);
    }

private:
    const YMoment cM;
    const double cRho;
    const double cG;
    const double cS;
    const double cMeanChord;
    const double cMomentInertia;
};

#endif