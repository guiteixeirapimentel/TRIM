#ifndef F3_H
#define F3_H
#include "YMoment.h"

#include "Alfa.h"
#include "AlfaDot.h"

class f3 : public MultiVarFunction
{
public:
    f3(YMoment M, double rho, double g, double S, double meanchord, double momentInertia)
    :
    MultiVarFunction(13),cM(M), cRho(rho), cG(g), cS(S), cMeanChord(meanchord), cMomentInertia(momentInertia){}
    ~f3(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT
    inline double operator() (const std::vector<double>& args) const override
    {
        const double V = sqrt((args[0]*args[0])+(args[1]*args[1]));
        const double dynPress = 0.5 * cRho * (V*V);

        const double alfa = cAlfa(args);
        const double alfadothat = cAlfadot({args[0], args[1], args[4], args[5]})*cMeanChord/(2.0*V);
        const double qhat = args[2]*cMeanChord/(2.0*V);

        return args[6] - 
        (cM({alfa, alfadothat, qhat, args[8], dynPress, cS, cMeanChord})/cMomentInertia);
    }

private:
    const YMoment cM;
    const Alfa cAlfa;
    const AlfaDot cAlfadot;
    const double cRho;
    const double cG;
    const double cS;
    const double cMeanChord;
    const double cMomentInertia;
};

#endif