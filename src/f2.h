#ifndef F2_H
#define F2_H
#include "ZForce.h"

#include "AlfaDot.h"


class f2 : public MultiVarFunction
{
public:
    f2(ZForce Z, double rho, double g, double S, double meanChord, double mass)
    :
    MultiVarFunction(13), cZ(Z), cRho(rho), cG(g), cS(S), cMeanChord(meanChord), 
    cMass(mass){};
    ~f2(){};

    // 0->u; 1->w; 2->q; 3->theta; 4->udot; 5->wdot; 6->qdot; 7->thetadot; 
    // 8->deltaE; 9->deltaT;
    inline double operator() (const std::vector<double>& args) const override
    {
        const double V = sqrt((args[0]*args[0])+(args[1]*args[1]));
        const double dynPress = 0.5 * cRho * (V*V);

        const double alfadothat = cAlfadot({args[0], args[1], args[4], args[5]})*cMeanChord/(2.0*V);
        const double qhat = args[2]*cMeanChord/(2.0*V);

        return args[5] - ((1/cMass)*
        (cZ({args[11], alfadothat,
        qhat, args[8], dynPress, cS })))
            - (cG * cos(args[3])) - (args[2]*args[0]);

    }
    
private:
    const ZForce cZ;
    const AlfaDot cAlfadot;
    const double cRho;
    const double cG;
    const double cS;
    const double cMeanChord;
    const double cMass;
};

#endif