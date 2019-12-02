#ifndef F2_H
#define F2_H
#include "ZForce.h"
#include "State.h"
#include "StateDot.h"

#include "Alfa.h"
#include "AlfaDot.h"

class f2
{
public:
    f2(ZForce Z, double rho, double g, double S)
    :
    cZ(Z), cRho(rho), cG(g), cS(S){}
    ~f2(){}

    inline double operator() (double meanchord,
    const State& X, const StateDot& Xdot, double mass, double deltaE) const
    {
        const double V = sqrt((X.u*X.u)+(X.w*X.w));
        const double dynPress = 0.5 * cRho * (V*V);

        return Xdot.wdot - ((1/mass)*
        cZ(cAlfa(X.u, X.w), cAlfaDot(X.u, X.w, Xdot.udot, Xdot.wdot)*meanchord/(2.0*V),
        X.q*meanchord/(2.0*V), deltaE, dynPress, cS))
        + (cG * cos(X.theta)) - (X.q*X.u);

    }

    double JacobianX1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.u = Xref.u - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, X, Xdotref, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX2(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.w = Xref.w - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, X, Xdotref, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX3(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.q = Xref.q - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, X, Xdotref, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX4(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.theta = Xref.theta - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, X, Xdotref, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.udot = Xdotref.udot - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, Xref, Xdot, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot2(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.wdot = Xdotref.wdot - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, Xref, Xdot, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot3(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.qdot = Xdotref.qdot - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, Xref, Xdot, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot4(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.thetadot = Xdotref.thetadot - delta;
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, Xref, Xdot, massref, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianU1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double massref, double deltaEref) const
    {
        double delta = 1e-10;        
        double res = operator()(meanchord, Xref, Xdotref, massref, deltaEref)- 
        operator()(meanchord, Xref, Xdotref, massref, deltaEref - delta);
        res /= delta;  

        return res;
    }
    
private:
    const ZForce cZ;
    const Alfa cAlfa;
    const AlfaDot cAlfaDot;
    const double cRho;
    const double cG;
    const double cS;
};

#endif