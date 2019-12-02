#ifndef F3_H
#define F3_H
#include "YMoment.h"
#include "State.h"
#include "StateDot.h"

#include "Alfa.h"
#include "AlfaDot.h"

class f3
{
public:
    f3(YMoment M, double rho, double g, double S)
    :
    cM(M), cRho(rho), cG(g), cS(S){}
    ~f3(){};

    inline double operator() (double meanchord,
    const State& X, const StateDot& Xdot, double momentInertia, double deltaE) const
    {
        const double V = sqrt((X.u*X.u)+(X.w*X.w));
        const double dynPress = 0.5 * cRho * (V*V);

        return Xdot.qdot - 
        (cM(cAlfa(X.u, X.w), cAlfaDot(X.u, X.w, Xdot.udot, Xdot.wdot)*meanchord/(2.0*V),
        X.q*meanchord/(2.0*V), deltaE, dynPress, cS, meanchord)/momentInertia);
    }

    double JacobianX1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.u = Xref.u - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, X, Xdotref, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX2(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.w = Xref.w - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, X, Xdotref, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX3(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.q = Xref.q - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, X, Xdotref, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianX4(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.theta = Xref.theta - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, X, Xdotref, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.udot = Xdotref.udot - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, Xref, Xdot, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot2(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.wdot = Xdotref.wdot - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, Xref, Xdot, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot3(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.qdot = Xdotref.qdot - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, Xref, Xdot, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianXdot4(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.thetadot = Xdotref.thetadot - delta;
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, Xref, Xdot, momentInertia, deltaEref);
        res /= delta;  

        return res;
    }

    double JacobianU1(double meanchord,
    const State& Xref, const StateDot& Xdotref, double momentInertia, double deltaEref) const
    {
        double delta = 1e-10;        
        double res = operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref)- 
        operator()(meanchord, Xref, Xdotref, momentInertia, deltaEref - delta);
        res /= delta;  

        return res;
    }
private:
    const YMoment cM;
    const Alfa cAlfa;
    const AlfaDot cAlfaDot;
    const double cRho;
    const double cG;
    const double cS;
};

#endif