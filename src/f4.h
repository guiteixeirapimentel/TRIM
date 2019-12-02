#ifndef F4_H
#define F4_H
#include "State.h"
#include "StateDot.h"

class f4
{
public:
    f4(){};
    ~f4(){};

    inline double operator() (const State& X, const StateDot& Xdot) const
    {
        return Xdot.thetadot - X.q;
    }

    double JacobianX1(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.u = Xref.u - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(X, Xdotref);
        res /= delta;  

        return res;
    }

    double JacobianX2(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.w = Xref.w - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(X, Xdotref);
        res /= delta;  

        return res;
    }

    double JacobianX3(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.q = Xref.q - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(X, Xdotref);
        res /= delta;  

        return res;
    }

    double JacobianX4(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        State X = Xref;
        
        X.theta = Xref.theta - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(X, Xdotref);
        res /= delta;  

        return res;
    }

    double JacobianXdot1(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.udot = Xdotref.udot - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(Xref, Xdot);
        res /= delta;  

        return res;
    }

    double JacobianXdot2(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.wdot = Xdotref.wdot - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(Xref, Xdot);
        res /= delta;  

        return res;
    }

    double JacobianXdot3(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.qdot = Xdotref.qdot - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(Xref, Xdot);
        res /= delta;  

        return res;
    }

    double JacobianXdot4(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;
        StateDot Xdot = Xdotref;
        
        Xdot.thetadot = Xdotref.thetadot - delta;
        double res = operator()(Xref, Xdotref)- 
        operator()(Xref, Xdot);
        res /= delta;  

        return res;
    }

    double JacobianU1(const State& Xref, const StateDot& Xdotref) const
    {
        double delta = 1e-10;        
        double res = operator()(Xref, Xdotref)- 
        operator()(Xref, Xdotref);
        res /= delta;  

        return res;
    }
};

#endif