#ifndef STATEDOT_H
#define STATEDOT_H

class StateDot
{
public:
    StateDot(double udot, double wdot, double qdot, double thetadot)
    :udot(udot), wdot(wdot), qdot(qdot), thetadot(thetadot){};
    ~StateDot(){};

public:
    double udot;
    double wdot;
    double qdot;
    double thetadot;
};

#endif