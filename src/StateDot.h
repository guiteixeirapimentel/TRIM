#ifndef STATEDOT_H
#define STATEDOT_H

class StateDot
{
public:
    StateDot(double udot, double wdot, double qdot, double thetadot)
    {cStatesDot.resize(4); cStatesDot[0] = udot; cStatesDot[1] = wdot; 
    cStatesDot[2] = qdot; cStatesDot[3] = thetadot;};
    ~StateDot(){};

public:
    // 0->udot; 1->wdot; 2->qdot; 3->thetadot;
    std::vector<double> cStatesDot;
};

#endif