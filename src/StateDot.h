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
    // 0->u; 1->w; 2->q; 3->theta;
    std::vector<double> cStatesDot;
};

#endif