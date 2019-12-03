#ifndef STATE_H
#define STATE_H
#include <vector>

class State
{
public:
    State(double u, double w, double q, double theta)
    {cStates.resize(4); cStates[0] = u; cStates[1] = w; 
    cStates[2] = q; cStates[3] = theta;};
    ~State(){};

public:
    // 0->u; 1->w; 2->q; 3->theta;
    std::vector<double> cStates;

};

#endif