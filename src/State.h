#ifndef STATE_H
#define STATE_H

class State
{
public:
    State(double u, double w, double q, double theta)
    :u(u), w(w), q(q), theta(theta){};
    ~State(){};

public:
    double u;
    double w;
    double q;
    double theta;
};

#endif