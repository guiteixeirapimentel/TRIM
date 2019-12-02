#ifndef ALFADOT_H
#define ALFADOT_H
#include <cmath>

class AlfaDot
{
public:
    AlfaDot(){};
    ~AlfaDot(){};

    double operator() (double u, double w, double udot, double wdot) const
    {
        const double V = sqrt((u*u)+(w*w));

        return ((wdot/V) - ((w*((u*udot) + (w*wdot))/pow((u*u)+(w*w), 1.5)))) / 
        (1 + ((w*w)/((w*w)+(u*u))));
    }
};

#endif