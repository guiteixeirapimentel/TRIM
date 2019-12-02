#ifndef ALFA_H
#define ALFA_H
#include <cmath>

class Alfa
{
public:
    Alfa(){};
    ~Alfa(){};

    double operator() (double u, double w) const
    {
        return atan2(w, sqrt((w*w)+(u*u)));
    }
};

#endif