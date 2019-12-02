#ifndef THRUSTFORCE_H
#define THRUSTFORCE_H

class ThrustForce
{
public:
    ThrustForce(double value):cValue(value){}
    ~ThrustForce(){}

    inline double operator() () const
    {
        return cValue;
    }

private:
    const double cValue;
};

#endif