#ifndef DRAGCOEFFICIENT_H
#define DRAGCOEFFICIENT_H

class DragCoefficient
{
public:
    DragCoefficient(double CD0, double CDAlfa)
    :
    cCD0(CD0),
    cCDAlfa(CDAlfa){}
    ~DragCoefficient(){}

    inline double operator ()(double alfa) const
    {
        return cCD0 + (cCDAlfa*alfa);
    }

private:
    const double cCD0;
    const double cCDAlfa;
};

#endif