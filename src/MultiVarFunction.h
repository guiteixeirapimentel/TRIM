#ifndef MULTIVARFUNCTION_H
#define MULTIVARFUNCTION_H
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>

class MultiVarFunction
{
public:
    MultiVarFunction(size_t nArgs): cNArgs(nArgs){};
    ~MultiVarFunction() {};

    double virtual operator() (const std::vector<double>& args) const = 0;

    double GetPartialDerivative(const std::vector<double>& ref, size_t argId, 
        const double precision = 1e-7, double firsth = 1e-7)
    {
        double err = 10000.0;
        double h = firsth;

        double res;
        while(err > precision)
        {
            std::vector<double> refmthree = ref;
            refmthree[argId] -= 3.0*h; 

            std::vector<double> refmtwo = ref;
            refmtwo[argId] -= 2.0*h; 

            std::vector<double> refmone = ref;
            refmone[argId] -= h; 

            std::vector<double> refpone = ref;
            refpone[argId] += h; 

            std::vector<double> refptwo = ref;
            refptwo[argId] += 2.0*h; 

            std::vector<double> refpthree = ref;
            refpthree[argId] += 3.0*h; 

            const double fourPoints = ((this->operator()(refmtwo)) - 
            (8.0*this->operator()(refmone))
            + (8.0 * this->operator()(refpone))-(this->operator()(refptwo))
            ) / (12.0 * h);

            const double onePoint = (this->operator()(ref) - this->operator()(refmone))/h;

            const double sixPoints = ((-this->operator()(refmthree))+
            (9.0*this->operator()(refmtwo)) - (45.0*this->operator()(refmone))
            + (45.0 * this->operator()(refpone))-(9.0*this->operator()(refptwo))
            +(this->operator()(refpthree))) / (60.0 * h);
 
            err = fabs((fourPoints - sixPoints)/sixPoints);

                        
            res = sixPoints;

            h *= 1e-1;
            if(onePoint == 0.0)
            {
                err = 0.0;
                res = 0.0;
            }
        }

        return res;
    }
public:
    const size_t cNArgs;
};

#endif