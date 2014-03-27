#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <functional>

class Optimization
{
public:
    //Brent's method: Find 'x' between a and b for continuous function 'f' such that f(x) = 0
    //in other words, return a point between a and b where the given function crosses the x axis
    //fa must be equal to f(a), fb must be equal to f(b), fa and fb must have opposite signs
    //the algorithm returns when abs(f(x)) is less than 'tolerance'
    static double brentsMethod(const std::function<double (double)> &f, double a, double fa, double b, double fb, double tolerance = .0001);
private:
    Optimization();
};

#endif // OPTIMIZATION_H
