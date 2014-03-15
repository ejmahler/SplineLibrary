#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <functional>

class Optimization
{
public:
    static double brentsMethod(const std::function<double (double)> &f, double a, double fa, double b, double fb, double tolerance = .0001);
private:
    Optimization();
};

#endif // OPTIMIZATION_H
