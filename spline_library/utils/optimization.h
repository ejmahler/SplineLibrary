#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <algorithm>

#include "spline_library/utils/utils.h"

class Optimization
{
public:
    //Brent's method: Find 'x' between a and b for continuous function 'f' such that f(x) = 0
    //in other words, return a point between a and b where the given function crosses the x axis
    //fa must be equal to f(a), fb must be equal to f(b), fa and fb must have opposite signs
    //the algorithm returns when abs(f(x)) is less than 'tolerance'
    template<class Function>
    static double brentsMethod(Function f, double a, double fa, double b, double fb, double tolerance = .0001);
private:
    Optimization();
};

template<class Function>
double Optimization::brentsMethod(Function f, double a, double fa, double b, double fb, double tolerance)
{
    double currentGuess = a;
    double currentGuessValue = fa;
    double contrapoint = b;
    double contrapointValue = fb;

    // http://en.wikipedia.org/wiki/Brent%27s_method
    // this algorithm is sort of a hybrid between the secant method and bisection method
    // the bisection method is too slow but is guaranteed to converge, the secant method is fast but not guaranteed to converge
    // this combines the good of both and leaves out the bad of both, but it's much more complicated to read

    //if a is a better guess than b, swap them - we want b to be better than a
    if(abs(contrapointValue) < abs(currentGuessValue))
    {
        std::swap(contrapoint, currentGuess);
        std::swap(contrapointValue, currentGuessValue);
    }

    bool mflag = true;
    double prevGuess = contrapoint;

    double prevGuessValue = contrapointValue;

    //this gets set at the end of the loop, but is only referenced if mflag is false
    //and mflag starts out true, so i can't be referenced while uninitialized
    double oldGuess;

    double minDelta = 0.001;

    while(currentGuessValue > tolerance && abs(currentGuess - contrapoint) > tolerance)
    {
        //s will be the next guess for the actual t value
        double nextGuess;

        if(contrapointValue != prevGuessValue && currentGuessValue != prevGuessValue)
        {
            //if f(a),f(b) and f(c) are all distinct, use inverse quadractic interpolation
            nextGuess = contrapoint * currentGuessValue * prevGuessValue / ((contrapointValue - currentGuessValue) * (contrapointValue - prevGuessValue));
            nextGuess += currentGuess * contrapointValue * prevGuessValue / ((currentGuessValue - contrapointValue) * (currentGuessValue - prevGuessValue));
            nextGuess += prevGuess * contrapointValue * currentGuessValue / ((prevGuessValue - contrapointValue) * (prevGuessValue - currentGuessValue));
        }
        else
        {
            //use secant method
            nextGuess = currentGuess - currentGuessValue * (currentGuess - contrapoint) / (currentGuessValue - contrapointValue);
        }


        //determine if we can use the s that we jsut calculated. if not we have to use bisection method :(
        //condition numbers correspond to the pseudocode in the wiki article
        if(     (nextGuess < (3 * contrapoint + currentGuess) / 4 || nextGuess  > currentGuess) //condition 1
                || (mflag && (abs(nextGuess - currentGuess) >= abs(currentGuess - prevGuess)/2)) //condition 2
                || (!mflag && (abs(nextGuess - currentGuess) >= abs(prevGuess - oldGuess)/2)) //condition 3
                || (mflag && (abs(currentGuess - prevGuess) < minDelta)) //condition 4
                || (!mflag && (abs(nextGuess - currentGuess) < minDelta)) //condition 5
                )
        {
            //one of these was true so we have to use bisection
            nextGuess = (contrapoint + currentGuess) / 2;
            mflag = true;
        }
        else
        {
            mflag = false;
        }

        double sValue = f(nextGuess);

        oldGuess = prevGuess;
        prevGuess = currentGuess;

        //choose the next a: if a and s have the same sign then s is the new a
        //if a and s have different signs then we keep the existing a
        if(sign(sValue) == sign(contrapointValue))
        {
            contrapoint = nextGuess;
            contrapointValue = sValue;
        }
        else
        {
            currentGuess = nextGuess;
            currentGuessValue = sValue;
        }

        //if a is a better guess than b, swap them - we want b to be better than a
        if(abs(contrapointValue) < abs(currentGuessValue))
        {
            std::swap(contrapoint, currentGuess);
            std::swap(contrapointValue, currentGuessValue);
        }
    }

    return currentGuess;
}

#endif // OPTIMIZATION_H
