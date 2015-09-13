#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <algorithm>
#include <cmath>

#include "spline_library/utils/utils.h"

class Optimization
{
public:
    //Brent's method: Find 'x' between a and b for continuous function 'f' such that f(x) = 0
    //in other words, return a point between a and b where the given function crosses the x axis
    //fa must be equal to f(a), fb must be equal to f(b), fa and fb must have opposite signs
    //the algorithm returns when abs(f(x)) is less than 'tolerance'
    template<class Function, typename floating_t>
    static floating_t brentsMethod(Function f, floating_t a, floating_t fa, floating_t b, floating_t fb, floating_t tolerance = .0001);
private:
    Optimization() = default;
};

template<class Function, typename floating_t>
floating_t Optimization::brentsMethod(Function f, floating_t a, floating_t fa, floating_t b, floating_t fb, floating_t tolerance)
{
    floating_t currentGuess = a;
    floating_t currentGuessValue = fa;
    floating_t contrapoint = b;
    floating_t contrapointValue = fb;

    // http://en.wikipedia.org/wiki/Brent%27s_method
    // this algorithm is sort of a hybrid between the secant method and bisection method
    // the bisection method is too slow but is guaranteed to converge, the secant method is fast but not guaranteed to converge
    // this combines the good of both and leaves out the bad of both, but it's much more complicated to read

    //if a is a better guess than b, swap them - we want b to be better than a
    if(std::abs(contrapointValue) < std::abs(currentGuessValue))
    {
        std::swap(contrapoint, currentGuess);
        std::swap(contrapointValue, currentGuessValue);
    }

    bool mflag = true;
    floating_t prevGuess = contrapoint;

    floating_t prevGuessValue = contrapointValue;

    //this gets set at the end of the loop, but is only referenced if mflag is false
    //and mflag starts out true, so i can't be referenced while uninitialized
    floating_t oldGuess;

    floating_t minDelta = 0.001;

    while(std::abs(currentGuessValue) > tolerance && std::abs(currentGuess - contrapoint) > tolerance)
    {
        //s will be the next guess for the actual t value
        floating_t nextGuess;

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
                || (mflag && (std::abs(nextGuess - currentGuess) >= std::abs(currentGuess - prevGuess)/2)) //condition 2
                || (!mflag && (std::abs(nextGuess - currentGuess) >= std::abs(prevGuess - oldGuess)/2)) //condition 3
                || (mflag && (std::abs(currentGuess - prevGuess) < minDelta)) //condition 4
                || (!mflag && (std::abs(nextGuess - currentGuess) < minDelta)) //condition 5
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

        floating_t sValue = f(nextGuess);

        oldGuess = prevGuess;
        prevGuess = currentGuess;
        prevGuessValue = currentGuessValue;

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
        if(std::abs(contrapointValue) < std::abs(currentGuessValue))
        {
            std::swap(contrapoint, currentGuess);
            std::swap(contrapointValue, currentGuessValue);
        }
    }

    return currentGuess;
}

#endif // OPTIMIZATION_H
