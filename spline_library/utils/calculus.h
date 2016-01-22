#pragma once

#include <cmath>

class SplineLibraryCalculus {
private:
    SplineLibraryCalculus() = default;

public:
    //numerically compute the integral of f from a to b, using the adaptive simpsons method
    template<class Function, typename floating_t>
    inline static floating_t adaptiveSimpsonsIntegral(Function f, floating_t a, floating_t b) {
        floating_t mid = (a+b)/2;
        return adaptiveSimpsonsHelper(f, a, f(a), b, f(b), mid, f(mid));
    }

private:
    //recursive helper method for adaptive simpsons method
    template<class Function, typename floating_t>
    static floating_t adaptiveSimpsonsHelper(Function f, floating_t a, floating_t aValue, floating_t b, floating_t bValue, floating_t mid, floating_t midValue);


    template<typename floating_t>
    inline static floating_t simpsons(floating_t a, floating_t aValue, floating_t b, floating_t bValue, floating_t midValue) {
        return ((b - a) / 6) * (aValue + 4 * midValue + bValue);
    }
};



template<class Function, typename floating_t>
floating_t SplineLibraryCalculus::adaptiveSimpsonsHelper(Function f,
                                                                floating_t a, floating_t aValue,
                                                                floating_t b, floating_t bValue,
                                                                floating_t mid, floating_t midValue)
{
    //compute integral for whole stretch
    floating_t mainIntegral = simpsons(a, aValue, b, bValue, midValue);

    //compute integral for left half
    floating_t leftMid = (a + mid)/2;
    floating_t leftMidValue = f(leftMid);
    floating_t leftIntegral = simpsons(a, aValue, mid, midValue, leftMidValue);

    //compute integral for right half
    floating_t rightMid = (b + mid)/2;
    floating_t rightMidValue = f(rightMid);
    floating_t rightIntegral = simpsons(mid, midValue, b, bValue, rightMidValue);

    //if mainIntegral is roughly equal to leftIntegral + rightIntegral, return
    floating_t sideIntegrals = leftIntegral + rightIntegral;
    if(std::abs((sideIntegrals - mainIntegral) / mainIntegral) < .001 || (b - a) < .001)
    {
        return sideIntegrals;
    }
    else
    {
        return adaptiveSimpsonsHelper(f, a, aValue, mid, midValue, leftMid, leftMidValue) +
               adaptiveSimpsonsHelper(f, mid, midValue, b, bValue, rightMid, rightMidValue);
    }
}
