#pragma once

#include <cmath>
#include <array>

class SplineLibraryCalculus {
private:
    SplineLibraryCalculus() = default;

public:
    //use the gauss-legendre quadrature algorithm to numerically integrate f from a to b
    //as of this writing, hardcoded to use 13 points
    template<class Function, typename floating_t>
    inline static floating_t gaussLegendreQuadratureIntegral(Function f, floating_t a, floating_t b);

};

template<class Function, typename floating_t>
floating_t SplineLibraryCalculus::gaussLegendreQuadratureIntegral(Function f, floating_t a, floating_t b)
{
    const size_t NUM_POINTS = 13;

    //these are precomputed :( It would be cool to compute these at compile time, but apparently
    //it's not easy to compute the points/weights just given the number of points.
    //it involves computing every root of a polynomial. which can obviously be done, but not in a reasonable amount of code
    std::array<floating_t, NUM_POINTS> quadraturePoints = {
         0.0000000000000000,
        -0.2304583159551348,
         0.2304583159551348,
        -0.4484927510364469,
         0.4484927510364469,
        -0.6423493394403402,
         0.6423493394403402,
        -0.8015780907333099,
         0.8015780907333099,
        -0.9175983992229779,
         0.9175983992229779,
        -0.9841830547185881,
         0.9841830547185881
    };

    std::array<floating_t, NUM_POINTS> quadratureWeights = {
        0.2325515532308739,
        0.2262831802628972,
        0.2262831802628972,
        0.2078160475368885,
        0.2078160475368885,
        0.1781459807619457,
        0.1781459807619457,
        0.1388735102197872,
        0.1388735102197872,
        0.0921214998377285,
        0.0921214998377285,
        0.0404840047653159,
        0.0404840047653159
    };

    floating_t halfDiff = (b - a) / 2;
    floating_t halfSum = (a + b) / 2;

    floating_t sum(0);
    for(size_t i = 0; i < NUM_POINTS; i++)
    {
        sum += quadratureWeights[i] * f(halfDiff * quadraturePoints[i] + halfSum);
    }
    return halfDiff * sum;
}
