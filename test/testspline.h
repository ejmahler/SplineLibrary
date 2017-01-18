#pragma once

#include <array>
#include <QObject>
#include "spline_library/vector.h"

class TestSpline : public QObject
{
    Q_OBJECT
public:
    explicit TestSpline(QObject *parent = 0);

signals:

private slots:

    //test each spline's basic methods like maxT, segmentCount
    void testMethods_data(void);
    void testMethods(void);

    //use numeric integration to verify the first, second, and third derivatives
    void testDerivatives_data(void);
    void testDerivatives(void);

    //use numeric integration to verify the first, second, and third derivatives
    //this one is for splines without continuous curvature
    void testDerivativesNonC2_data(void);
    void testDerivativesNonC2(void);

private:
    //use the gauss-legendre quadrature algorithm to numerically integrate f from a to b
    //this one uses 7 points atm just because 7 is accurate enough for these unit tests
    //if it turns out to not be accurate enough, the point count can be upped
    //this is different from the one in calculus.h because we need this to return a vector2 for unit test purposes
    template<class Function>
    Vector2 gaussLegendreQuadratureIntegral(Function f, float a, float b) const;
};

template<class Function>
Vector2 TestSpline::gaussLegendreQuadratureIntegral(Function f, float a, float b) const
{
    std::array<float, 7> quadraturePoints = {
         0.0000000000000000f,
         0.4058451513773972f,
        -0.4058451513773972f,
         0.7415311855993945f,
        -0.7415311855993945f,
         0.9491079123427585f,
        -0.9491079123427585f
    };

    std::array<float, 7> quadratureWeights = {
        0.4179591836734694f,
        0.3818300505051189f,
        0.3818300505051189f,
        0.2797053914892766f,
        0.2797053914892766f,
        0.1294849661688697f,
        0.1294849661688697f
    };

    float halfDiff = (b - a) / 2;
    float halfSum = (a + b) / 2;

    Vector2 sum;
    for(int i = 0; i < 7; i++)
    {
        sum += quadratureWeights[i] * f(halfDiff * quadraturePoints[i] + halfSum);
    }
    return halfDiff  * sum;
}








