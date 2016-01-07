#ifndef CUBIC_HERMITE_SPLINE_COMMON
#define CUBIC_HERMITE_SPLINE_COMMON


#include <vector>

#include "../../utils/spline_setup.h"

template<class InterpolationType, typename floating_t>
class CubicHermiteSplineCommon
{
public:
    struct alignas(8) CubicHermiteSplinePoint
    {
        InterpolationType position, tangent;
    };

    inline CubicHermiteSplineCommon(void) = default;
    inline CubicHermiteSplineCommon(std::vector<CubicHermiteSplinePoint> points, std::vector<floating_t> knots)
        :points(std::move(points)), knots(std::move(knots))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return computePosition(knotIndex, localT, tDiff);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex, localT, tDiff),
                    computeTangent(knotIndex, localT, tDiff)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex, localT, tDiff),
                    computeTangent(knotIndex, localT, tDiff),
                    computeCurvature(knotIndex, localT, tDiff)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex, localT, tDiff),
                    computeTangent(knotIndex, localT, tDiff),
                    computeCurvature(knotIndex, localT, tDiff),
                    computeWiggle(knotIndex, tDiff)
                    );
    }


private: //methods
    inline InterpolationType computePosition(size_t index, floating_t t, floating_t tDiff) const
    {
        auto oneMinusT = 1 - t;

        auto basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
        auto basis10 = t * oneMinusT * oneMinusT;

        auto basis11 = t * t * -oneMinusT;
        auto basis01 = t * t * (3 - 2*t);

        return
                basis00 * points[index].position +
                basis10 * tDiff * points[index].tangent +

                basis11 * tDiff * points[index + 1].tangent +
                basis01 * points[index + 1].position;
    }

    inline InterpolationType computeTangent(size_t index, floating_t t, floating_t tDiff) const
    {
        auto oneMinusT = 1 - t;

        auto d_basis00 = 6 * t * (t - 1);
        auto d_basis10 = (1 - 3*t) * oneMinusT;

        auto d_basis11 = t * (3 * t - 2);
        auto d_basis01 = -d_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the derivative of the position function and nothing else
        //if you know why please let me know
        return (
                d_basis00 * points[index].position +
                d_basis10 * tDiff * points[index].tangent +

                d_basis11 * tDiff * points[index + 1].tangent +
                d_basis01 * points[index + 1].position
                ) / tDiff;
    }

    inline InterpolationType computeCurvature(size_t index, floating_t t, floating_t tDiff) const
    {
        auto d2_basis00 = 6 * (2 * t - 1);
        auto d2_basis10 = 2 * (3 * t - 2);

        auto d2_basis11 = 2 * (3 * t - 1);
        auto d2_basis01 = -d2_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                d2_basis00 * points[index].position +
                d2_basis10 * tDiff * points[index].tangent +

                d2_basis11 * tDiff * points[index + 1].tangent +
                d2_basis01 * points[index + 1].position
                ) / (tDiff * tDiff);
    }

    inline InterpolationType computeWiggle(size_t index, floating_t tDiff) const
    {
        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                    12 * (points[index].position - points[index + 1].position) + 6 * tDiff * (points[index].tangent + points[index + 1].tangent)
                ) / (tDiff * tDiff * tDiff);
    }

private: //data
    std::vector<CubicHermiteSplinePoint> points;
    std::vector<floating_t> knots;
};

#endif // CUBIC_HERMITE_SPLINE_COMMON

