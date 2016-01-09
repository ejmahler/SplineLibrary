#ifndef UNIFORM_CRSPLINE_COMMON_H
#define UNIFORM_CRSPLINE_COMMON_H


#include <vector>

#include "../../utils/spline_setup.h"

template<class InterpolationType, typename floating_t>
class UniformCRSplineCommon
{
public:

    inline UniformCRSplineCommon(void) = default;
    inline UniformCRSplineCommon(std::vector<InterpolationType> points)
        :points(std::move(points))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);
        floating_t localT = globalT - knotIndex;

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        return computePosition(knotIndex + 1, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);
        floating_t localT = globalT - knotIndex;

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);
        floating_t localT = globalT - knotIndex;

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT),
                    computeCurvature(knotIndex + 1, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);
        floating_t localT = globalT - knotIndex;

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT),
                    computeCurvature(knotIndex + 1, localT),
                    computeWiggle(knotIndex + 1)
                    );
    }


private: //methods
    inline InterpolationType computePosition(size_t index, floating_t t) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto oneMinusT = 1 - t;

        auto basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
        auto basis10 = t * oneMinusT * oneMinusT;

        auto basis11 = t * t * -oneMinusT;
        auto basis01 = t * t * (3 - 2*t);

        return
                basis00 * points[index] +
                basis10 * beforeTangent +

                basis11 * afterTangent +
                basis01 * points[index + 1];
    }

    inline InterpolationType computeTangent(size_t index, floating_t t) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto oneMinusT = 1 - t;

        auto d_basis00 = 6 * t * (t - 1);
        auto d_basis10 = (1 - 3*t) * oneMinusT;

        auto d_basis11 = t * (3 * t - 2);
        auto d_basis01 = -d_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the derivative of the position function and nothing else
        //if you know why please let me know
        return
                d_basis00 * points[index] +
                d_basis10 * beforeTangent +

                d_basis11 * afterTangent +
                d_basis01 * points[index + 1];
    }

    inline InterpolationType computeCurvature(size_t index, floating_t t) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto d2_basis00 = 6 * (2 * t - 1);
        auto d2_basis10 = 2 * (3 * t - 2);

        auto d2_basis11 = 2 * (3 * t - 1);
        auto d2_basis01 = -d2_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return
                d2_basis00 * points[index] +
                d2_basis10 * beforeTangent +

                d2_basis11 * afterTangent +
                d2_basis01 * points[index + 1];
    }

    inline InterpolationType computeWiggle(size_t index) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return 12 * (points[index] - points[index + 1]) + 6 * (beforeTangent + afterTangent);
    }

    inline InterpolationType computeTangentAtIndex(size_t i) const
    {
        return (points[i + 1] - points[i - 1]) / 2;
    }

private: //data
    std::vector<InterpolationType> points;
};

#endif // UNIFORM_CRSPLINE_COMMON_H

