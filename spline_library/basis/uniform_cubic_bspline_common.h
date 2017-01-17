#pragma once

#include <vector>

#include "../utils/calculus.h"

template<class InterpolationType, typename floating_t>
class UniformCubicBSplineCommon
{
public:

    inline UniformCubicBSplineCommon(void) = default;
    inline UniformCubicBSplineCommon(std::vector<InterpolationType> points)
        :points(std::move(points))
    {}

    inline size_t segmentCount(void) const
    {
        return points.size() - 3;
    }

    inline size_t segmentForT(floating_t t) const
    {
        size_t segmentIndex = size_t(t);
        if(segmentIndex > segmentCount() - 1)
            return segmentCount() - 1;
        else
            return segmentIndex;
    }

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return computePosition(knotIndex, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT),
                    computeCurvature(knotIndex, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT),
                    computeCurvature(knotIndex, localT),
                    computeWiggle(knotIndex)
                    );
    }

    inline floating_t segmentLength(size_t index, floating_t from, floating_t to) const
    {
        auto segmentFunction = [this, index](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, t);
            return tangent.length();
        };

        return SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, from, to);
    }

private: //methods
    inline InterpolationType computePosition(size_t index, floating_t t) const
    {
        return (
                    points[index] * ((1 - t) * (1 - t) * (1 - t)) +
                    points[index + 1] * (t * t * 3 * (t - 2) + 4) +
                    points[index + 2] * (t * (t * (-3 * t + 3) + 3) + 1) +
                    points[index + 3] * (t * t * t)
                ) / floating_t(6);
    }

    inline InterpolationType computeTangent(size_t index, floating_t t) const
    {
        return (
                    points[index] * (-(1 - t) * (1 - t)) +
                    points[index + 1] * (t * (3 * t - 4)) +
                    points[index + 2] * ((3 * t + 1) * (1 - t)) +
                    points[index + 3] * (t * t)
                ) / floating_t(2);
    }

    inline InterpolationType computeCurvature(size_t index, floating_t t) const
    {
        return (
                    points[index] * (1 - t) +
                    points[index + 1] * (3 * t - 2) +
                    points[index + 2] * (1 - 3 * t) +
                    points[index + 3] * (t)
                );
    }

    inline InterpolationType computeWiggle(size_t index) const
    {
        return floating_t(3) * (points[index + 1] - points[index + 2]) + (points[index + 3] - points[index]);
    }

private: //data
    std::vector<InterpolationType> points;
};
