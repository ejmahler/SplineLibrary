#ifndef CUBICHERMITESPLINE_H
#define CUBICHERMITESPLINE_H

#include <unordered_map>
#include <cassert>

#include "spline_library/spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline_common.h"

#include "spline_library/utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class CubicHermiteSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    CubicHermiteSpline(const std::vector<InterpolationType> &points, const std::vector<InterpolationType> &tangents, floating_t alpha = 0.0);
    CubicHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0);

//methods
public:
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

//data
private:
    CubicHermiteSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
CubicHermiteSpline<InterpolationType,floating_t>::CubicHermiteSpline(
        const std::vector<InterpolationType> &points,
        const std::vector<InterpolationType> &tangents,
        floating_t alpha
        )
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    int size = points.size();
    int firstTangent = 0;
    int numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

    //pre-arrange the data needed for interpolation
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename CubicHermiteSplineCommon<InterpolationType, floating_t>::CubicHermiteSplineSegment> positionData(numSegments + 1);
    for(int i = 0; i < numSegments + 1; i++)
    {
        knots[i] = indexToT[i];
        points[i].position = points[i];
        points[i].tangent = tangents[i];
    }

    common = CubicHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
}

template<class InterpolationType, typename floating_t>
CubicHermiteSpline<InterpolationType,floating_t>::CubicHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 4);

    int size = points.size();
    int firstTangent = 1;
    int numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

    //compute the tangents
    std::vector<InterpolationType> tangents(size);
    for(int i = firstTangent; i < firstTangent + numSegments + 1; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = points.at(i - 1);
        InterpolationType pCurrent = points.at(i);
        InterpolationType pNext = points.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        tangents[i] =
                            pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                            + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                         //plus a little something extra - this is derived from the pyramid contruction
                         //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                         //yielding the standard catmull-rom formula
                            - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }

    //pre-arrange the data needed for interpolation
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename CubicHermiteSplineCommon<InterpolationType, floating_t>::CubicHermiteSplinePoint> positionData(numSegments + 1);
    for(int i = firstTangent; i < firstTangent + numSegments + 1; i++)
    {
        knots[i - firstTangent] = indexToT[i];
        positionData[i - firstTangent].position = points[i];
        positionData[i - firstTangent].tangent = tangents[i];
    }

    common = CubicHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
}

#endif // CUBICHERMITESPLINE_H
