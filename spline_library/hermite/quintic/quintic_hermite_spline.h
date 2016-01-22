#pragma once

#include <unordered_map>
#include <cassert>

#include "../../spline.h"
#include "quintic_hermite_spline_common.h"

#include "../../utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class QuinticHermiteSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    QuinticHermiteSpline(const std::vector<InterpolationType> &points,
                         const std::vector<InterpolationType> &tangents,
                         const std::vector<InterpolationType> &curvatures,
                         floating_t alpha = 0.0
                         );
    QuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0f);

//methods
public:
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t arcLength(floating_t a, floating_t b) const override { if(a > b) std::swap(a,b); return common.getLength(a,b); }
    floating_t totalLength(void) const override { return common.getTotalLength(); }

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

//data
private:
    QuinticHermiteSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
QuinticHermiteSpline<InterpolationType,floating_t>::QuinticHermiteSpline(
        const std::vector<InterpolationType> &points,
        const std::vector<InterpolationType> &tangents,
        const std::vector<InterpolationType> &curvatures,
        floating_t alpha
       )
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    int size = points.size();
    int padding = 0;
    int numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, padding);
    maxT = indexToT.at(numSegments);

    //pre-arrange the data needed for interpolation
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
    for(int i = 0; i < numSegments + 1; i++)
    {
        knots[i] = indexToT.at(i);

        positionData[i].position = points.at(i);
        positionData[i].tangent = tangents.at(i);
        positionData[i].curvature = curvatures.at(i);
    }

    common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
}

template<class InterpolationType, typename floating_t>
QuinticHermiteSpline<InterpolationType,floating_t>::QuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 6);

    int size = points.size();
    int firstTangent = 1;
    int lastTangent = size - 2;
    int firstCurvature = 2;
    int numSegments = size - 5;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, firstCurvature);
    maxT = indexToT.at(firstCurvature + numSegments);


    //compute the tangents
    std::vector<InterpolationType> tangents(size);
    for(int i = firstTangent; i <= lastTangent; i++)
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

    //compute the curvatures
    std::vector<InterpolationType> curves(size);
    int lastCurvature = firstCurvature + numSegments;
    for(int i = firstCurvature; i <= lastCurvature; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = tangents.at(i - 1);
        InterpolationType pCurrent = tangents.at(i);
        InterpolationType pNext = tangents.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        curves[i] =
                  pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

             //plus a little something extra - this is derived from the pyramid contruction
             //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
             //yielding the standard catmull-rom formula
                - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }


    //pre-arrange the data needed for interpolation
    int lastSegment = firstCurvature + numSegments;
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
    for(int i = firstCurvature; i < lastSegment + 1; i++)
    {
        knots[i - firstCurvature] = indexToT.at(i);

        positionData[i - firstCurvature].position = points.at(i);
        positionData[i - firstCurvature].tangent = tangents.at(i);
        positionData[i - firstCurvature].curvature = curves.at(i);
    }
    common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
}
