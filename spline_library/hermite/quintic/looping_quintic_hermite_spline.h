#pragma once

#include <unordered_map>

#include "../../spline.h"
#include "quintic_hermite_spline_common.h"

#include "../../utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class LoopingQuinticHermiteSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points,
                                const std::vector<InterpolationType> &tangents,
                                const std::vector<InterpolationType> &curvatures,
                                floating_t alpha = 0.0
                                );

    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t arcLength(floating_t a, floating_t b) const override;
    floating_t totalLength(void) const override { return common.getTotalLength(); }

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return true; }

    size_t segmentCount(void) const override { return common.segmentCount(); }
    floating_t segmentT(size_t segmentIndex) const override { return common.segmentT(segmentIndex); }
    floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const override { return common.segmentLength(segmentIndex, a, b); }

//data
private:
    QuinticHermiteSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingQuinticHermiteSpline<InterpolationType,floating_t>::LoopingQuinticHermiteSpline(
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
    int numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

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
LoopingQuinticHermiteSpline<InterpolationType,floating_t>::LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 3);

    int size = points.size();
    int numSegments = size;

    //compute the T values for each point
    int padding = 2;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //compute the tangents
    std::vector<InterpolationType> tangents(size + 1);
    for(int i = 0; i < size + 1; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = points.at((i - 1 + size)%size);
        InterpolationType pCurrent = points.at((i + size)%size);
        InterpolationType pNext = points.at((i + 1)%size);

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
    std::vector<InterpolationType> curves(size + 1);
    for(int i = 0; i < size + 1; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = tangents.at((i - 1 + size)%size);
        InterpolationType pCurrent = tangents.at(i);
        InterpolationType pNext = tangents.at((i + 1)%size);

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
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename QuinticHermiteSplineCommon<InterpolationType, floating_t>::QuinticHermiteSplinePoint> positionData(numSegments + 1);
    for(int i = 0; i < numSegments + 1; i++)
    {
        knots[i] = indexToT[i];
        positionData[i].position = points[i%size];
        positionData[i].tangent = tangents[i];
        positionData[i].curvature = curves[i];
    }

    common = QuinticHermiteSplineCommon<InterpolationType, floating_t>(std::move(positionData), std::move(knots));
}


template<class InterpolationType, typename floating_t>
InterpolationType LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getPosition(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getTangent(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getCurvature(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getWiggle(wrappedT);
}

template<class InterpolationType, typename floating_t>
floating_t LoopingQuinticHermiteSpline<InterpolationType,floating_t>::arcLength(floating_t a, floating_t b) const
{
    floating_t wrappedA =  SplineSetup::wrapGlobalT(a, maxT);
    floating_t wrappedB =  SplineSetup::wrapGlobalT(b, maxT);

    return common.getLength(wrappedA, wrappedB);
}
