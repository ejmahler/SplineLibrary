#ifndef LOOPINGNATURALSPLINE_H
#define LOOPINGNATURALSPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"
#include "spline_library/natural/natural_spline_kernel.h"

#include "spline_library/utils/linearalgebra.h"
#include "spline_library/utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class LoopingNaturalSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    LoopingNaturalSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override;
    floating_t getMaxT(void) const override;

    bool isLooping(void) const override;

//data
private:
    //a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
    std::vector<NaturalSplineKernel::InterpolationData<InterpolationType, floating_t>> segmentData;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingNaturalSpline<InterpolationType,floating_t>::LoopingNaturalSpline(const std::vector<InterpolationType> &points, floating_t alpha)
    :Spline<InterpolationType,floating_t>(points)
{

    std::unordered_map<int, floating_t> indexToT_Raw;
    std::unordered_map<int, InterpolationType> pointMap;

    size_t size = points.size();
    int numSegments = size;

    //compute the T values for each point
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, 1);
    maxT = indexToT.at(size);

    //now that we know the t values, we need to prepare the tridiagonal matrix calculation
    //note that there several ways to formulate this matrix - i chose the following:
    // http://www-hagen.informatik.uni-kl.de/~alggeom/pdf/ws1213/alggeom_script_ws12_02.pdf

    //the tridiagonal matrix's main diagonal will be neighborDeltaT, and the secondary diagonals will be deltaT
    //the list of values to solve for will be neighborDeltaPoint

    //create an array of the differences in T between one point and the next
    std::vector<floating_t> upperDiagonal(size);
    for(size_t i = 0; i < size; i++)
    {
        floating_t delta = indexToT.at(i + 1) - indexToT.at(i);
        upperDiagonal[i] = delta;
    }

    //create an array that stores 2 * (deltaT.at(i - 1) + deltaT.at(i))
    //when i = 0, wrap i - 1 back around to the end of the list
    std::vector<floating_t> diagonal(size);
    for(size_t i = 0; i < size; i++)
    {
        floating_t neighborDelta = 2 * (upperDiagonal.at((i - 1 + size)%size) + upperDiagonal.at(i));
        diagonal[i] = neighborDelta;
    }

    //create an array of displacement between each point, divided by delta t
    std::vector<InterpolationType> deltaPoint(size);
    for(size_t i = 0; i < size; i++)
    {
        InterpolationType displacement = points.at((i + 1)%size) - points.at(i);
        deltaPoint[i] = displacement / upperDiagonal.at(i);
    }

    //create an array that stores 3 * (deltaPoint(i - 1) + deltaPoint(i))
    //when i = 0, wrap i - 1 back around to the end of the list
    std::vector<InterpolationType> inputVector(size);
    for(size_t i = 0; i < size; i++)
    {
        InterpolationType neighborDelta = 3 * (deltaPoint.at(i) - deltaPoint.at((i - 1 + size) % size));
        inputVector[i] = neighborDelta;
    }

    //solve the cyclic tridiagonal system to get the curvature at each point
    std::vector<InterpolationType> curvatures = LinearAlgebra::solveCyclicSymmetricTridiagonal(
                std::move(diagonal),
                std::move(upperDiagonal),
                std::move(inputVector)
                );

    //we now have the curvature for every point
    //use this curvature to determine a,b,c,and d to build each segment
    for(int i = 0; i < numSegments; i++) {

        NaturalSplineKernel::InterpolationData<InterpolationType, floating_t> segment;
        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        floating_t currentDeltaT = segment.t1 - segment.t0;
        InterpolationType currentPoint = points.at(i);
        InterpolationType nextPoint = points.at((i + 1)%size);
        InterpolationType currentCurvature = curvatures.at(i);
        InterpolationType nextCurvature = curvatures.at((i + 1)%size);

        segment.a = currentPoint;
        segment.b = (nextPoint - currentPoint) / currentDeltaT - (currentDeltaT / 3) * (nextCurvature + 2*currentCurvature);
        segment.c = currentCurvature;
        segment.d = (nextCurvature - currentCurvature) / (3 * currentDeltaT);

        segment.tDistanceInverse = 1 / currentDeltaT;

        segmentData.push_back(segment);
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingNaturalSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, segmentData.size());
    if(globalT < 0)
        globalT += segmentData.size();

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingNaturalSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, segmentData.size());
    if(globalT < 0)
        globalT += segmentData.size();

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingNaturalSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, segmentData.size());
    if(globalT < 0)
        globalT += segmentData.size();

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingNaturalSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, segmentData.size());
    if(globalT < 0)
        globalT += segmentData.size();

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT),
                segment.computeWiggle()
                );
}

template<class InterpolationType, typename floating_t>
floating_t LoopingNaturalSpline<InterpolationType,floating_t>::getT(int index) const
{
    return indexToT.at(index);
}

template<class InterpolationType, typename floating_t>
floating_t LoopingNaturalSpline<InterpolationType,floating_t>::getMaxT(void) const
{
    return maxT;
}

template<class InterpolationType, typename floating_t>
bool LoopingNaturalSpline<InterpolationType,floating_t>::isLooping(void) const
{
    return true;
}

#endif // LOOPINGNATURALSPLINE_H
