#ifndef NATURALSPLINE_H
#define NATURALSPLINE_H

#include "spline_library/spline.h"
#include "spline_library/natural/natural_spline_kernel.h"

#include "spline_library/utils/linearalgebra.h"
#include "spline_library/utils/spline_setup.h"

#include <unordered_map>
#include <cassert>

template<class InterpolationType, typename floating_t=float>
class NaturalSpline final : public Spline<InterpolationType, floating_t>
{
public:
    enum EndConditions { Natural, NotAKnot };

//constructors
public:
    NaturalSpline(const std::vector<InterpolationType> &points,
                  bool includeEndpoints = true,
                  floating_t alpha = 0.0,
                  EndConditions endConditions = Natural);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override;
    floating_t getMaxT(void) const override;

    bool isLooping(void) const override;

//methods
private:
    std::vector<InterpolationType> computeCurvaturesNatural(void) const;
    std::vector<InterpolationType> computeCurvaturesNotAKnot(void) const;

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
NaturalSpline<InterpolationType,floating_t>::NaturalSpline(const std::vector<InterpolationType> &points,
                                                           bool includeEndpoints,
                                                           floating_t alpha,
                                                           EndConditions endConditions)
    :Spline<InterpolationType,floating_t>(points)
{
    size_t size = points.size();
    int firstPoint;

    int numSegments;
    if(includeEndpoints)
    {
        assert(points.size() >= 3);
        numSegments = size - 1;
        firstPoint = 0;
    }
    else
    {
        assert(points.size() >= 4);
        numSegments = size - 3;
        firstPoint = 1;
    }

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, firstPoint);
    maxT = indexToT.at(firstPoint + numSegments);

    //next we compute curvatures
    std::vector<InterpolationType> curvatures;
    if(endConditions == Natural)
        curvatures = computeCurvaturesNatural();
    else
        curvatures = computeCurvaturesNotAKnot();

    //we now have 0 curvature for index 0 and n - 1, and the final (usually nonzero) curvature for every other point
    //use this curvature to determine a,b,c,and d to build each segment
    for(int i = firstPoint; i < firstPoint + numSegments; i++) {

        NaturalSplineKernel::InterpolationData<InterpolationType, floating_t> segment;
        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        floating_t currentDeltaT = segment.t1 - segment.t0;
        InterpolationType currentPoint = points.at(i);
        InterpolationType nextPoint = points.at(i + 1);
        InterpolationType currentCurvature = curvatures.at(i);
        InterpolationType nextCurvature = curvatures.at(i + 1);

        segment.a = currentPoint;
        segment.b = (nextPoint - currentPoint) / currentDeltaT - (currentDeltaT / 3) * (nextCurvature + 2*currentCurvature);
        segment.c = currentCurvature;
        segment.d = (nextCurvature - currentCurvature) / (3 * currentDeltaT);

        segment.tDistanceInverse = 1 / currentDeltaT;

        segmentData.push_back(segment);
    }
}

template<class InterpolationType, typename floating_t>
std::vector<InterpolationType> NaturalSpline<InterpolationType,floating_t>::computeCurvaturesNatural(void) const
{

    //now that we know the t values, we need to prepare the tridiagonal matrix calculation
    //note that there several ways to formulate this matrix - for the "natural boundary conditions" i chose the following:
    // http://www-hagen.informatik.uni-kl.de/~alggeom/pdf/ws1213/alggeom_script_ws12_02.pdf

    //the tridiagonal matrix's main diagonal will be neighborDeltaT, and the secondary diagonals will be deltaT
    //the list of values to solve for will be neighborDeltaPoint

    size_t loop_limit = this->getOriginalPoints().size() - 1;

    //create an array of the differences in T between one point and the next
    std::vector<floating_t> upperDiagonal(loop_limit);
    for(size_t i = 0; i < loop_limit; i++)
    {
        floating_t delta = indexToT.at(i + 1) - indexToT.at(i);
        upperDiagonal[i] = delta;
    }

    //create an array that stores 2 * (deltaT.at(i - 1) + deltaT.at(i))
    std::vector<floating_t> diagonal(loop_limit - 1);
    for(size_t i = 1; i < loop_limit; i++)
    {
        floating_t neighborDelta = 2 * (upperDiagonal.at(i - 1) + upperDiagonal.at(i));
        diagonal[i - 1] = neighborDelta;
    }

    //create an array of displacement between each point, divided by delta t
    std::vector<InterpolationType> deltaPoint(loop_limit);
    for(size_t i = 0; i < loop_limit; i++)
    {
        InterpolationType displacement = this->getOriginalPoints().at(i + 1) - this->getOriginalPoints().at(i);
        deltaPoint[i] = displacement / upperDiagonal.at(i);
    }

    //create an array that stores 3 * (deltaPoint(i - 1) + deltaPoint(i))
    std::vector<InterpolationType> inputVector(loop_limit - 1);
    for(size_t i = 1; i < loop_limit; i++)
    {
        InterpolationType neighborDelta = 3 * (deltaPoint.at(i) - deltaPoint.at(i - 1));
        inputVector[i - 1] = neighborDelta;
    }

    //the first element in upperDiagonal is garbage, so remove it
    upperDiagonal.erase(upperDiagonal.begin());

    //solve the tridiagonal system to get the curvature at each point
    std::vector<InterpolationType> curvatures = LinearAlgebra::solveSymmetricTridiagonal(
                std::move(diagonal),
                std::move(upperDiagonal),
                std::move(inputVector)
                );

    //we didn't compute the first or last curvature, which will be 0
    curvatures.insert(curvatures.begin(), InterpolationType());
    curvatures.push_back(InterpolationType());

    return curvatures;
}


template<class InterpolationType, typename floating_t>
std::vector<InterpolationType> NaturalSpline<InterpolationType,floating_t>::computeCurvaturesNotAKnot(void) const
{

    //now that we know the t values, we need to prepare the tridiagonal matrix calculation
    //note that there several ways to formulate this matrix; for "not a knot" i chose the following:
    // http://sepwww.stanford.edu/data/media/public/sep//sergey/128A/answers6.pdf

    //the tridiagonal matrix's main diagonal will be neighborDeltaT, and the secondary diagonals will be deltaT
    //the list of values to solve for will be neighborDeltaPoint

    size_t size = this->getOriginalPoints().size() - 1;

    //create an array of the differences in T between one point and the next
    std::vector<floating_t> deltaT(size);
    for(size_t i = 0; i < size; i++)
    {
        deltaT[i] = indexToT.at(i + 1) - indexToT.at(i);
    }

    //the main diagonal of the tridiagonal will be 2 * (deltaT[i] + deltaT[i + 1])
    float mainDiagonalSize = size - 1;
    std::vector<floating_t> mainDiagonal(mainDiagonalSize);
    for(size_t i = 0; i < mainDiagonalSize; i++)
    {
        mainDiagonal[i] = 2 * (deltaT[i] + deltaT[i + 1]);
    }

    //the upper diagonal will just be deltaT[i + 1]
    float secondaryDiagonalSize = size - 2;
    std::vector<floating_t> upperDiagonal(secondaryDiagonalSize);
    for(size_t i = 0; i < secondaryDiagonalSize; i++)
    {
        upperDiagonal[i] = deltaT[i + 1];
    }

    //the lower diagonal is just a copy of the upper diagonal
    std::vector<floating_t> lowerDiagonal = upperDiagonal;

    //create an array of displacement between each point, divided by delta t
    std::vector<InterpolationType> deltaPoint(size);
    for(size_t i = 0; i < size; i++)
    {
        InterpolationType displacement = this->getOriginalPoints().at(i + 1) - this->getOriginalPoints().at(i);
        deltaPoint[i] = displacement / deltaT.at(i);
    }

    //create an array that stores 3 * (deltaPoint(i - 1) + deltaPoint(i))
    std::vector<InterpolationType> inputVector(mainDiagonalSize);
    for(size_t i = 0; i < mainDiagonalSize; i++)
    {
        inputVector[i] = 3 * (deltaPoint.at(i + 1) - deltaPoint.at(i));
    }

    //the first and last of the values in maindiagonalare different than normal
    mainDiagonal[0] = 3*deltaT[0] + 2*deltaT[1] + deltaT[0]*deltaT[0]/deltaT[1];
    mainDiagonal[mainDiagonalSize - 1] = 3*deltaT[size - 1] + 2*deltaT[size - 2] + deltaT[size - 1]*deltaT[size - 1]/deltaT[size - 2];

    //the first value in the upper diagonal is different than normal
    upperDiagonal[0] = deltaT[1] - deltaT[0]*deltaT[0]/deltaT[1];

    //the last value in the upper diagonal is different than normal
    lowerDiagonal[secondaryDiagonalSize - 1] = deltaT[size - 2] - deltaT[size - 1]*deltaT[size - 1]/deltaT[size - 2];

    //solve the tridiagonal system to get the curvature at each point
    std::vector<InterpolationType> curvatures = LinearAlgebra::solveTridiagonal(
                std::move(mainDiagonal),
                std::move(upperDiagonal),
                std::move(lowerDiagonal),
                std::move(inputVector)
                );

    //we didn't compute the first or last curvature, which will be calculated based on the others
    curvatures.insert(curvatures.begin(), curvatures[0] * (1 + deltaT[0]/deltaT[1]) - curvatures[1] * (deltaT[0]/deltaT[1]));
    curvatures.push_back(curvatures.back() * (1 + deltaT[size - 1]/deltaT[size - 2])
            - curvatures[curvatures.size() - 2] * (deltaT[size - 1]/deltaT[size - 2]));

    return curvatures;
}

template<class InterpolationType, typename floating_t>
InterpolationType NaturalSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    NaturalSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    NaturalSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
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
    NaturalSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
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
floating_t NaturalSpline<InterpolationType,floating_t>::getT(int index) const
{
    return indexToT.at(index);
}

template<class InterpolationType, typename floating_t>
floating_t NaturalSpline<InterpolationType,floating_t>::getMaxT(void) const
{
    return maxT;
}

template<class InterpolationType, typename floating_t>
bool NaturalSpline<InterpolationType,floating_t>::isLooping(void) const
{
    return false;
}

#endif // NATURALSPLINE_H
