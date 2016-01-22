#ifndef NATURALSPLINE_H
#define NATURALSPLINE_H

#include "spline_library/spline.h"
#include "spline_library/natural/natural_spline_common.h"

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
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t arcLength(floating_t a, floating_t b) const override { if(a > b) std::swap(a, b); return common.getLength(a, b); }
    floating_t totalLength(void) const override { return common.getTotalLength(); }

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

//methods
private:
    std::vector<InterpolationType> computeCurvaturesNatural(void) const;
    std::vector<InterpolationType> computeCurvaturesNotAKnot(void) const;

//data
private:
    NaturalSplineCommon<InterpolationType, floating_t> common;

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
    std::vector<floating_t> knots(numSegments + 1);
    std::vector<typename NaturalSplineCommon<InterpolationType, floating_t>::NaturalSplineSegment> segments(numSegments + 1);
    for(int i = firstPoint; i < numSegments + firstPoint + 1; i++) {

        knots[i - firstPoint] = indexToT.at(i);
        segments[i - firstPoint].a = points.at(i);
        segments[i - firstPoint].c = curvatures.at(i);
    }

    common = NaturalSplineCommon<InterpolationType, floating_t>(std::move(segments), std::move(knots));
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

#endif // NATURALSPLINE_H
