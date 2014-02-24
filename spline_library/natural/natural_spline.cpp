#include "natural_spline.h"

#include <cassert>
#include <cmath>

#include "spline_library/utils/linearalgebra.h"
#include "spline_library/utils/t_calculator.h"

NaturalSpline::NaturalSpline()
{
}

NaturalSpline::NaturalSpline(const std::vector<Vector3D> &points, bool includeEndpoints, double alpha)
{
    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;

    size_t size = points.size();
    int firstPoint;

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
    indexToT = TCalculator::computeTValues(points, alpha, firstPoint);
    maxT = indexToT.at(firstPoint + numSegments);

    //now that we know the t values, we need to prepare the tridiagonal matrix calculation
    //note that there several ways to formulate this matrix - i chose the following:
    // http://www-hagen.informatik.uni-kl.de/~alggeom/pdf/ws1213/alggeom_script_ws12_02.pdf

    //the tridiagonal matrix's main diagonal will be neighborDeltaT, and the secondary diagonals will be deltaT
    //the list of values to solve for will be neighborDeltaPoint

    size_t loop_limit = size - 1;

    //create an array of the differences in T between one point and the next
    std::vector<double> upperDiagonal(loop_limit);
    for(size_t i = 0; i < loop_limit; i++)
    {
        double delta = indexToT.at(i + 1) - indexToT.at(i);
        upperDiagonal[i] = delta;
    }

    //create an array that stores 2 * (deltaT.at(i - 1) + deltaT.at(i))
    std::vector<double> diagonal(loop_limit - 1);
    for(size_t i = 1; i < loop_limit; i++)
    {
        double neighborDelta = 2 * (upperDiagonal.at(i - 1) + upperDiagonal.at(i));
        diagonal[i - 1] = neighborDelta;
    }

    //create an array of displacement between each point, divided by delta t
    std::vector<Vector3D> deltaPoint(loop_limit);
    for(size_t i = 0; i < loop_limit; i++)
    {
        Vector3D displacement = points.at(i + 1) - points.at(i);
        deltaPoint[i] = displacement / upperDiagonal.at(i);
    }

    //create an array that stores 3 * (deltaPoint(i - 1) + deltaPoint(i))
    std::vector<Vector3D> inputVector(loop_limit - 1);
    for(size_t i = 1; i < loop_limit; i++)
    {
        Vector3D neighborDelta = 3 * (deltaPoint.at(i) - deltaPoint.at(i - 1));
        inputVector[i - 1] = neighborDelta;
    }

    //the first element in upperDiagonal is garbage, so remove it
    upperDiagonal.erase(upperDiagonal.begin());

    //solve the tridiagonal system to get the curvature at each point
    std::vector<Vector3D> curvatures = LinearAlgebra::solveSymmetricTridiagonal(diagonal, upperDiagonal, inputVector);

    //we didn't compute the first or last curvature, which will be 0
    curvatures.insert(curvatures.begin(), Vector3D());
    curvatures.push_back(Vector3D());

    //we now have 0 curvature for index 0 and n - 1, and the final (usually nonzero) curvature for every other point
    //use this curvature to determine a,b,c,and d to build each segment
    for(int i = firstPoint; i < firstPoint + numSegments; i++) {

        InterpolationData segment;
        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        double currentDeltaT = segment.t1 - segment.t0;
        Vector3D currentPoint = points.at(i);
        Vector3D nextPoint = points.at(i + 1);
        Vector3D currentCurvature = curvatures.at(i);
        Vector3D nextCurvature = curvatures.at(i + 1);

        segment.a = currentPoint;
        segment.b = (nextPoint - currentPoint) / currentDeltaT - (currentDeltaT / 3) * (nextCurvature + 2*currentCurvature);
        segment.c = currentCurvature;
        segment.d = (nextCurvature - currentCurvature) / (3 * currentDeltaT);

        segment.tDistanceInverse = 1 / currentDeltaT;

        segmentData.push_back(segment);
    }
}


Vector3D NaturalSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0);

    return computePosition(t, segment);
}

Spline::InterpolatedPT NaturalSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0);

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC NaturalSpline::getCurvature(double x) const
{
    int index = getSegmentIndex(x);
    InterpolationData segment = segmentData.at(index);
    double t = (x - segment.t0);

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW NaturalSpline::getWiggle(double x) const
{
    int index = getSegmentIndex(x);
    InterpolationData segment = segmentData.at(index);
    double t = (x - segment.t0);

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(segment)
                );
}

int NaturalSpline::getSegmentIndex(double x) const
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(x < segmentData[0].t0)
        return 0;
    if(x > segmentData[numSegments - 1].t1)
        return numSegments - 1;

    //perform a binary search on segmentData
    int currentMin = 0;
    int currentMax = segmentData.size() - 1;
    int currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x

    while(x < segmentData[currentIndex].t0 || x > segmentData[currentIndex].t1)
    {
        //if t0 is greater than x, search the left half of the array
        if(segmentData[currentIndex].t0 > x)
        {
            currentMax = currentIndex - 1;
        }

        //the only other possibility is that t1 is less than x, so search the right half of the array
        else
        {
            currentMin = currentIndex + 1;
        }
        currentIndex = (currentMin + currentMax) / 2;
    }
    return currentIndex;
}

double NaturalSpline::getT(int index) const
{
    return indexToT.at(index);
}

double NaturalSpline::getMaxT(void) const
{
    return maxT;
}

const std::vector<Vector3D> &NaturalSpline::getPoints(void) const
{
    return points;
}

bool NaturalSpline::isLooping(void) const
{
    return false;
}
