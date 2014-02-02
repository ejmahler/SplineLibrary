#include "looping_natural_spline.h"

#include <cassert>
#include <cmath>

#include "../utils/linearsolver.h"
#include "../utils/t_calculator.h"

LoopingNaturalSpline::LoopingNaturalSpline(const std::vector<Vector3D> &points, double alpha)
{
    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;
    std::unordered_map<int, Vector3D> pointMap;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    indexToT = TCalculator::computeLoopingTValues(points, alpha, 1);
    maxT = indexToT.at(size);

    //now that we know the t values, we need to prepare the tridiagonal matrix calculation
    //note that there several ways to formulate this matrix - i chose the following:
    // http://www-hagen.informatik.uni-kl.de/~alggeom/pdf/ws1213/alggeom_script_ws12_02.pdf

    //the tridiagonal matrix's main diagonal will be neighborDeltaT, and the secondary diagonals will be deltaT
    //the list of values to solve for will be neighborDeltaPoint

    //create an array of the differences in T between one point and the next
    std::vector<double> upperDiagonal;
    for(int i = 0; i < points.size(); i++)
    {
        double delta = indexToT.at(i + 1) - indexToT.at(i);
        upperDiagonal.push_back(delta);
    }

    //create an array that stores 2 * (deltaT.at(i - 1) + deltaT.at(i))
    //when i = 0, wrap i - 1 back around to the end of the list
    std::vector<double> diagonal;
    for(int i = 0; i < points.size(); i++)
    {
        double neighborDelta = 2 * (upperDiagonal.at((i - 1 + size)%size) + upperDiagonal.at(i));
        diagonal.push_back(neighborDelta);
    }

    //create an array of displacement between each point, divided by delta t
    std::vector<Vector3D> deltaPoint;
    for(int i = 0; i < points.size(); i++)
    {
        Vector3D displacement = points.at((i + 1)%size) - points.at(i);
        deltaPoint.push_back(displacement / upperDiagonal.at(i));
    }

    //create an array that stores 3 * (deltaPoint(i - 1) + deltaPoint(i))
    //when i = 0, wrap i - 1 back around to the end of the list
    std::vector<Vector3D> inputVector;
    for(int i = 0; i < points.size(); i++)
    {
        Vector3D neighborDelta = 3 * (deltaPoint.at(i) - deltaPoint.at((i - 1 + size) % size));
        inputVector.push_back(neighborDelta);
    }

    //solve the cyclic tridiagonal system to get the curvature at each point
    std::vector<Vector3D> curvatures = LinearSolver::solveCyclicSymmetricTridiagonal(diagonal, upperDiagonal, inputVector);

    //we now have the curvature for every point
    //use this curvature to determine a,b,c,and d to build each segment
    for(int i = 0; i < numSegments; i++) {

        InterpolationData segment;
        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        double currentDeltaT = segment.t1 - segment.t0;
        Vector3D currentPoint = points.at(i);
        Vector3D nextPoint = points.at((i + 1)%size);
        Vector3D currentCurvature = curvatures.at(i);
        Vector3D nextCurvature = curvatures.at((i + 1)%size);

        segment.a = currentPoint;
        segment.b = (nextPoint - currentPoint) / currentDeltaT - (currentDeltaT / 3) * (nextCurvature + 2*currentCurvature);
        segment.c = currentCurvature;
        segment.d = (nextCurvature - currentCurvature) / (3 * currentDeltaT);

        segment.tDistanceInverse = 1 / currentDeltaT;

        segmentData.push_back(segment);
    }
}


Vector3D LoopingNaturalSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0);

    return computePosition(t, segment);
}

Spline::InterpolatedPT LoopingNaturalSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0);

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC LoopingNaturalSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0);

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW LoopingNaturalSpline::getWiggle(double x) const
{
    int index = getSegmentIndex(x);
    InterpolationData segment = segmentData.at(index);
    double t = (x - segment.t0);

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

bool LoopingNaturalSpline::isLooping(void) const
{
    return true;
}