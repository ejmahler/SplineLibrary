#include "cubic_hermite_spline.h"

#include <cmath>
#include <cassert>

CubicHermiteSpline::CubicHermiteSpline()
{

}

CubicHermiteSpline::CubicHermiteSpline(const std::vector<Vector3D> &points, const std::vector<Vector3D> &tangents, double alpha)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;

    int size = points.size();

    numSegments = size - 1;

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;

    //compute the t values of the other points
    for(int i = 1; i < size; i++)
    {
        double distance = (points.at(i) - points.at(i - 1)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(numSegments);

    //now that we have all of our t values and indexes figured out, normalize the t values by dividing them by maxT and multplying by numSegments
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D CubicHermiteSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

InterpolatedPT CubicHermiteSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

InterpolatedPTC CubicHermiteSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

int CubicHermiteSpline::getSegmentIndex(double x) const
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

double CubicHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicHermiteSpline::getMaxT(void) const
{
    return maxT;
}

int CubicHermiteSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &CubicHermiteSpline::getPoints(void) const
{
    return points;
}

bool CubicHermiteSpline::isLooping(void) const
{
    return false;
}
