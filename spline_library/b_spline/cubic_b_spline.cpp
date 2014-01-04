#include "cubic_b_spline.h"


#include <cmath>
#include <cassert>

CubicBSpline::CubicBSpline()
{

}

CubicBSpline::CubicBSpline(const std::vector<Vector3D> &points)
{
    assert(points.size() >= 4);
    double alpha = 0.0;

    this->points = points;

    std::unordered_map<int, double> indexToT_Raw;
    std::unordered_map<int, Vector3D> pointMap;

    int size = points.size();

    numSegments = size - 3;

    int firstSegment = 1;
    int lastSegment = firstSegment + numSegments;

    //we know points[1] will have a t value of 0
    indexToT_Raw[1] = 0;
    pointMap[1] = points[1];

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = 0; i >= 0; i--)
    {
        //points[1] is a control point, so give it a negative t value, so that the first actual point can have a t value of 0
        double distance = (points.at(i) - points.at(i + 1)).length();
        indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
        pointMap[i] = points[i];
    }

    //compute the t values of the other points
    for(int i = 2; i < size; i++)
    {
        double distance = (points.at(i) - points.at(i - 1)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

        pointMap[i] = points.at(i);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(lastSegment);

    //now that we have all of our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(lastSegment);


    //pre-arrange the data needed for interpolation
    for(int i = firstSegment; i < lastSegment; i++)
    {
        InterpolationData segment;

        segment.t1 = indexToT.at(i);
        segment.t2 = indexToT.at(i + 1);

        segment.p0 = pointMap.at(i - 1);
        segment.p1 = pointMap.at(i);
        segment.p2 = pointMap.at(i + 1);
        segment.p3 = pointMap.at(i + 2);

        segmentData.push_back(segment);
    }
}

Vector3D CubicBSpline::getPosition(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return computePosition(t, segment);
}

Spline::InterpolatedPT CubicBSpline::getTangent(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC CubicBSpline::getCurvature(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW CubicBSpline::getWiggle(double x) const
{
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = x - segment.t1;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

int CubicBSpline::getSegmentIndex(double x) const
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(x < segmentData[0].t1)
        return 0;
    if(x > segmentData[numSegments - 1].t2)
        return numSegments - 1;

    //perform a binary search on segmentData
    int currentMin = 0;
    int currentMax = segmentData.size() - 1;
    int currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x

    while(x < segmentData[currentIndex].t1 || x > segmentData[currentIndex].t2)
    {
        //if t0 is greater than x, search the left half of the array
        if(segmentData[currentIndex].t1 > x)
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

double CubicBSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicBSpline::getMaxT(void) const
{
    return maxT;
}

int CubicBSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &CubicBSpline::getPoints(void) const
{
    return points;
}

bool CubicBSpline::isLooping(void) const
{
    return false;
}
