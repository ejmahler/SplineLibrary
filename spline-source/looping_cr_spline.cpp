#include "looping_cr_spline.h"


#include <cmath>
#include <cassert>

LoopingCRSpline::LoopingCRSpline(const std::vector<Vector3D> &points, double alpha)
    :points(points)
{
    assert(points.size() >= 6);

    std::unordered_map<int, double> indexToT_Raw;
    std::unordered_map<int, Vector3D> pointMap;

    int size = points.size();
    numSegments = size;

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;
    pointMap[0] = points[0];

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -1; i--)
    {
        //points[1] is a control point, so give it a nagative t value, so that the first actual point can have a t value of 0
        double distance = (points.at((i + size)%size) - points.at((i + 1 + size)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
        pointMap[i] = points.at((i + size)%size);
    }

    //compute the t values of the other points
    for(int i = 1; i < size + 2; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);

        pointMap[i] = points.at((i + size)%size);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = numSegments * it->second / maxTRaw;
    }
    maxT = indexToT.at(size);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = 0; i < size + 1; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = pointMap.at(i - 1);
        Vector3D pCurrent = pointMap.at(i);
        Vector3D pNext = pointMap.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        tangentMap[i] =
                            pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                            + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

                         //plus a little something extra - this is derived from the pyramid contruction
                         //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
                         //yielding the standard catmull-rom formula
                            - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = pointMap.at(i);
        segment.p1 = pointMap.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

LoopingCRSpline::~LoopingCRSpline()
{

}

Vector3D LoopingCRSpline::getPosition(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

InterpolatedPT LoopingCRSpline::getTangent(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
        computePosition(t, segment),
        computeTangent(t, segment)
        );
}

InterpolatedPTC LoopingCRSpline::getCurvature(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    //find the interpolation data for this t value
    InterpolationData segment = segmentData.at(getSegmentIndex(x));
    double t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
        computePosition(t, segment),
        computeTangent(t, segment),
        computeCurvature(t, segment)
        );
}

double LoopingCRSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LoopingCRSpline::getMaxT(void) const
{
    return maxT;
}

int LoopingCRSpline::getNumSegments(void) const
{
    return numSegments;
}

const std::vector<Vector3D> &LoopingCRSpline::getPoints(void) const
{
    return points;
}


bool LoopingCRSpline::isLoop(void) const
{
    return true;
}
