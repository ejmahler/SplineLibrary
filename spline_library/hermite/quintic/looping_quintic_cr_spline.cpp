#include "looping_quintic_cr_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LoopingQuinticCRSpline::LoopingQuinticCRSpline(const std::vector<Vector3D> &points, double alpha)
{
    assert(points.size() >= 3);

    this->points = points;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 2;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = -1; i < size + 2; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = points.at((i - 1 + size)%size);
        Vector3D pCurrent = points.at((i + size)%size);
        Vector3D pNext = points.at((i + 1 + size)%size);

        //the tangent is the standard catmull-rom spline tangent calculation
        tangentMap[i] =
                  pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

             //plus a little something extra - this is derived from the pyramid contruction
             //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
             //yielding the standard catmull-rom formula
                - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }

    //compute the curvatures
    std::map<int, Vector3D> curveMap;
    for(int i = 0; i < size + 1; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = tangentMap.at(i - 1);
        Vector3D pCurrent = tangentMap.at(i);
        Vector3D pNext = tangentMap.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        curveMap[i] =
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

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curveMap.at(i) * tDistance * tDistance;
        segment.c1 = curveMap.at(i + 1) * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}
