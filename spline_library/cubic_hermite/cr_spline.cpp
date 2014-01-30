
#include "cr_spline.h"
#include "../utils/t_calculator.h"

#include <cmath>
#include <cassert>

CRSpline::CRSpline(const std::vector<Vector3D> &points, double alpha)
{
    assert(points.size() >= 4);

    this->points = points;

    int size = points.size();
    int firstTangent = 1;
    numSegments = size - 3;

    //compute the T values for each point
    indexToT = TCalculator::computeTValues(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = firstTangent; i < firstTangent + numSegments + 1; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = points.at(i - 1);
        Vector3D pCurrent = points.at(i);
        Vector3D pNext = points.at(i + 1);

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
    for(int i = firstTangent; i < firstTangent + numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}
