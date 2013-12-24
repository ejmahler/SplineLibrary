#include "quintic_hermite_spline.h"

/*
Vector3D QuinticHermiteSpline::computePosition2(double t, const InterpolationData &segment) const
{
    double a = segment.x0;
    double b = segment.x1;
    double c = segment.x2;
    double d = segment.x3;
    double e = segment.x4;
    double f = segment.x5;

    Vector3D W01 = segment.P0 * (b - t) / (b - a) + segment.P1 * (t - a) / (b - a);
    Vector3D W12 = segment.P1 * (c - t) / (c - b) + segment.P2 * (t - b) / (c - b);
    Vector3D W23 = segment.P2 * (d - t) / (d - c) + segment.P3 * (t - c) / (d - c);
    Vector3D W34 = segment.P3 * (e - t) / (e - d) + segment.P4 * (t - d) / (e - d);
    Vector3D W45 = segment.P4 * (f - t) / (f - e) + segment.P5 * (t - e) / (f - e);

    Vector3D X1 = W01 * (c - t) / (c - a) + W12 * (t - a) / (c - a);
    Vector3D X2 = W12 * (d - t) / (d - b) + W23 * (t - b) / (d - b);
    Vector3D X3 = W23 * (e - t) / (e - c) + W34 * (t - c) / (e - c);
    Vector3D X4 = W34 * (f - t) / (f - d) + W45 * (t - d) / (f - d);

    Vector3D Y12 = X1 * (c - t) / (c - b) + X2 * (t - b) / (c - b);
    Vector3D Y23 = X2 * (d - t) / (d - c) + X3 * (t - c) / (d - c);
    Vector3D Y34 = X3 * (e - t) / (e - d) + X4 * (t - d) / (e - d);

    Vector3D Z2 = Y12 * (d - t) / (d - b) + Y23 * (t - b) / (d - b);
    Vector3D Z3 = Y23 * (e - t) / (e - c) + Y34 * (t - c) / (e - c);

    return Z2 * (d - t) / (d - c) + Z3 * (t - c) / (d - c);
}*/

int QuinticHermiteSpline::getSegmentIndex(double x) const
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
