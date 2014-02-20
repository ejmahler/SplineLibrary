#include "splinelengthcalculator.h"

#include "spline.h"

#include <cmath>

#define LENGTH_INTERVAL 0.1
#define MINIMUM_SEGMENTS 5
#define RECURSIVE_MINIMUM_INTERVAL .00001

SplineLengthCalculator::SplineLengthCalculator(const std::shared_ptr<Spline> &spline)
    :
      spline(spline),
      maxT(spline->getMaxT()),

      //the spline length will be lazy-computed when it's needed
      splineLength(-1)
{

}

double SplineLengthCalculator::findLength(double beginT, double endT, bool useShortestPath) const
{
    //if this is a looping spline and the calles had requested the shortest path, the behavior will be slightly different
    if(useShortestPath && spline->isLooping())
    {
        //use fmod on both args to make sure we're not going around the loop multiple times
        beginT = fmod(beginT, maxT);
        endT = fmod(endT, maxT);

        //make sure they're both positive
        beginT = (beginT < 0) ? (beginT + maxT) : beginT;
        endT = (endT < 0) ? (endT + maxT) : endT;

        //make sure beginT is less than endT
        double actualBeginT = std::min(beginT, endT);
        double actualEndT = std::max(beginT,endT);

        //compute length
        double computedLength = computeLength(actualBeginT, actualEndT);

        //lazy compute the spline length
        if(splineLength < 0)
        {
            //modifying a variable inside of a const function can be non thread safe, but in this case thread safety is not a concern
            //we're using std::atomic so the worst possible thing that could happen is this variable could be computed multiple times by simultaneous threads
            //but if that happens, it'll be the same both times so who cares
            splineLength = findLength(0, maxT, false);
        }

        //that was one direction around, the other direction will be (splineLength - computedLength). return the smaller of the two
        return std::min(computedLength, splineLength - computedLength);
    }
    else
    {
        //make sure beginT is less than endT
        double actualBeginT = std::min(beginT, endT);
        double actualEndT = std::max(beginT,endT);

        return computeLength(actualBeginT, actualEndT);
    }
}

double SplineLengthCalculator::computeLength(double beginT, double endT) const
{
    double tDistance = endT - beginT;

    //find the number of segments we're going to use
    double interval = LENGTH_INTERVAL;
    double numSegments = tDistance / interval;
    if(tDistance / interval < MINIMUM_SEGMENTS)
    {
        numSegments = MINIMUM_SEGMENTS;
    }
    else
    {
        //the nuber of segments probably isn't a whole number, so round it up
        numSegments = ceil(numSegments);
    }

    //now that we know the number of segments, find the actual interval
    interval = tDistance / numSegments;

    //for each segment, compute the length of a circle arc passing though that segment's endpoints
    auto previousData = spline->getTangent(beginT);
    previousData.tangent = previousData.tangent.normalized();
    double currentT = beginT + interval;

    double lengthSum = 0;

    while(currentT <= endT)
    {
        auto currentData = spline->getTangent(currentT);
        currentData.tangent = currentData.tangent.normalized();

        lengthSum += computeLengthHelper(
                    currentT - interval,    previousData.position, previousData.tangent,
                    currentT,               currentData.position, currentData.tangent);

        previousData = currentData;
        currentT += interval;
    }

    return lengthSum;
}

double SplineLengthCalculator::computeLengthHelper(
        double beginT, const Vector3D &beginPosition, const Vector3D &beginTangentNormalized,
        double endT, const Vector3D &endPosition, const Vector3D &endTangentNormalized) const
{
    //compute the angle between the tangents
    double cosTangentAngle = Vector3D::dotProduct(beginTangentNormalized, endTangentNormalized);

    //if the cos(angle) between the velocities is almost 1, there is essentially a straight line between the two points
    //if that's the case, just return the distance between the two points
    if(cosTangentAngle > .9999 || (endT - beginT < RECURSIVE_MINIMUM_INTERVAL))
    {
        return (endPosition - beginPosition).length();
    }
    else
    {
        double middleT = (beginT + endT) * 0.5;
        auto middleData = spline->getTangent(middleT);

        Vector3D middleTangentNormalized = middleData.tangent.normalized();

        double firstHalf = computeLengthHelper(beginT, beginPosition, beginTangentNormalized,
                                                      middleT, middleData.position, middleTangentNormalized);

        double secondHalf = computeLengthHelper(middleT, middleData.position, middleTangentNormalized,
                                                      endT, endPosition, endTangentNormalized);

        return firstHalf + secondHalf;
    }
}
