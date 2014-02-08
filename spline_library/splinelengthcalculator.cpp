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

      //compute the length of the whole spline, to optimize "shortest path" length calculations later
      splineLength(findLengthPrecise(0, maxT, false))
{

}

double SplineLengthCalculator::findLengthFast(double beginT, double endT, bool useShortestPath) const
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
        double computedLength = computeLengthFast(actualBeginT, actualEndT);

        //that was one direction around, the other direction will be (splineLength - computedLength). return the smaller of the two
        return std::min(computedLength, splineLength - computedLength);
    }
    else
    {
        //make sure beginT is less than endT
        double actualBeginT = std::min(beginT, endT);
        double actualEndT = std::max(beginT,endT);

        return computeLengthFast(actualBeginT, actualEndT);
    }
}

double SplineLengthCalculator::findLengthPrecise(double beginT, double endT, bool useShortestPath) const
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
        double computedLength = computeLengthPrecise(actualBeginT, actualEndT);

        //that was one direction around, the other direction will be (splineLength - computedLength). return the smaller of the two
        return std::min(computedLength, splineLength - computedLength);
    }
    else
    {
        //make sure beginT is less than endT
        double actualBeginT = std::min(beginT, endT);
        double actualEndT = std::max(beginT,endT);

        return computeLengthPrecise(actualBeginT, actualEndT);
    }
}




double SplineLengthCalculator::computeLengthFast(double beginT, double endT) const
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
    double currentT = beginT + interval;

    double lengthSum = 0;

    while(currentT <= endT)
    {
        auto currentData = spline->getTangent(currentT);

        lengthSum += computeArcLength(previousData.position, previousData.tangent, currentData.position, currentData.tangent);

        previousData = currentData;
        currentT += interval;
    }

    return lengthSum;
}

double SplineLengthCalculator::computeArcLength(
        const Vector3D &beginPosition, const Vector3D &beginTangent,
        const Vector3D &endPosition, const Vector3D &endTangent) const
{

    Vector3D displacement = endPosition - beginPosition;

    double displacementLength = displacement.length();

    //compute the angle between the tangents
    Vector3D endTangentNormalized = endTangent.normalized();
    Vector3D beginTangentNormalized = beginTangent.normalized();
    double cosTangentAngle = Vector3D::dotProduct(beginTangentNormalized, endTangentNormalized);

    //if the cos(angle) between the velocities is almost 1, there is essentially a straight line between the two points
    //if that's the case, just return the distance between the two points
    if(cosTangentAngle > .99999)
    {
        return displacementLength;
    }
    else
    {
        //construct an isoceles tringle whose base is the sample displacement, and whose opposite angle is the angle between the two tangents
        double arcAngle = acos(cosTangentAngle);
        double radius = displacementLength / (2 * sin(arcAngle * 0.5));

        return arcAngle * radius;
    }
}





double SplineLengthCalculator::computeLengthPrecise(double beginT, double endT) const
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
    double currentT = beginT + interval;

    double lengthSum = 0;

    while(currentT <= endT)
    {
        auto currentData = spline->getTangent(currentT);

        lengthSum += computeLengthPreciseHelper(
                    currentT - interval,    previousData.position, previousData.tangent,
                    currentT,               currentData.position, currentData.tangent);

        previousData = currentData;
        currentT += interval;
    }

    return lengthSum;
}

double SplineLengthCalculator::computeLengthPreciseHelper(
        double beginT, const Vector3D &beginPosition, const Vector3D &beginTangent,
        double endT, const Vector3D &endPosition, const Vector3D &endTangent) const
{
    //compute the angle between the tangents
    Vector3D endTangentNormalized = endTangent.normalized();
    Vector3D beginTangentNormalized = beginTangent.normalized();
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

        double firstHalf = computeLengthPreciseHelper(beginT, beginPosition, beginTangent,
                                                      middleT, middleData.position, middleData.tangent);

        double secondHalf = computeLengthPreciseHelper(middleT, middleData.position, middleData.tangent,
                                                      endT, endPosition, endTangent);

        return firstHalf + secondHalf;
    }
}
