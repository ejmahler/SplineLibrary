#include "splinelengthcalculator.h"

#include "spline.h"

#include <cmath>

#define LENGTH_INTERVAL 0.25
#define RECURSIVE_MINIMUM_INTERVAL .00001

SplineLengthCalculator::SplineLengthCalculator(const std::shared_ptr<Spline> &spline)
    :
      spline(spline),
      maxT(spline->getMaxT()),

      //the spline length will be lazy-computed when it's needed
      atomic_splineLength(-1)
{

}

double SplineLengthCalculator::findLength(double beginT, double endT, bool useShortestPath, double eps) const
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
        double computedLength = computeLength(actualBeginT, actualEndT, eps);

        double splineLength = atomic_splineLength.load(std::memory_order_consume);

        //lazy compute the spline length
        if(splineLength < 0)
        {
            //modifying a variable inside of a const function can be non thread safe, but in this case thread safety is not a concern
            //we're using std::atomic so the worst possible thing that could happen is this variable could be computed multiple times by simultaneous threads
            //but if that happens, it'll be the same both times so who cares
            splineLength = computeLength(0, maxT, eps);
            atomic_splineLength.store(splineLength, std::memory_order_release);
        }

        //that was one direction around, the other direction will be (splineLength - computedLength). return the smaller of the two
        return std::min(computedLength, splineLength - computedLength);
    }
    else
    {
        //make sure beginT is less than endT
        double actualBeginT = std::min(beginT, endT);
        double actualEndT = std::max(beginT,endT);

        return computeLength(actualBeginT, actualEndT, eps);
    }
}

double SplineLengthCalculator::computeLength(double beginT, double endT, double eps) const
{
    double tDistance = endT - beginT;

    //find the number of segments we're going to use
    int numSegments = ceil(tDistance / LENGTH_INTERVAL);

    //now that we know the number of segments, find the actual interval
    double interval = tDistance / numSegments;

    //for each segment, compute the length of a circle arc passing though that segment's endpoints
    auto previousData = spline->getPosition(beginT);
    double lengthSum = 0;

    for(int i = 0; i < numSegments; i++)
    {
        double currentT = beginT + (i + 1) * interval;

        auto currentData = spline->getPosition(currentT);
        lengthSum += computeLengthHelper(
                    currentT - interval,    previousData,
                    currentT,               currentData, eps);

        previousData = currentData;
    }

    return lengthSum;
}

double SplineLengthCalculator::computeLengthHelper(
        double beginT, const Vector3D &beginPosition,
        double endT, const Vector3D &endPosition, double eps) const
{
    //compute the midpoint between the start and end
    double midT = (beginT + endT) * 0.5;
    Vector3D midPosition = spline->getPosition(midT);

    //compute the length squared for start->mid plus mid->end
    double fullLength = (endPosition - beginPosition).lengthSquared();
    double halfLengths = ((endPosition - midPosition).lengthSquared() + (midPosition - beginPosition).lengthSquared()) * 2;

    //if the difference bween the total lengh squared and the length squared for the two halves is small enough, return the length of the whole segment
    //eps is a percentage
    double test = fullLength / halfLengths;
    double limit = 1 - eps;
    if(fullLength / halfLengths > limit || (endT - beginT < RECURSIVE_MINIMUM_INTERVAL))
    {
        return sqrt(fullLength);
    }
    else
    {
        //there is too much of a difference, so recursively compute the lengths
        double firstHalf = computeLengthHelper(beginT, beginPosition, midT, midPosition, eps);
        double secondHalf = computeLengthHelper(midT, midPosition, endT, endPosition, eps);

        return firstHalf + secondHalf;
    }
}
