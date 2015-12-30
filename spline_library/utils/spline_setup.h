#ifndef T_CALCULATOR_H
#define T_CALCULATOR_H

#include <unordered_map>
#include <vector>

namespace SplineSetup
{
    //compute the T values for the given points, with the given alpha.
    //the distance in T between adjacent points is the magitude of the distance, raised to the power alpha

    template<class floating_t>
    floating_t computeTForDistanceSq(floating_t distanceSq, floating_t alpha);

    //if padding is > 0, the element whose index is equal to padding wil lbe zero and everything to the left will be negative
    //additionally, the T values will be normalized to a maximum of (points.size() - 2 * padding)
    //the implication of the above two lines is that the leftmost and rightmost elements in the array will be "excluded" from the final calculation,
    //with the padding indicating how many to exclude
    template<class SegmentType, typename floating_t>
    std::unordered_map<int, floating_t> computeTValues(const std::vector<SegmentType> &points, floating_t alpha, int padding);


    template<class InterpolationType, typename floating_t>
    std::unordered_map<int, floating_t> computeBSplineKnots(const std::vector<InterpolationType> &points, floating_t alpha, int splineDimension);

    template<class InterpolationType, typename floating_t>
    std::unordered_map<int, floating_t> computeLoopingBSplineKnots(const std::vector<InterpolationType> &points,
                                                                         floating_t alpha,
                                                                         int padding);


    //compute the T values for the given points, with the given alpha.
    //if padding is zero, this method will return points.size() + 1 points - the "extra" point is because the first point in the list is represented twice
    //once for the beginning and a second time when we wrap around at the end

    //if padding is > 0, this method will also compute "extra" T values before the beginning and after the end
    //these won't actually add any extra information, but may simplify calculations that have to wrap around the loop
    template<class SegmentType, typename floating_t>
    std::unordered_map<int, floating_t> computeLoopingTValues(const std::vector<SegmentType> &points, floating_t alpha, int padding);



    //given a list of spline segments and a t value, return the index of the segment the t value falls within
    //the segment type only needs a 't0' variable denoting the start of the segment
    //and a 't1' variable denoting the end of the segment
    template<class SegmentType, typename floating_t>
    const SegmentType& getSegmentForT(const std::vector<SegmentType> &segmentData, floating_t t);

    //given a list of knots and a t value, return the index of the knot the t value falls within
    template<typename floating_t>
    size_t getIndexForT(const std::vector<floating_t> &knotData, floating_t t);
}

template<class SegmentType, typename floating_t>
const SegmentType& SplineSetup::getSegmentForT(const std::vector<SegmentType> &segmentData, floating_t t)
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(t < segmentData[0].t0)
        return segmentData[0];
    if(t > segmentData[segmentData.size() - 1].t1)
        return segmentData[segmentData.size() - 1];

    //perform a binary search on segmentData
    size_t currentMin = 0;
    size_t currentMax = segmentData.size() - 1;
    size_t currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x

    while(t < segmentData[currentIndex].t0 || t > segmentData[currentIndex].t1)
    {
        //if t0 is greater than x, search the left half of the array
        if(segmentData[currentIndex].t0 > t)
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
    return segmentData[currentIndex];
}

template<typename floating_t>
size_t SplineSetup::getIndexForT(const std::vector<floating_t> &knotData, floating_t t)
{
    //we want to find the segment whos t0 and t1 values bound x

    //if no segments bound x, return -1
    if(t < knotData.front())
        return 0;
    if(t >= knotData.back())
        return knotData.size() - 1;

    //perform a binary search on segmentData
    size_t currentMin = 0;
    size_t currentMax = knotData.size() - 2;
    size_t currentIndex = (currentMin + currentMax) / 2;

    //keep looping as long as this segment does not bound x
    while(t < knotData[currentIndex] || t >= knotData[currentIndex + 1])
    {
        //if t0 is greater than x, search the left half of the array
        if(t < knotData[currentIndex])
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

template<class floating_t>
floating_t SplineSetup::computeTForDistanceSq(floating_t distanceSq, floating_t alpha)
{
    //if these points are right on top of each other, don't bother with the power calculation
    if(distanceSq < .0001)
    {
        return 0;
    }
    else
    {
        //multiply alpha by 0.5 so that we tke the square root of distanceSq
        //ie: (distanceSq)^(0.5) = distance
        //and result = distance ^ alpha
        //so result = (distanceSq^alpha)^(0.5) = (distanceSq)^(0.5*alpha)
        //this way we don't have to do a pow AND a sqrt
        return pow(distanceSq, alpha * 0.5);
    }
}

template<class SegmentType, typename floating_t>
std::unordered_map<int, floating_t> SplineSetup::computeTValues(
        const std::vector<SegmentType> &points,
        floating_t alpha,
        int padding)
{
    int size = points.size();
    int endPaddingIndex = size - 1 - padding;
    int desiredMaxT = size - 2 * padding - 1;

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size);

    //we know points[padding] will have a t value of 0
    indexToT[padding] = 0;

    //loop backwards from padding to give the earlier points negative t values
    for(int i = padding - 1; i >= 0; i--)
    {
        //Points inside the padding will not be interpolated
        //so give it a negative t value, so that the first actual point can have a t value of 0
        indexToT[i] = indexToT.at(i + 1) + computeTForDistanceSq((points.at(i) - points.at(i + 1)).lengthSquared(), alpha);
    }

    //compute the t values of the other points
    for(int i = padding + 1; i < size; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTForDistanceSq((points.at(i) - points.at(i - 1)).lengthSquared(), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(endPaddingIndex);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    return indexToT;
}

template<class InterpolationType, typename floating_t>
std::unordered_map<int, floating_t> SplineSetup::computeBSplineKnots(const std::vector<InterpolationType> &points,
                                                                     floating_t alpha,
                                                                     int padding)
{
    int size = points.size();

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size + padding * 2);

    //compute the t values each point
    indexToT[0] = 0;
    for(int i = 1; i < size; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTForDistanceSq((points.at(i) - points.at(i - 1)).lengthSquared(), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(size - 1);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t desiredMaxT = size - 1;
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    //add padding in addition to the points - 2 on each end
    //we calculate the padding by taking the difference in t between the nearest real T values
    for(int i = size; i < padding + size; i++)
    {
        floating_t tDiff = indexToT.at(i - 1) - indexToT.at(i - 2);
        indexToT[i] = indexToT.at(i - 1) + tDiff;
    }
    for(int i = -1; i >= -padding; i--)
    {
        floating_t tDiff = indexToT.at(i + 2) - indexToT.at(i + 1);
        indexToT[i] = indexToT.at(i + 1) - tDiff;
    }

    return indexToT;
}

template<class InterpolationType, typename floating_t>
std::unordered_map<int, floating_t> SplineSetup::computeLoopingBSplineKnots(const std::vector<InterpolationType> &points,
                                                                     floating_t alpha,
                                                                     int padding)
{
    int size = points.size();

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size + padding * 2);

    //compute the t values each point
    indexToT[0] = 0;
    for(int i = 1; i < size + 1; i++)
    {
        indexToT[i] = indexToT.at(i - 1) + computeTForDistanceSq((points.at(i%size) - points.at(i - 1)).lengthSquared(), alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t desiredMaxT = size;
    floating_t multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    //add padding in addition to the points - 2 on each end
    //we calculate the padding by basically wraping the difference in T values
    for(int i = 1; i < padding + 1; i++)
    {
        floating_t tDiff = indexToT.at(i) - indexToT.at(i - 1);
        indexToT[i + size] = indexToT.at(i + size - 1) + tDiff;
    }
    for(int i = -1; i >= -padding; i--)
    {
        floating_t tDiff = indexToT.at(i + size + 1) - indexToT.at(i + size);
        indexToT[i] = indexToT.at(i + 1) - tDiff;
    }

    return indexToT;
}

template<class SegmentType, typename floating_t>
std::unordered_map<int, floating_t> SplineSetup::computeLoopingTValues(
        const std::vector<SegmentType> &points,
        floating_t alpha,
        int padding)
{
    int size = points.size();

    std::unordered_map<int, floating_t> indexToT;
    indexToT.reserve(size + padding * 2);

    //we know points[0] will have a t value of 0
    indexToT[0] = 0;

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -padding; i--)
    {
        floating_t distance = (points.at(i + size) - points.at((i + 1 + size)%size)).length();
        indexToT[i] = indexToT.at(i + 1) - pow(distance, alpha);
    }

    //compute the t values of the other points
    for(int i = 1; i <= size + padding; i++)
    {
        floating_t distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT[i] = indexToT.at(i - 1) + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    floating_t maxTRaw = indexToT.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    floating_t multiplier = size / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    return indexToT;
}

#endif // T_CALCULATOR_H
