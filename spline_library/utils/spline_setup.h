#ifndef T_CALCULATOR_H
#define T_CALCULATOR_H

#include <unordered_map>
#include <vector>
#include "../vector3d.h"

namespace SplineSetup
{
    //compute the T values for the given points, with the given alpha.
    //the distance in T between adjacent points is the magitude of the distance, raised to the power alpha

    //if padding is > 0, the element whose index is equal to padding wil lbe zero and everything to the left will be negative
    //additionally, the T values will be normalized to a maximum of (points.size() - 2 * padding)
    //the implication of the above two lines is that the leftmost and rightmost elements in the array will be "excluded" from the final calculation,
    //with the padding indicating how many to exclude
    std::unordered_map<int, double> computeTValues(const std::vector<Vector3D> &points, double alpha, int padding);

    //compute the T values for the given points, with the given alpha.
    //if padding is zero, this method will return points.size() + 1 points - the "extra" point is because the first point in the list is represented twice
    //once for the beginning and a second time when we wrap around at the end

    //if padding is > 0, this method will also compute "extra" T values before the beginning and after the end
    //these won't actually add any extra information, but may simplify calculations that have to wrap around the loop
    std::unordered_map<int, double> computeLoopingTValues(const std::vector<Vector3D> &points, double alpha, int padding);




    //given a list of spline segments and a t value, return the index of the segment the t value falls within
    //the segment type only needs a 't0' variable denoting the start of the segment
    //and a 't1' variable denoting the end of the segment
    template<class SegmentType>
    const SegmentType& getSegmentForT(const std::vector<SegmentType> &segmentData, double t)
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
}

#endif // T_CALCULATOR_H
