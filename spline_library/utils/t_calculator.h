#ifndef T_CALCULATOR_H
#define T_CALCULATOR_H

#include <unordered_map>
#include <vector>
#include "../vector3d.h"

class TCalculator
{
private:
    TCalculator();

public:

    //compute the T values for the given points, with the given alpha.
    //the distance in T between adjacent points is the magitude of the distance, raised to the power alpha

    //if padding is > 0, the element whose index is equal to padding wil lbe zero and everything to the left will be negative
    //additionally, the T values will be normalized to a maximum of (points.size() - 2 * padding)
    //the implication of the above two lines is that the leftmost and rightmost elements in the array will be "excluded" from the final calculation,
    //with the padding indicating how many to exclude
    static std::unordered_map<int, double> computeTValues(const std::vector<Vector3D> &points, double alpha, int padding);

    //compute the T values for the given points, with the given alpha.
    //if padding is zero, this method will return points.size() + 1 points - the "extra" point is because the first point in the list is represented twice
    //once for the beginning and a second time when we wrap around at the end

    //if padding is > 0, this method will also compute "extra" T values before the beginning and after the end
    //these won't actually add any extra information, but may simplify calculations that have to wrap around the loop
    static std::unordered_map<int, double> computeLoopingTValues(const std::vector<Vector3D> &points, double alpha, int padding);
};

#endif // T_CALCULATOR_H
