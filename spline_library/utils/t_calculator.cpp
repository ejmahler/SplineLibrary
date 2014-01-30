#include "t_calculator.h"

#include <cmath>

TCalculator::TCalculator()
{
}

std::unordered_map<int, double> TCalculator::computeTValues(const std::vector<Vector3D> &points, double alpha, int padding)
{
    std::unordered_map<int, double> indexToT_Raw, indexToT;

    int size = points.size();
    int endPaddingIndex = size - 1 - padding;
    int unpaddedSize = size - 2 * padding;

    //we know points[padding] will have a t value of 0
    indexToT_Raw[padding] = 0;

    //loop backwards from padding to give the earlier points negative t values
    for(int i = padding - 1; i >= 0; i--)
    {
        //points[1] is a control point, so give it a negative t value, so that the first actual point can have a t value of 0
        double distance = (points.at(i) - points.at(i + 1)).length();
        indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
    }

    //compute the t values of the other points
    for(int i = padding + 1; i < size; i++)
    {
        double distance = (points.at(i) - points.at(i - 1)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(endPaddingIndex);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = unpaddedSize * it->second / maxTRaw;
    }

    return indexToT;
}


std::unordered_map<int, double> TCalculator::computeLoopingTValues(const std::vector<Vector3D> &points, double alpha, int padding)
{
    std::unordered_map<int, double> indexToT_Raw, indexToT;

    int size = points.size();

    //we know points[0] will have a t value of 0
    indexToT_Raw[0] = 0;

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -padding; i--)
    {
        double distance = (points.at(i + size) - points.at((i + 1 + size)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i + 1] - pow(distance, alpha);
    }

    //compute the t values of the other points
    for(int i = 1; i <= size + padding; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT_Raw[i] = indexToT_Raw[i - 1] + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT_Raw.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing tem by maxT
    for(auto it = indexToT_Raw.begin(); it != indexToT_Raw.end(); it++)
    {
        indexToT[it->first] = size * it->second / maxTRaw;
    }

    return indexToT;
}
