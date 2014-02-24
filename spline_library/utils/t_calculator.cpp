#include "t_calculator.h"

#include <cmath>

TCalculator::TCalculator()
{
}

std::unordered_map<int, double> TCalculator::computeTValues(const std::vector<Vector3D> &points, double alpha, int padding)
{
    int size = points.size();
    int endPaddingIndex = size - 1 - padding;
    int desiredMaxT = size - 2 * padding - 1;

    std::unordered_map<int, double> indexToT;
    indexToT.reserve(size);

    //we know points[padding] will have a t value of 0
    indexToT[padding] = 0;

    //loop backwards from padding to give the earlier points negative t values
    for(int i = padding - 1; i >= 0; i--)
    {
        //points[1] is a control point, so give it a negative t value, so that the first actual point can have a t value of 0
        double distance = (points.at(i) - points.at(i + 1)).length();
        indexToT[i] = indexToT.at(i + 1) - pow(distance, alpha);
    }

    //compute the t values of the other points
    for(int i = padding + 1; i < size; i++)
    {
        double distance = (points.at(i) - points.at(i - 1)).length();
        indexToT[i] = indexToT.at(i - 1) + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT.at(endPaddingIndex);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    double multiplier = desiredMaxT / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    return indexToT;
}


std::unordered_map<int, double> TCalculator::computeLoopingTValues(const std::vector<Vector3D> &points, double alpha, int padding)
{
    int size = points.size();

    std::unordered_map<int, double> indexToT;
    indexToT.reserve(size + padding * 2);

    //we know points[0] will have a t value of 0
    indexToT[0] = 0;

    //loop backwards from 0 to give the earlier points negative t values
    for(int i = -1; i >= -padding; i--)
    {
        double distance = (points.at(i + size) - points.at((i + 1 + size)%size)).length();
        indexToT[i] = indexToT.at(i + 1) - pow(distance, alpha);
    }

    //compute the t values of the other points
    for(int i = 1; i <= size + padding; i++)
    {
        double distance = (points.at(i%size) - points.at((i - 1)%size)).length();
        indexToT[i] = indexToT.at(i - 1) + pow(distance, alpha);
    }

    //we want to know the t value of the last segment so that we can normalize them all
    float maxTRaw = indexToT.at(size);

    //now that we have all ouf our t values and indexes figured out, normalize the t values by dividing them by maxT
    double multiplier = size / maxTRaw;
    for(auto &entry: indexToT)
    {
        entry.second *= multiplier;
    }

    return indexToT;
}
