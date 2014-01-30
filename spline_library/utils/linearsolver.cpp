#include "linearsolver.h"

#include "../vector3d.h"

LinearSolver::LinearSolver()
{
}

//given two arrays of doubles, compute the dot product of the two arrays
double LinearSolver::vectorDotProduct(const std::vector<double> &left, const std::vector<double> &right)
{
    double sum = 0;
    for(int i = 0; i < left.size(); i++)
    {
        sum += left.at(i) * right.at(i);
    }
    return sum;
}

//given an array of doubles and an array of vectors, compute the dot product of the doubles with each of the components of the vectors
//ie, output.x will be inputDouble (dot) inputvector.x, output.y will be inputDouble (dot) inputvector.y, output.z will be inputDouble (dot) inputvector.z
Vector3D LinearSolver::vectorDotProduct(const std::vector<double> &left, const std::vector<Vector3D> &right)
{
    Vector3D sum;
    for(int i = 0; i < left.size(); i++)
    {
        sum += left.at(i) * right.at(i);
    }
    return sum;
}
