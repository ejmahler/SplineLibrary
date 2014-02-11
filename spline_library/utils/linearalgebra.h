#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>

class Vector3D;

class LinearAlgebra
{
private:
    LinearAlgebra();

public:

    //solve the given tridiagonal matrix system, with the assumption that the lower diagonal and upper diagonal (ie secondaryDiagonal) are identical.
    //in other words, assume that the matrix is symmetric
    template<class T>
    static std::vector<T> solveSymmetricTridiagonal(const std::vector<double> &mainDiagonal, const std::vector<double> &secondaryDiagonal, const std::vector<T> &inputVector);

    //solve the given cyclic tridiagonal matrix system, with the assumption that the lower diagonal and upper diagonal (ie secondaryDiagonal) are identical
    //in other words, assume that the matrix is symmetric
    template<class T>
    static std::vector<T> solveCyclicSymmetricTridiagonal(const std::vector<double> &mainDiagonal, const std::vector<double> &secondaryDiagonal, const std::vector<T> &inputVector);

private:
    //compute the dot product of the two arrays of data
    template<class S, class T>
    static S vectorDotProduct(const std::vector<S> &left, const std::vector<T> &right);
};

template<class T>
std::vector<T> LinearAlgebra::solveSymmetricTridiagonal(const std::vector<double> &mainDiagonalReadOnly, const std::vector<double> &secondaryDiagonalReadOnly, const std::vector<T> &inputVectorReadOnly)
{
    //use the thomas algorithm to solve the tridiagonal matrix
    // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    std::vector<T> outputVector(inputVectorReadOnly.size());

    //create modifiable copies of the input vectors
    std::vector<double> mainDiagonal(mainDiagonalReadOnly);
    std::vector<T> inputVector(inputVectorReadOnly);

    //forward sweep
    for(int i = 1; i < inputVector.size(); i++)
    {
        double m = secondaryDiagonalReadOnly.at(i - 1) / mainDiagonal.at(i - 1);
        mainDiagonal[i] -= m * secondaryDiagonalReadOnly.at(i - 1);
        inputVector[i] -= m * inputVector.at(i - 1);
    }

    //back substitution
    int finalIndex = inputVector.size() - 1;
    outputVector[finalIndex] = inputVector.at(finalIndex) / mainDiagonal.at(finalIndex);
    for(int i = finalIndex - 1; i >= 0; i--)
    {
        outputVector[i] = (inputVector.at(i) - secondaryDiagonalReadOnly.at(i) * outputVector.at(i + 1)) / mainDiagonal.at(i);
    }

    return outputVector;
}

template<class T>
std::vector<T> LinearAlgebra::solveCyclicSymmetricTridiagonal(const std::vector<double> &mainDiagonal, const std::vector<double> &secondaryDiagonal, const std::vector<T> &inputVector)
{
    //apply the sherman-morrison algorithm to the cyclic tridiagonal matrix so that we can use the standard tridiagonal algorithm
    //we're getting this algorithm from http://www.cs.princeton.edu/courses/archive/fall11/cos323/notes/cos323_f11_lecture06_linsys2.pdf
    //basically, we're going to solve two different non-cyclic versions of this system and then combine the results

    int size = inputVector.size();


    //the value at the upper right and lower left of the input matrix. it's at the end of the secondary diagonal array because almost all
    //cyclic tridiagonal papers treat it as an extension of the secondary diagonals
    double cornerValue = secondaryDiagonal.at(size - 1);

    //gamma value - doesn't affect actual output (the algorithm makes sure it cancels out), but a good choice for this value can reduce floating point errors
    double gamma = -mainDiagonal[0];

    //corrective vector U: should be all 0, except for gamma in the first element, and cornerValue at the end
    std::vector<double> correctionInputU(size);
    correctionInputU[0] = gamma;
    correctionInputU[size - 1] = cornerValue;

    //corrective vector V: should be all 0, except for 1 in the first element, and cornerValue/gamma at the end
    std::vector<double> correctionV(size);
    correctionV[0] = 1;
    correctionV[size - 1] = cornerValue/gamma;

    //modify the main diagonal of the matrix to account for the correction vector
    std::vector<double> modifiedMainDiagonal(mainDiagonal);
    modifiedMainDiagonal[0] -= gamma;
    modifiedMainDiagonal[size - 1] -= cornerValue * cornerValue / gamma;

    //solve the modified system for the input vector
    std::vector<T> initialOutput = solveSymmetricTridiagonal(modifiedMainDiagonal, secondaryDiagonal, inputVector);

    //solve the modified system for the correction vector
    std::vector<double> correctionOutput = solveSymmetricTridiagonal(modifiedMainDiagonal, secondaryDiagonal, correctionInputU);

    //compute the corrective T to apply to each initial output
    T factor = vectorDotProduct(initialOutput, correctionV) / (1 + vectorDotProduct(correctionV, correctionOutput));

    //use the correction factor to modify the result
    for(int i = 0; i < size; i++)
    {
        initialOutput[i] -= factor * correctionOutput.at(i);
    }

    return initialOutput;
}

//given two arrays of doubles, compute the dot product of the two arrays
template<class S, class T>
S LinearAlgebra::vectorDotProduct(const std::vector<S> &left, const std::vector<T> &right)
{
    S sum = S();
    for(int i = 0; i < left.size(); i++)
    {
        sum += left.at(i) * right.at(i);
    }
    return sum;
}

#endif // LINEARSOLVER_H
