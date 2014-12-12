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
    //this method is destructive to avoid having to make copies. the contents of the input arrays is undefined after the method returns.
    template<class T>
    static std::vector<T> solveSymmetricTridiagonal(std::vector<double> &mainDiagonal, std::vector<double> &secondaryDiagonal, std::vector<T> &inputVector);

    //solve the given cyclic tridiagonal matrix system, with the assumption that the lower diagonal and upper diagonal (ie secondaryDiagonal) are identical
    //in other words, assume that the matrix is symmetric
    //this method is destructive to avoid having to make copies. the contents of the input arrays is undefined after the method returns.
    template<class T>
    static std::vector<T> solveCyclicSymmetricTridiagonal(std::vector<double> &mainDiagonal, std::vector<double> &secondaryDiagonal, std::vector<T> &inputVector);
};

template<class T>
std::vector<T> LinearAlgebra::solveSymmetricTridiagonal(std::vector<double> &mainDiagonal, std::vector<double> &secondaryDiagonal, std::vector<T> &inputVector)
{
    //use the thomas algorithm to solve the tridiagonal matrix
    // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    std::vector<T> outputVector(inputVector.size());

    //forward sweep
    for(size_t i = 1; i < inputVector.size(); i++)
    {
        double m = secondaryDiagonal.at(i - 1) / mainDiagonal.at(i - 1);
        mainDiagonal[i] -= m * secondaryDiagonal.at(i - 1);
        inputVector[i] -= m * inputVector.at(i - 1);
    }

    //back substitution
    int finalIndex = inputVector.size() - 1;
    outputVector[finalIndex] = inputVector.at(finalIndex) / mainDiagonal.at(finalIndex);
    for(int i = finalIndex - 1; i >= 0; i--)
    {
        outputVector[i] = (inputVector.at(i) - secondaryDiagonal.at(i) * outputVector.at(i + 1)) / mainDiagonal.at(i);
    }

    return outputVector;
}

template<class T>
std::vector<T> LinearAlgebra::solveCyclicSymmetricTridiagonal(std::vector<double> &mainDiagonal, std::vector<double> &secondaryDiagonal, std::vector<T> &inputVector)
{
    //apply the sherman-morrison algorithm to the cyclic tridiagonal matrix so that we can use the standard tridiagonal algorithm
    //we're getting this algorithm from http://www.cs.princeton.edu/courses/archive/fall11/cos323/notes/cos323_f11_lecture06_linsys2.pdf
    //basically, we're going to solve two different non-cyclic versions of this system and then combine the results

    size_t size = inputVector.size();


    //the value at the upper right and lower left of the input matrix. it's at the end of the secondary diagonal array because almost all
    //cyclic tridiagonal papers treat it as an extension of the secondary diagonals
    double cornerValue = secondaryDiagonal.at(size - 1);

    //gamma value - doesn't affect actual output (the algorithm makes sure it cancels out), but a good choice for this value can reduce floating point errors
    double gamma = -mainDiagonal.at(0);
    double cornerMultiplier = cornerValue/gamma;

    //corrective vector U: should be all 0, except for gamma in the first element, and cornerValue at the end
    std::vector<double> correctionInputU(size);
    correctionInputU[0] = gamma;
    correctionInputU[size - 1] = cornerValue;

    //corrective vector V: should be all 0, except for 1 in the first element, and cornerValue/gamma at the end
    std::vector<double> correctionV(size);
    correctionV[0] = 1;
    correctionV[size - 1] = cornerMultiplier;

    //modify the main diagonal of the matrix to account for the correction vector
    mainDiagonal[0] -= gamma;
    mainDiagonal[size - 1] -= cornerValue * cornerMultiplier;

    //create a copy of the main diangonal, because the non-cyclic tridiagonal algorithm will destroy the contents, but we need them twice
    std::vector<double> mainDiagonalCopy(mainDiagonal);

    //solve the modified system for the input vector
    //NOTE: even though the public api for the noncyclic method says that the secondary diagonal may be modified, we know that it isn't modified
    //so it's safe to pass it in without copying
    std::vector<T> initialOutput = solveSymmetricTridiagonal(mainDiagonal, secondaryDiagonal, inputVector);

    //solve the modified system for the correction vector
    std::vector<double> correctionOutput = solveSymmetricTridiagonal(mainDiagonalCopy, secondaryDiagonal, correctionInputU);

    //compute the corrective T to apply to each initial output
    //this involves a couple dot products, but all of the elements on the correctionV vector are 0 except the first and last
    //so just compute those directly instead of looping through and multplying a bunch of 0s
    //T factor = vectorDotProduct(initialOutput, correctionV) / (1 + vectorDotProduct(correctionV, correctionOutput));
    T factor = (initialOutput.at(0) + initialOutput.at(size - 1) * cornerMultiplier) / (1 + correctionOutput.at(0) + correctionOutput.at(size - 1) * cornerMultiplier);

    //use the correction factor to modify the result
    for(size_t i = 0; i < size; i++)
    {
        initialOutput[i] -= factor * correctionOutput.at(i);
    }

    return initialOutput;
}

#endif // LINEARSOLVER_H
