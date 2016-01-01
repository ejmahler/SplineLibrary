#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>

class Vector3D;

class LinearAlgebra
{
private:
    LinearAlgebra() = default;

public:

    //solve the given tridiagonal matrix system, with the assumption that the lower diagonal and upper diagonal (ie secondaryDiagonal) are identical.
    //in other words, assume that the matrix is symmetric
    template<class OutputType, typename floating_t>
    static std::vector<OutputType> solveSymmetricTridiagonal(

            std::vector<floating_t> mainDiagonal,
            std::vector<floating_t> secondaryDiagonal,
            std::vector<OutputType> inputVector);

    //solve the given tridiagonal matrix system
    template<class OutputType, typename floating_t>
    static std::vector<OutputType> solveTridiagonal(

            std::vector<floating_t> mainDiagonal,
            std::vector<floating_t> upperDiagonal,
            std::vector<floating_t> lowerDiagonal,
            std::vector<OutputType> inputVector);

    //solve the given cyclic tridiagonal matrix system, with the assumption that the lower diagonal and upper diagonal (ie secondaryDiagonal) are identical
    //in other words, assume that the matrix is symmetric
    template<class OutputType, typename floating_t>
    static std::vector<OutputType> solveCyclicSymmetricTridiagonal(

            std::vector<floating_t> mainDiagonal,
            std::vector<floating_t> secondaryDiagonal,
            std::vector<OutputType> inputVector);
};

template<class OutputType, typename floating_t>
std::vector<OutputType> LinearAlgebra::solveTridiagonal(

        std::vector<floating_t> mainDiagonal,
        std::vector<floating_t> upperDiagonal,
        std::vector<floating_t> lowerDiagonal,
        std::vector<OutputType> inputVector)
{
    //use the thomas algorithm to solve the tridiagonal matrix
    // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

    //forward sweep
    upperDiagonal[0] /= mainDiagonal[0];
    inputVector[0] /= mainDiagonal[0];
    for(size_t i = 1; i < inputVector.size(); i++)
    {
        upperDiagonal[i] /= (mainDiagonal[i] - lowerDiagonal[i - 1] * upperDiagonal[i - 1]);
        inputVector[i] = (inputVector[i] - lowerDiagonal[i - 1] * inputVector[i - 1]) /
                (mainDiagonal[i] - lowerDiagonal[i - 1] * upperDiagonal[i - 1]);
    }

    //back substitution
    size_t finalIndex = inputVector.size() - 1;
    for(size_t i = finalIndex - 1; i >= 0; i--)
    {
        inputVector[i] -= upperDiagonal[i] * inputVector[i + 1];
    }

    return inputVector;
}

template<class OutputType, typename floating_t>
std::vector<OutputType> LinearAlgebra::solveSymmetricTridiagonal(

        std::vector<floating_t> mainDiagonal,
        std::vector<floating_t> secondaryDiagonal,
        std::vector<OutputType> inputVector)
{
    //use the thomas algorithm to solve the tridiagonal matrix
    // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    std::vector<OutputType> outputVector(inputVector.size());

    //forward sweep
    for(size_t i = 1; i < inputVector.size(); i++)
    {
        floating_t m = secondaryDiagonal.at(i - 1) / mainDiagonal.at(i - 1);
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

template<class OutputType, typename floating_t>
std::vector<OutputType> LinearAlgebra::solveCyclicSymmetricTridiagonal(

        std::vector<floating_t> mainDiagonal,
        std::vector<floating_t> secondaryDiagonal,
        std::vector<OutputType> inputVector)
{
    //apply the sherman-morrison algorithm to the cyclic tridiagonal matrix so that we can use the standard tridiagonal algorithm
    //we're getting this algorithm from http://www.cs.princeton.edu/courses/archive/fall11/cos323/notes/cos323_f11_lecture06_linsys2.pdf
    //basically, we're going to solve two different non-cyclic versions of this system and then combine the results

    size_t size = inputVector.size();


    //the value at the upper right and lower left of the input matrix. it's at the end of the secondary diagonal array because almost all
    //cyclic tridiagonal papers treat it as an extension of the secondary diagonals
    floating_t cornerValue = secondaryDiagonal.at(size - 1);

    //gamma value - doesn't affect actual output (the algorithm makes sure it cancels out), but a good choice for this value can reduce floating point errors
    floating_t gamma = -mainDiagonal.at(0);
    floating_t cornerMultiplier = cornerValue/gamma;

    //corrective vector U: should be all 0, except for gamma in the first element, and cornerValue at the end
    std::vector<floating_t> correctionInputU(size);
    correctionInputU[0] = gamma;
    correctionInputU[size - 1] = cornerValue;

    //modify the main diagonal of the matrix to account for the correction vector
    mainDiagonal[0] -= gamma;
    mainDiagonal[size - 1] -= cornerValue * cornerMultiplier;

    //solve the modified system for the input vector
    std::vector<OutputType> initialOutput = solveSymmetricTridiagonal(
                mainDiagonal,
                secondaryDiagonal,
                std::move(inputVector)
                );

    //solve the modified system for the correction vector
    std::vector<floating_t> correctionOutput = solveSymmetricTridiagonal(
                std::move(mainDiagonal),
                std::move(secondaryDiagonal),
                std::move(correctionInputU)
                );

    //compute the corrective OutputType to apply to each initial output
    //this involves a couple dot products, but all of the elements on the correctionV vector are 0 except the first and last
    //so just compute those directly instead of looping through and multplying a bunch of 0s
    OutputType factor = (initialOutput.at(0) + initialOutput.at(size - 1) * cornerMultiplier) / (1 + correctionOutput.at(0) + correctionOutput.at(size - 1) * cornerMultiplier);

    /*std::vector<floating_t> correctionV(size);
    correctionV[0] = 1;
    correctionV[size - 1] = cornerMultiplier;
    OutputType factor = vectorDotProduct(initialOutput, correctionV) / (1 + vectorDotProduct(correctionV, correctionOutput));*/

    //use the correction factor to modify the result
    for(size_t i = 0; i < size; i++)
    {
        initialOutput[i] -= factor * correctionOutput.at(i);
    }

    return initialOutput;
}

#endif // LINEARSOLVER_H
