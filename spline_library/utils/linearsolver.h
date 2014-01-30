#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
//#include <armadillo>

class Vector3D;

class LinearSolver
{
private:
    LinearSolver();

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
    //dot product methods: we need to implement one of these for every "cyclic symmetric tridiagonal" input types, plus one for doubles

    //given two arrays of doubles, compute the dot product of the two arrays
    static double vectorDotProduct(const std::vector<double> &left, const std::vector<double> &right);

    //given an array of doubles and an array of vectors, compute the dot product of the doubles with each of the components of the vectors
    //ie, output.x will be inputDouble (dot) inputvector.x, output.y will be inputDouble (dot) inputvector.y, output.z will be inputDouble (dot) inputvector.z
    static Vector3D vectorDotProduct(const std::vector<double> &left, const std::vector<Vector3D> &right);
};

template<class T>
std::vector<T> LinearSolver::solveSymmetricTridiagonal(const std::vector<double> &mainDiagonalReadOnly, const std::vector<double> &secondaryDiagonalReadOnly, const std::vector<T> &inputVectorReadOnly)
{
    //use the thomas algorithm to solve the tridiagonal matrix
    // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    std::vector<T> outputVector(inputVectorReadOnly.size());

    //create modifiable copies of the input vectors
    std::vector<double> mainDiagonal(mainDiagonalReadOnly);
    std::vector<T> inputVector(inputVectorReadOnly);

    //forward sweep
    for(int i = 2; i < inputVector.size(); i++)
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
std::vector<T> LinearSolver::solveCyclicSymmetricTridiagonal(const std::vector<double> &mainDiagonal, const std::vector<double> &secondaryDiagonal, const std::vector<T> &inputVector)
{
    //apply the sherman-morrison algorithm to the cyclic tridiagonal matrix so that we can use the standard tridiagonal algorithm
    //we're getting this algorithm from http://www.cs.princeton.edu/courses/archive/fall11/cos323/notes/cos323_f11_lecture06_linsys2.pdf
    //basically, we're going to solve two different non-cyclic versions of this system and then combine the results

    int size = inputVector.size();


    //the value at the upper right and lower left of the input matrix. it's at the end of the secondary diagonal array because almost all
    //cyclic tridiagonal papers treat it as an extension of the secondary diagonals
    double cornerValue = secondaryDiagonal.at(size - 1);

    //gamma value - doesn't affect actual output (the algorithm makes sure it cancels out), but a good choice for this value can reduce floating point errors
    double gamma = -mainDiagonal[0]*100;

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
    T factor = vectorDotProduct(correctionV, initialOutput) / (1 + vectorDotProduct(correctionV, correctionOutput));

    //use the correction input and output to modify the result
    for(int i = 0; i < size; i++)
    {
        initialOutput[i] -= factor * correctionOutput[i];
    }

    /*
    //since we know the algorithm is wrong, use armidillo to check the actual result against wat the algorithm has returned
    arma::mat matrix(size, size);
    matrix.fill(0);
    for(int i = 0; i < size; i++)
    {
        matrix.at(i,i) = mainDiagonal.at(i);
        matrix.at(i, (i + 1)%size) = secondaryDiagonal.at(i);
        matrix.at((i + 1)%size, i) = secondaryDiagonal.at(i);
    }

    arma::vec xVector(size);
    for(int i = 0; i < size; i++)
    {
        xVector.at(i) = inputVector.at(i).x();
    }

    arma::vec yVector(size);
    for(int i = 0; i < size; i++)
    {
        yVector.at(i) = inputVector.at(i).y();
    }

    arma::mat inverse = matrix.i();
    arma::vec xOutput = inverse * xVector;
    arma::vec yOutput = inverse * yVector;

    std::vector<T> outputVector2(size);
    for(int i = 0; i < size; i++)
    {
        outputVector2[i] = T(xOutput.at(i), yOutput.at(i), 0);
    }


    T out1 = initialOutput[0];
    T out2 = outputVector2[0];
    */

    return initialOutput;
}

#endif // LINEARSOLVER_H
