#pragma once

#include <QObject>
#include <QtTest/QtTest>

#include <memory>
#include <random>

#include "spline_library/vector.h"
#include "spline_library/spline.h"

Q_DECLARE_METATYPE(Vector2)
Q_DECLARE_METATYPE(Vector3)
Q_DECLARE_METATYPE(std::shared_ptr<Spline<Vector2>>)


//perform a linear interpolation between a and b
template<class InterpolationType, typename floating_t>
InterpolationType lerp(InterpolationType a, InterpolationType b, floating_t t)
{
    return a * (floating_t(1) - t) + b * t;
}

//we need to pad out the ends of the data differently depending on spline type
//this way all of the splines will have the same arc length, so it'll be easier to test
template<class T>
std::vector<T> addPadding(std::vector<T> list, size_t paddingSize)
{
    list.reserve(list.size() + paddingSize * 2);
    for(size_t i = 0; i < paddingSize; i++) {
        list.insert(list.begin(), list[0] - (list[1] - list[0]));
    }
    for(size_t i = 0; i < paddingSize; i++) {
        list.push_back(list[list.size() - 1] + (list[list.size() - 1] - list[list.size() - 2]));
    }
    return list;
}

template<class T>
void compareFloatsLenient(T actual, T expected, T tol)
{
    auto error = std::abs(actual - expected) / expected;
    if(error > tol) {
        std::string errorMessage = QString("Compared floats were different. Actual: %1, Expected: %2").arg(QString::number(actual), QString::number(expected)).toStdString();
        QFAIL(errorMessage.data());
    }
};

//we have a bunch of functions to help create test data. the alternative is copy/pasting a list of snarled one-liners into 10 different tests
//we're putting all these functions in a class so we can typedef the whole class,
//so that every infocation doesn't need to supply template parameters
template<typename floating_t>
class SplineCreator
{
public:
    typedef std::shared_ptr<Spline<Vector<2,floating_t>, floating_t>> SplinePtr;
    typedef Vector<2, floating_t> T;

    //several functions that create instances of splines
    static SplinePtr createUniformCR(std::vector<T> data) {
        return std::make_shared<UniformCRSpline<T, floating_t>>(addPadding(data,1));
    }
    static SplinePtr createCatmullRom(std::vector<T> data, floating_t alpha) {
        auto padded = addPadding(data,1);
        return std::make_shared<CubicHermiteSpline<T, floating_t>>(padded, alpha);
    }
    static SplinePtr createCubicHermite(std::vector<T> data, floating_t alpha) {
        auto tangents = makeTangents(data);
        return std::make_shared<CubicHermiteSpline<T, floating_t>>(data, tangents, alpha);
    }
    static SplinePtr createQuinticCatmullRom(std::vector<T> data, floating_t alpha) {
        auto padded = addPadding(data,2);
        return std::make_shared<QuinticHermiteSpline<T, floating_t>>(padded, alpha);
    }
    static SplinePtr createQuinticHermite(std::vector<T> data, floating_t alpha) {
        auto tangents = makeTangents(data);
        auto curves = makeTangents(tangents);
        return std::make_shared<QuinticHermiteSpline<T, floating_t>>(data, tangents, curves, alpha);
    }
    static SplinePtr createNatural(std::vector<T> data, bool includeEndpoints, floating_t alpha) {
        if(!includeEndpoints) {
            data = addPadding(data, 1);
        }
        return std::make_shared<NaturalSpline<T, floating_t>>(data, includeEndpoints, alpha);
    }
    static SplinePtr createNotAKnot(std::vector<T> data, bool includeEndpoints, floating_t alpha) {
        if(!includeEndpoints) {
            data = addPadding(data, 1);
        }
        return std::make_shared<NaturalSpline<T, floating_t>>(data, includeEndpoints, alpha, NaturalSpline<Vector<2,floating_t>, floating_t>::NotAKnot);
    }
    static SplinePtr createUniformBSpline(std::vector<T> data) {
        return std::make_shared<UniformCubicBSpline<T, floating_t>>(addPadding(data,1));
    }
    static SplinePtr createGenericBSpline(std::vector<T> data, size_t degree) {
        auto padded = addPadding(data, (degree - 1)/2);
        return std::make_shared<GenericBSpline<T, floating_t>>(padded, degree);
    }



    //several functions that generate data points
    static std::vector<T> generateRandomData(size_t size, size_t seed=10) {
        std::minstd_rand gen;
        gen.seed(seed);
        std::uniform_real_distribution<floating_t> distribution(2,5);

        std::vector<T> result(size);
        result[0] = T({distribution(gen), distribution(gen)});
        for(size_t i = 1; i < size; i++)
        {
            result[i] = result[i - 1] + T({distribution(gen), distribution(gen)});
        }
        return result;
    }

    static std::vector<Vector<2, floating_t>> generateTriangleNumberData(size_t size) {
        std::vector<Vector<2, floating_t>> result(size);
        size_t currentTriangleNumber = 0;
        for(size_t i = 0; i < size; i++)
        {
            currentTriangleNumber += i;
            floating_t value = currentTriangleNumber;
            result[i] = T({value, value});
        }
        return result;
    }
private:

    //Given a list of points, compute an equal-sized list of tangents to use in a cubic or quintic hermite spline
    //use the finite difference algorithm
    static std::vector<T> makeTangents(std::vector<T> points) {
        std::vector<T> tangents(points.size());

        //one-sided difference at the start
        tangents[0] = points[1] - points[0];

        //normal finite difference in the middle
        for(size_t i = 1; i < points.size() - 1; i++)
        {
            tangents[i] = 0.5f * ((points[i+1] - points[i]) + (points[i] - points[i-1]));
        }

        //one-sided difference at the start
        tangents[tangents.size() - 1] = points[points.size() - 1] - points[points.size() - 2];
        return tangents;
    }
};

typedef SplineCreator<float> TestDataFloat;
typedef SplineCreator<float> TestDataDouble;





