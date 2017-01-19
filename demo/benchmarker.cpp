#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QVector3D>
#include <QTime>

#include "spline_library/arclength.h"
#include "spline_library/basis/generic_b_spline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"
#include "spline_library/natural/natural_spline.h"

#include <boost/math/tools/roots.hpp>

Benchmarker::Benchmarker(QObject *parent)
    :QObject(parent), repeats(100)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    //create a distribution of random numbers, and a parameterles function to generate them
    std::uniform_real_distribution<FloatingT> distribution(10,15);
    auto randomSource = [this, &distribution]() {
        return distribution(this->gen);
    };

    //create the functions used to generate the splines
    auto crSpline = [this, randomSource](size_t size) {
        auto points = randomPoints_SmallVariance<VectorT, D, FloatingT>(randomSource, size);
        SplinePtr result = std::make_unique<UniformCRSpline<VectorT,FloatingT>>(points);
        return result;
    };

    auto genericBSpline = [this, randomSource](size_t size) {
        auto points = randomPoints_SmallVariance<VectorT, D, FloatingT>(randomSource, size);
        SplinePtr result = std::make_unique<GenericBSpline<VectorT,FloatingT>>(points, 7);
        return result;
    };

    auto natural = [this, randomSource](size_t size) {
        auto points = randomPoints_SmallVariance<VectorT, D, FloatingT>(randomSource, size);
        SplinePtr result = std::make_unique<NaturalSpline<VectorT,FloatingT>>(points, true, 0.5f);
        return result;
    };

    QMap<QString, float> results;
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, crSpline, "FAST: uniform_cr[10]",    10000, 12);
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, crSpline, "FAST: uniform_cr[1000]",  1000, 1002);
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, genericBSpline, "FAST: bspline[10]",    1000, 16);
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, genericBSpline, "FAST: bspline[1000]",  100, 1006);
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, natural, "FAST: natural[10]",    10000, 10);
    timeSplineMemberFunction(results, &Benchmarker::testSolveArcLength, natural, "FAST: natural[1000]",  1000, 1000);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::timeSplineMemberFunction(
        QMap<QString, float>& results,
        void(Benchmarker::*testFunction)(QString, int, const Spline<VectorT, FloatingT>&) ,
        std::function<SplinePtr(size_t)> splineFunction,
        QString message, int queries, size_t size) {

    emit setProgressText(message);
    emit setProgressRange(0, repeats);

    int totalElapsed = 0;
    gen.seed(10);
    QTime t;

    for(int i = 0; i < repeats; i++) {
        if(canceled) return;

        emit setProgressValue(i);
        auto spline = splineFunction(size);

        t.start();
        (this->*testFunction)(message, queries, *spline.get());
        totalElapsed += t.elapsed();
    }
    emit setProgressValue(repeats);

    results[message] = 1000 * float(totalElapsed) / (repeats * queries);
}

void Benchmarker::testSolveArcLength(QString message, int queries, const Spline<VectorT, FloatingT> &spline)
{
    std::uniform_real_distribution<FloatingT> lengthDist(1, spline.totalLength()/2);

    for(int q = 0; q < queries; q++)
    {
        //set up the test
        FloatingT desired = lengthDist(gen);

        ArcLength::partition(spline, desired);
    }
}

