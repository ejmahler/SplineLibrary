#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QVector3D>
#include <QTime>

#include "spline_library/utils/arclength.h"
#include "spline_library/splines/generic_b_spline.h"
#include "spline_library/splines/uniform_cr_spline.h"
#include "spline_library/splines/natural_spline.h"

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
        auto points = randomPoints_Uniform<VectorT, D, FloatingT>(randomSource, size);
        std::unique_ptr<SplineType> result = std::make_unique<LoopingUniformCRSpline<VectorT,FloatingT>>(points);
        return result;
    };

    auto genericBSpline = [this, randomSource](size_t size) {
        auto points = randomPoints_Uniform<VectorT, D, FloatingT>(randomSource, size);
        std::unique_ptr<SplineType> result = std::make_unique<LoopingGenericBSpline<VectorT,FloatingT>>(points, 7);
        return result;
    };

    QMap<QString, float> results;
    timeSplineMemberFunction(results, &Benchmarker::testArcLength, crSpline, "uniform_cr[10]",    10000, 12);
    timeSplineMemberFunction(results, &Benchmarker::testArcLength, crSpline, "uniform_cr[1000]",  1000, 1002);
    timeSplineMemberFunction(results, &Benchmarker::testArcLength, genericBSpline, "bspline[10]",    1000, 16);
    timeSplineMemberFunction(results, &Benchmarker::testArcLength, genericBSpline, "bspline[1000]",  100, 1006);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::timeSplineMemberFunction(
        QMap<QString, float>& results,
        void(Benchmarker::*testFunction)(int, const SplineType&) ,
        std::function<std::unique_ptr<SplineType>(size_t)> splineFunction,
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
        (this->*testFunction)(queries, *spline);
        totalElapsed += t.elapsed();
    }
    emit setProgressValue(repeats);

    results[message] = 1000 * float(totalElapsed) / (repeats * queries);
}

void Benchmarker::testArcLength(int queries, const LoopingSpline<VectorT, FloatingT> &spline)
{
    std::uniform_real_distribution<FloatingT> dist(0, spline.getMaxT());

    for(int q = 0; q < queries; q++)
    {
        FloatingT a = dist(gen);
        FloatingT b = dist(gen);

        float result = spline.cyclicArcLength(a,b);
    }
}

