#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QVector3D>
#include <QTime>

#include "spline_library/basis/uniform_cubic_bspline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"

Benchmarker::Benchmarker(QObject *parent)
    :QObject(parent), bigDistribution(0, 1000), smallDistribution(0,1), smallVarianceDistribution(4, 10)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    QMap<QString, float> results;
    results["uniform_cr[100]"] = timeFunction(&Benchmarker::uniformCR,      1000,  1000000, 100);
    results["uniform_cr[10000]"] = timeFunction(&Benchmarker::uniformCR,    10000, 100000, 100000);
    results["bspline[100]"] = timeFunction(&Benchmarker::bspline,           1000,  1000000, 100);
    results["bspline[10000]"] = timeFunction(&Benchmarker::bspline,         10000, 100000, 100000);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::uniformCR(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running uniform CR, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<UniformCRSpline<QVector2D>>(randomPoints2D_SmallVariance(size));

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float a = randomFloat(max);

            s->getPosition(a);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::bspline(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running bspline, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<CubicHermiteSpline<QVector2D>>(randomPoints2D_SmallVariance(size));

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float a = randomFloat(max);

            s->getPosition(a);
        }
    }
    emit setProgressValue(repeat);
}

std::vector<QVector2D> Benchmarker::randomPoints2D_Uniform(size_t size)
{
    std::vector<QVector2D> result(size);
    for(size_t s = 0; s < size; s++)
    {
        result[s] = QVector2D(bigDistribution(gen), bigDistribution(gen));
    }
    return result;
}

std::vector<QVector2D> Benchmarker::randomPoints2D_SmallVariance(size_t size)
{
    std::vector<QVector2D> result(size);
    result[0] = QVector2D(smallVarianceDistribution(gen), smallVarianceDistribution(gen));
    for(size_t s = 1; s < size; s++)
    {
        result[s] = result[s - 1] + QVector2D(smallVarianceDistribution(gen), smallVarianceDistribution(gen));
    }
    return result;
}

float Benchmarker::randomFloat(float max)
{
    return smallDistribution(gen) * max;
}

