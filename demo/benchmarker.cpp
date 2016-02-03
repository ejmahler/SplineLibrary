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

    auto randomGenerator = std::bind(&Benchmarker::randomPoints2D_Uniform, this, std::placeholders::_1);

    QMap<QString, float> results;
    results["uniform_cr[100]"] = timeFunction(&Benchmarker::testPrecision,  "Testing random precision",         100, 1000, 100,randomGenerator);
    results["uniform_cr[1000]"] = timeFunction(&Benchmarker::testPrecision, "Testing random precision",         100, 1000, 1000,randomGenerator);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::testPrecision(QString message, int repeat, int queries, size_t size, std::function<std::vector<VectorT> (size_t)> pointGenerator)
{
    gen.seed(10);

    emit setProgressText(message + QString(", size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        UniformCRSpline<VectorT, double> s(pointGenerator(size));

        double max = s.getMaxT();
        for(int q = 0; q < queries; q++)
        {
            double a = randomFloat(max);
            double b = randomFloat(max);

            s.arcLength(a, b);
        }
    }
    emit setProgressValue(repeat);
}




std::vector<VectorT> Benchmarker::randomPoints2D_Uniform(size_t size)
{
    std::vector<VectorT> result(size);
    for(size_t s = 0; s < size; s++)
    {
        result[s] = VectorT({bigDistribution(gen), bigDistribution(gen)});
    }
    return result;
}

std::vector<VectorT> Benchmarker::randomPoints2D_SmallVariance(size_t size)
{
    std::vector<VectorT> result(size);
    result[0] = VectorT({smallVarianceDistribution(gen), smallVarianceDistribution(gen)});
    for(size_t s = 1; s < size; s++)
    {
        result[s] = result[s - 1] + VectorT({smallVarianceDistribution(gen), smallVarianceDistribution(gen)});
    }
    return result;
}

double Benchmarker::randomFloat(double max)
{
    return smallDistribution(gen) * max;
}

