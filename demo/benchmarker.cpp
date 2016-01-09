#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QVector3D>
#include <QTime>

#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/uniform_cr_spline.h"

Benchmarker::Benchmarker(QObject *parent)
    :QObject(parent), bigDistribution(0, 1000), smallDistribution(0,1)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    QMap<QString, float> results;
    results["Cubic Hermite[100]"] = timeFunction(&Benchmarker::cubicHermiteQuery,   10000,  100000, 100);
    results["Cubic Hermite[10000]"] = timeFunction(&Benchmarker::cubicHermiteQuery, 10000,  10000, 100000);
    results["Uniform CR[100]"] = timeFunction(&Benchmarker::uniformCRQuery,   10000,  100000, 100);
    results["Uniform CR[100000]"] = timeFunction(&Benchmarker::uniformCRQuery,10000,  10000, 100000);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::cubicHermiteQuery(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running cubic hermite queries, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<CubicHermiteSpline<QVector2D>>(randomPoints2D_Uniform(size));

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::uniformCRQuery(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running uniform CR queries, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<UniformCRSpline<QVector2D, double>>(randomPoints2D_Uniform(size));

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat);
}

std::vector<QVector3D> Benchmarker::randomPoints3D_Uniform(size_t size)
{
    std::vector<QVector3D> result(size);
    for(size_t s = 0; s < size; s++)
    {
        result[s] = QVector3D(bigDistribution(gen), bigDistribution(gen), bigDistribution(gen));
    }
    return result;
}

std::vector<Vector3D> Benchmarker::randomPoints3DDouble_Uniform(size_t size)
{
    std::vector<Vector3D> result(size);
    for(size_t s = 0; s < size; s++)
    {
        result[s] = Vector3D(bigDistribution(gen), bigDistribution(gen), bigDistribution(gen));
    }
    return result;
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

std::vector<QVector2D> Benchmarker::randomPoints2D_Unbalanced(size_t size)
{
    std::vector<QVector2D> result(size);
    for(size_t s = 0; s < size; s++)
    {
        if(s < size/2)
            result[s] = QVector2D(smallDistribution(gen), smallDistribution(gen));
        else
            result[s] = QVector2D(bigDistribution(gen), bigDistribution(gen));
    }
    return result;
}

float Benchmarker::randomFloat(float max)
{
    return smallDistribution(gen) * max;
}

