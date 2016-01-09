#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QVector3D>
#include <QTime>

#include "spline_library/natural/natural_spline.h"

Benchmarker::Benchmarker(QObject *parent)
    :QObject(parent), bigDistribution(0, 1000), smallDistribution(0,1)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    QMap<QString, float> results;
    results["Natural 3D[f][100]"] = timeFunction(&Benchmarker::naturalSpline3DQuery,   10000,  100000, 100);
    results["Natural 3D[f][10000]"] = timeFunction(&Benchmarker::naturalSpline3DQuery, 10000,   100000, 10000);
    results["Natural 3D[d][100]"] = timeFunction(&Benchmarker::naturalSplineDouble3DQuery,   10000,  100000, 100);
    results["Natural 3D[d][10000]"] = timeFunction(&Benchmarker::naturalSplineDouble3DQuery, 10000,   100000, 10000);
    results["Natural 2D[f][100]"] = timeFunction(&Benchmarker::naturalSpline2DQuery,   10000,  100000, 100);
    results["Natural 2D[f][10000]"] = timeFunction(&Benchmarker::naturalSpline2DQuery, 10000,   100000, 10000);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::naturalSpline3DQuery(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running 3D natural spline, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<NaturalSpline<QVector3D>>(randomPoints3D_Uniform(size), 0);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getWiggle(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::naturalSplineDouble3DQuery(int repeat, int queries,  size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running 3D double natural spline, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<NaturalSpline<Vector3D, double>>(randomPoints3DDouble_Uniform(size), 0);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getWiggle(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::naturalSpline2DQuery(int repeat, int queries, size_t size)
{
    gen.seed(10);

    emit setProgressText(QString("Running 2D natural spline, size ") + QString::number(size));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        auto s = std::make_shared<NaturalSpline<QVector2D>>(randomPoints2D_Uniform(size), 0);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getWiggle(t);
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

