#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QTime>

#include "spline_library/natural/natural_spline.h"
#include "spline_library/basis/cubic_b_spline.h"
#include "spline_library/basis/generic_b_spline.h"

Benchmarker::Benchmarker(QObject *parent)
    :QObject(parent), bigDistribution(0, 1000), smallDistribution(0,1)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    QMap<QString, float> results;
    results["unBalanced Natural[100]"] = timeFunction(&Benchmarker::naturalSplineQueryUnbalanced,   10000,  100000, 100, 1);
    results["unBalanced Natural"] = timeFunction(&Benchmarker::naturalSplineQueryUnbalanced,        1000,   100000, 10000, 1);
    results["Balanced Natural[100]"] = timeFunction(&Benchmarker::naturalSplineQueryBalanced,       10000,  100000, 100, 1);
    results["Balanced Natural"] = timeFunction(&Benchmarker::naturalSplineQueryBalanced,            1000,   100000, 10000, 1);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::naturalSplineQueryBalanced(int repeat, int queries, size_t size, float alpha)
{
    gen.seed(10);

    emit setProgressText(QString("Running Balanced Natural With Alpha") + QString::number(alpha));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        std::vector<QVector2D> inputPoints = randomPoints2D_Uniform(size);
        auto s = std::make_shared<NaturalSpline<QVector2D>>(inputPoints, false, alpha);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::naturalSplineQueryUnbalanced(int repeat, int queries, size_t size, float alpha)
{
    gen.seed(10);

    emit setProgressText(QString("Running Unbalanced Natural With Alpha") + QString::number(alpha));
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        std::vector<QVector2D> inputPoints = randomPoints2D_Unbalanced(size);
        auto s = std::make_shared<NaturalSpline<QVector2D>>(inputPoints, false, alpha);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::cubicBSplineQuery(int repeat, int queries, size_t size)
{
    gen.seed(10);

    emit setProgressText("Running Cubic B-Spline Queries");
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        std::vector<QVector2D> inputPoints = randomPoints2D_Uniform(size);
        auto s = std::make_shared<CubicBSpline<QVector2D>>(inputPoints);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat);
}

void Benchmarker::genericBSplineQuery(int repeat, int queries, size_t size)
{
    gen.seed(10);

    emit setProgressText("Running Generic B-Spline Queries");
    emit setProgressRange(0, repeat);

    for(int i = 0; i < repeat; i++)
    {
        emit setProgressValue(i);
        if(canceled)
            return;

        std::vector<QVector2D> inputPoints = randomPoints2D_Uniform(size);
        auto s = std::make_shared<GenericBSpline<QVector2D>>(inputPoints, 3);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
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

