#include "benchmarker.h"

#include <vector>
#include <memory>

#include <QVector2D>
#include <QTime>

#include "spline_library/basis/cubic_b_spline.h"
#include "spline_library/basis/generic_b_spline.h"

Benchmarker::Benchmarker(QObject *parent) : QObject(parent)
{

}

QMap<QString, float> Benchmarker::runBenchmark(void)
{
    canceled = false;

    QMap<QString, float> results;

    results["Cubic B-Spline Queries"] = timeFunction(&Benchmarker::cubicBSplineQuery, 10, 100000, 100000);
    results["Generic B-Spline Queries"] = timeFunction(&Benchmarker::genericBSplineQuery, 10, 100000, 100000);

    return results;
}

void Benchmarker::cancel(void)
{
    canceled = true;
}

void Benchmarker::cubicBSplineQuery(int repeat, int queries, size_t size)
{
    qsrand(10);

    emit setProgressText("Running Cubic B-Spline Queries");
    emit setProgressRange(0, repeat * queries);

    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<CubicBSpline<QVector2D>>(inputPoints);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            if(q % 100 == 0)
            {
                emit setProgressValue(i * queries + q);
                if(canceled)
                    return;
            }

            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat * queries);
}

void Benchmarker::genericBSplineQuery(int repeat, int queries, size_t size)
{
    emit setProgressText("Running Generic B-Spline Queries");
    emit setProgressRange(0, repeat * queries);

    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<GenericBSpline<QVector2D>>(inputPoints, 3);

        float max = s->getMaxT();
        for(int q = 0; q < queries; q++)
        {
            if(q % 100 == 0)
            {
                emit setProgressValue(i * queries + q);
                if(canceled)
                    return;
            }

            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
    emit setProgressValue(repeat * queries);
}

QVector2D Benchmarker::randomPoint2D(float maxComponent)
{
    float multiplier = maxComponent / RAND_MAX;
    return QVector2D(qrand() * multiplier, qrand() * multiplier);
}

float Benchmarker::randomFloat(float max)
{

}

