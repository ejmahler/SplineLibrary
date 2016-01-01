#include "benchmarker.h"

#include <vector>

#include <QVector2D>

#include "spline_library/basis/cubic_b_spline.h"
#include "spline_library/basis/generic_b_spline.h"

Benchmarker::Benchmarker(QObject *parent) : QObject(parent)
{

}

void Benchmarker::cubicBSplineConstructor(int repeat, size_t size)
{
    qsrand(10);

    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<CubicBSpline<QVector2D>>(inputPoints);
        s->getMaxT();
    }
}

void Benchmarker::genericBSplineConstructor(int repeat, size_t size)
{
    qsrand(10);

    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<GenericBSpline<QVector2D>>(inputPoints, 3);
        s->getMaxT();
    }
}

void Benchmarker::cubicBSplineQuery(int repeat, int queries,  size_t size)
{
    qsrand(10);

    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<CubicBSpline<QVector2D>>(inputPoints);

        float max = s->getMaxT();
        for(int i = 0; i < queries; i++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
}

void Benchmarker::genericBSplineQuery(int repeat, int queries, size_t size)
{
    for(int i = 0; i < repeat; i++)
    {
        std::vector<QVector2D> inputPoints(size);
        for(size_t s = 0; s < size; s++)
        {
            inputPoints[s] = randomPoint2D(1000);
        }
        auto s = std::make_shared<GenericBSpline<QVector2D>>(inputPoints, 3);

        float max = s->getMaxT();
        for(int i = 0; i < queries; i++)
        {
            float t = randomFloat(max);
            s->getPosition(t);
        }
    }
}

QVector2D Benchmarker::randomPoint2D(float maxComponent)
{
    float multiplier = maxComponent / RAND_MAX;
    return QVector2D(qrand() * multiplier, qrand() * multiplier);
}

float Benchmarker::randomFloat(float max)
{

}

