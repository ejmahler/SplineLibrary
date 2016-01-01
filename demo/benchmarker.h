#ifndef BENCHMARKER_H
#define BENCHMARKER_H

#include <QObject>

class Benchmarker : public QObject
{
    Q_OBJECT
public:
    explicit Benchmarker(QObject *parent = 0);

signals:

public slots:
    //**********
    //all of these functions can change based on whatever you want - i just needed a common place to put performance comparisons

    void cubicBSplineConstructor(int repeat, size_t size);
    void genericBSplineConstructor(int repeat, size_t size);

    void cubicBSplineQuery(int repeat, int queries,  size_t size);
    void genericBSplineQuery(int repeat, int queries, size_t size);

private:
    QVector2D randomPoint2D(float maxComponent);
    float randomFloat(float max);
};

#endif // BENCHMARKER_H
