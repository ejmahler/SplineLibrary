#pragma once

#include <QObject>
#include <QMap>
#include <QTime>

#include <vector>
#include <random>
#include <functional>

#include "spline_library/vector.h"

typedef Vector<2, double> VectorT;

class Benchmarker : public QObject
{
    Q_OBJECT
public:
    explicit Benchmarker(QObject *parent = 0);

signals:
    void setProgressText(const QString & text);
    void setProgressRange(int minimum, int maximum);
    void setProgressValue(int progress);

public slots:
    QMap<QString, float> runBenchmark(void);

    void cancel(void);

private:
    //**********
    //all of these functions can change based on whatever you want - i just needed a common place to put performance comparisons

    void testPrecision(QString message, int repeat, int queries, size_t size, std::function<std::vector<VectorT> (size_t)> pointGenerator);

private://support stuff

    std::vector<VectorT> randomPoints2D_Uniform(size_t size);
    std::vector<VectorT> randomPoints2D_SmallVariance(size_t size);
    double randomFloat(double max);

    template <class Function, typename... Args>
    float timeFunction(Function f, Args&&... a) {
        QTime t;
        t.start();
        (this->*f)(std::forward<Args>(a)...);
        int elapsedMilli = t.elapsed();
        return float(elapsedMilli) / 1000.0;
    }

private: //data
    std::minstd_rand gen;
    std::uniform_real_distribution<double> bigDistribution;
    std::uniform_real_distribution<double> smallDistribution;
    std::uniform_real_distribution<double> smallVarianceDistribution;
    bool canceled;
};
