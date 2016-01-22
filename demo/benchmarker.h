#ifndef BENCHMARKER_H
#define BENCHMARKER_H

#include <QObject>
#include <QMap>
#include <QTime>

#include <vector>
#include <random>

#include "spline_library/vector3d.h"

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

    void recursiveLength(int repeat, int queries,  size_t size);
    void builtinLength(int repeat, int queries,  size_t size);

private://support stuff

    std::vector<QVector2D> randomPoints2D_Uniform(size_t size);
    std::vector<QVector2D> randomPoints2D_SmallVariance(size_t size);
    float randomFloat(float max);

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
    std::uniform_real_distribution<float> bigDistribution;
    std::uniform_real_distribution<float> smallDistribution;
    std::uniform_real_distribution<float> smallVarianceDistribution;
    bool canceled;
};

#endif // BENCHMARKER_H
