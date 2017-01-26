#pragma once

#include <QObject>
#include <QMap>
#include <QTime>

#include <vector>
#include <random>
#include <functional>
#include <memory>

#include "spline_library/vector.h"
#include "spline_library/spline.h"

const size_t D = 2;
typedef float FloatingT;
typedef Vector<D, FloatingT> VectorT;
typedef LoopingSpline<VectorT, FloatingT> SplineType;

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
    void timeSplineMemberFunction(
            QMap<QString, float>& results,
            void(Benchmarker::*testFunction)(int, const SplineType&),
            std::function<std::unique_ptr<SplineType>(size_t)> splineFunction,
            QString message, int queries, size_t size);

    //**********
    //all of these functions can change based on whatever you want - i just needed a common place to put performance comparisons
    void testArcLength(int queries, const SplineType &spline);

private://support stuff

    template<class InterpolationType, size_t dimension, typename floating_t, class RandomSource>
    static InterpolationType makeRandomPoint(RandomSource randomSource)
    {
        std::array<floating_t, dimension> result{};

        for(size_t i = 0; i < dimension; i++) {
            result[i] = randomSource();
        }

        return result;
    }

    template<class InterpolationType, size_t dimension, typename floating_t, class RandomSource>
    static std::vector<InterpolationType> randomPoints_Uniform(RandomSource randomSource, size_t size)
    {
        std::vector<InterpolationType> result(size);
        for(size_t i = 0; i < size; i++)
        {
            result[i] = makeRandomPoint<InterpolationType, dimension, floating_t>(randomSource);
        }
        return result;
    }

    template<class InterpolationType, size_t dimension, typename floating_t, class RandomSource>
    static std::vector<InterpolationType> randomPoints_SmallVariance(RandomSource randomSource, size_t size)
    {
        std::vector<InterpolationType> result(size);
        result[0] = makeRandomPoint<InterpolationType, dimension, floating_t>(randomSource);
        for(size_t i = 1; i < size; i++)
        {
            result[i] = result[i - 1] + makeRandomPoint<InterpolationType, dimension, floating_t>(randomSource);
        }
        return result;
    }

private: //data
    std::minstd_rand gen;
    bool canceled;
    int repeats;
};
