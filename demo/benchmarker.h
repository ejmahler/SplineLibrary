#ifndef BENCHMARKER_H
#define BENCHMARKER_H

#include <QObject>
#include <QMap>
#include <QTime>

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

    void cubicBSplineQuery(int repeat, int queries,  size_t size);
    void genericBSplineQuery(int repeat, int queries, size_t size);

private://support stuff

    QVector2D randomPoint2D(float maxComponent);
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
    bool canceled;
};

#endif // BENCHMARKER_H
