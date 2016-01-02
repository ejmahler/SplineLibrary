#ifndef BENCHMARKER_H
#define BENCHMARKER_H

#include <QObject>
#include <QMap>
#include <QTime>

#include <vector>

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
    void naturalSplineQueryBalanced(int repeat, int queries, size_t size, float alpha);
    void naturalSplineQueryUnbalanced(int repeat, int queries, size_t size, float alpha);

private://support stuff

    std::vector<QVector2D> randomPoints2D_Uniform(size_t size, float maxComponent);
    std::vector<QVector2D> randomPoints2D_Unbalanced(size_t size, float maxComponent);
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
