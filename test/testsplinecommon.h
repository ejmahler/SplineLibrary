#pragma once

#include <QObject>

class TestSplineCommon : public QObject
{
    Q_OBJECT
public:
    explicit TestSplineCommon(QObject *parent = nullptr);

signals:

private slots:
    //test the computeTValuesWithInnerPadding method, which computes the T values for a non-looping spline
    void testInnerPadding_data(void);
    void testInnerPadding(void);
};
