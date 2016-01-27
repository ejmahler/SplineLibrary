#pragma once

#include <QObject>

class TestVector : public QObject
{
    Q_OBJECT
public:
    explicit TestVector(QObject *parent = 0);

signals:

private slots:
    void testConstructors(void);

    void testVectorArithmetic_data(void);
    void testVectorArithmetic(void);

    void testScalarArithmetic_data(void);
    void testScalarArithmetic(void);

    void testLengthOperations_data(void);
    void testLengthOperations(void);

    //verify that we can create a spline using Vector as the interpolation type and get the expected results
    void testSplineFunctionality_data(void);
    void testSplineFunctionality(void);
};
