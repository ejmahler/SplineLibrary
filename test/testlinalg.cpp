#include "testlinalg.h"

#include "spline_library/utils/linearalgebra.h"

#include <vector>

#include <QtTest/QtTest>
#include <QDebug>

TestLinAlg::TestLinAlg(QObject *parent) : QObject(parent)
{

}

void TestLinAlg::testTridiagonal_data(void)
{
    QTest::addColumn<std::vector<float>>("lower_diagonal");
    QTest::addColumn<std::vector<float>>("main_diagonal");
    QTest::addColumn<std::vector<float>>("upper_diagonal");
    QTest::addColumn<std::vector<float>>("input");
    QTest::addColumn<std::vector<float>>("expected_output");

    QTest::newRow("Identity") << std::vector<float> { 0.0f, 0.0f } << std::vector<float> { 1.0f, 1.0f, 1.0f } << std::vector<float> { 0.0f, 0.0f } << std::vector<float> { 5.0f, 5.0f, 5.0f } << std::vector<float> { 5.0f, 5.0f, 5.0f };
    QTest::newRow("Simple")
            << std::vector<float> { 1.0f, 1.0f, 1.0f }
            << std::vector<float> { 2.0f, 2.0f, 2.0f, 2.0f }
            << std::vector<float> { 1.0f, 1.0f, 1.0f }
            << std::vector<float> { 5.0f, 5.0f, 5.0f, 5.0f }
            << std::vector<float> { 2.0f, 1.0f, 1.0f, 2.0f };

    QTest::newRow("AllDifferent")
            << std::vector<float> { 1.0f, 2.0f, 3.0f }
            << std::vector<float> { 11.0f, 12.0f, 13.0f, 14.0f }
            << std::vector<float> { 4.0f, 5.0f, 6.0f }
            << std::vector<float> { 7.0f, 8.0f, 9.0f, 10.0f }
            << std::vector<float> { 0.455891f, 0.496299f, 0.317705f, 0.646206f };
}
void TestLinAlg::testTridiagonal(void)
{
    QFETCH(std::vector<float>, lower_diagonal);
    QFETCH(std::vector<float>, main_diagonal);
    QFETCH(std::vector<float>, upper_diagonal);
    QFETCH(std::vector<float>, input);
    QFETCH(std::vector<float>, expected_output);

    auto result = LinearAlgebra::solveTridiagonal(main_diagonal, upper_diagonal, lower_diagonal, input);

    for(size_t i = 0; i < result.size(); i++) {
        QCOMPARE(result[i], expected_output[i]);
    }
}

void TestLinAlg::testSymmetricTridiagonal_data(void)
{
    QTest::addColumn<std::vector<float>>("main_diagonal");
    QTest::addColumn<std::vector<float>>("secondary_diagonal");
    QTest::addColumn<std::vector<float>>("input");

    QTest::newRow("Identity") << std::vector<float> { 1.0f, 1.0f, 1.0f } << std::vector<float> { 0.0f, 0.0f } << std::vector<float> { 5.0f, 5.0f, 5.0f };
    QTest::newRow("Simple")
            << std::vector<float> { 2.0f, 2.0f, 2.0f, 2.0f }
            << std::vector<float> { 1.0f, 1.0f, 1.0f }
            << std::vector<float> { 5.0f, 5.0f, 5.0f, 5.0f };

    QTest::newRow("AllDifferent")
            << std::vector<float> { 11.0f, 12.0f, 13.0f, 14.0f }
            << std::vector<float> { 1.0f, 2.0f, 3.0f }
            << std::vector<float> { 7.0f, 8.0f, 9.0f, 10.0f };
}
void TestLinAlg::testSymmetricTridiagonal(void)
{
    QFETCH(std::vector<float>, main_diagonal);
    QFETCH(std::vector<float>, secondary_diagonal);
    QFETCH(std::vector<float>, input);

    //verify that the symmetric algorithm gives the same result as the non-symmetric one
    auto result = LinearAlgebra::solveSymmetricTridiagonal(main_diagonal, secondary_diagonal, input);
    auto expected = LinearAlgebra::solveTridiagonal(main_diagonal, secondary_diagonal, secondary_diagonal, input);

    for(size_t i = 0; i < result.size(); i++) {
        QCOMPARE(result[i], expected[i]);
    }
}

void TestLinAlg::testCyclicTridiagonal_data(void)
{
    QTest::addColumn<std::vector<float>>("main_diagonal");
    QTest::addColumn<std::vector<float>>("secondary_diagonal");
    QTest::addColumn<std::vector<float>>("input");
    QTest::addColumn<std::vector<float>>("expected_output");

    QTest::newRow("Identity") << std::vector<float> { 1.0f, 1.0f, 1.0f } << std::vector<float> { 0.0f, 0.0f, 0.0f } << std::vector<float> { 5.0f, 5.0f, 5.0f } << std::vector<float> { 5.0f, 5.0f, 5.0f };
    QTest::newRow("Simple")
            << std::vector<float> { 3.0f, 3.0f, 3.0f, 3.0f, 3.0f }
            << std::vector<float> { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f }
            << std::vector<float> { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f }
            << std::vector<float> { 0.2f, 0.2f, 0.2f, 0.2f, 0.2f };

    QTest::newRow("AllDifferent")
            << std::vector<float> { 10.0f, 11.0f, 12.0f, 13.0f, 14.0f }
            << std::vector<float> { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f }
            << std::vector<float> { 20.0f, 21.0f, 22.0f, 23.0f, 24.0f }
            << std::vector<float> { 1.41308f, 1.54923f, 1.27271f, 1.20969f, 0.863988f };
}
void TestLinAlg::testCyclicTridiagonal(void)
{
    QFETCH(std::vector<float>, main_diagonal);
    QFETCH(std::vector<float>, secondary_diagonal);
    QFETCH(std::vector<float>, input);
    QFETCH(std::vector<float>, expected_output);

    auto result = LinearAlgebra::solveCyclicSymmetricTridiagonal(main_diagonal, secondary_diagonal, input);

    for(size_t i = 0; i < result.size(); i++) {
        QCOMPARE(result[i], expected_output[i]);
    }
}
