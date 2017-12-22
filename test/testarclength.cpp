#include "testarclength.h"

#include "common.h"
#include "spline_library/utils/arclength.h"

#include "spline_library/utils/calculus.h"
#include "spline_library/splines/uniform_cubic_bspline.h"
#include "spline_library/splines/generic_b_spline.h"
#include "spline_library/splines/natural_spline.h"
#include "spline_library/splines/cubic_hermite_spline.h"
#include "spline_library/splines/uniform_cr_spline.h"
#include "spline_library/splines/quintic_hermite_spline.h"

#include <QtTest/QtTest>

TestArcLength::TestArcLength(QObject *parent) : QObject(parent)
{

}

void TestArcLength::testArcLengthTotalLength_data(void)
{
    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");

    auto data = TestDataFloat::generateRandomData(10);

    QTest::newRow("uniformCR") << TestDataFloat::createUniformCR(data);
    QTest::newRow("CubicHermiteAlpha") << TestDataFloat::createCubicHermite(data, 0.5f);
}

void TestArcLength::testArcLengthTotalLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);

    float arc = spline->arcLength(0, spline->getMaxT());
    float total = spline->totalLength();

    QCOMPARE(arc, total);
}


void TestArcLength::testKnownArcLength_data(void)
{
    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<float>("a");
    QTest::addColumn<float>("b");
    QTest::addColumn<float>("expectedLength");

    auto data = TestDataFloat::generateTriangleNumberData(10);

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {

        //add a test row for the whole spline
        std::string allName = QString("%1 (All)").arg(name).toStdString();
        QTest::newRow(allName.data()) << spline << 0.0f << spline->getMaxT() << (data[data.size() - 1] - data[0]).length();

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(2), spline->segmentT(3), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (Partial)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB << (spline->getPosition(partialA) - spline->getPosition(partialB)).length();

        //add a row where a and b are in the same segment, since most splines treat this as a special case
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (Same)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB << (spline->getPosition(sameSegmentA) - spline->getPosition(sameSegmentB)).length();
    };

    rowFunction("uniformCR", TestDataFloat::createUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createCubicHermite(data, 0.5f));
}

void TestArcLength::testKnownArcLength(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);
    QFETCH(float, expectedLength);

    float arcResult = spline->arcLength(a, b);
    QCOMPARE(arcResult, expectedLength);
}



void TestArcLength::testCyclicArcLength_data(void)
{
    QTest::addColumn<std::shared_ptr<LoopingSpline<Vector2>>>("spline");
    QTest::addColumn<float>("a");
    QTest::addColumn<float>("b");

    auto data = TestDataFloat::generateRandomData(10);

    auto rowFunction = [=](const char* name, std::shared_ptr<LoopingSpline<Vector2>> spline) {

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(2), spline->segmentT(3), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (DifferentSegment)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB;

        //add a row where a and b are in the same segment
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (SameSegment)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB;
    };

    rowFunction("uniformCR", TestDataFloat::createLoopingUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createLoopingCatmullRom(data, 0.5f));
}

void TestArcLength::testCyclicArcLength(void)
{
    QFETCH(std::shared_ptr<LoopingSpline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);

    float maxT = spline->getMaxT();

    //when a < b and both values are in range, the result should be the same as arcLength
    float arcResult = spline->arcLength(a, b);
    float cyclicArcResult = spline->cyclicArcLength(a, b);
    QCOMPARE(cyclicArcResult, arcResult);

    float totalLength = spline->totalLength();

    //verify that if we reverse a and b, we get totalLength - arcResult
    float reversedResult = spline->cyclicArcLength(b, a);
    QCOMPARE(reversedResult, totalLength - arcResult);

    //verify that we can independently move a and b out of range and get the same result
    float outOfRangeA = spline->cyclicArcLength(a + maxT, b);
    float outOfRangeB = spline->cyclicArcLength(a, b + maxT);
    float outOfRangeAReversed = spline->cyclicArcLength(b, a + maxT);
    float outOfRangeBReversed = spline->cyclicArcLength(b + maxT, a);

    QCOMPARE(outOfRangeA, arcResult);
    QCOMPARE(outOfRangeB, arcResult);
    QCOMPARE(outOfRangeAReversed, reversedResult);
    QCOMPARE(outOfRangeBReversed, reversedResult);
}

void TestArcLength::testSolve_data(void)
{
    auto data = TestDataFloat::generateRandomData(10);

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<float>("a");
    QTest::addColumn<float>("b");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(1), spline->segmentT(2), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (Partial)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB;

        //add a row where a and b are in the same segment, since this is a special case
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (Same)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB;
    };

    rowFunction("uniformCR", TestDataFloat::createUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createCubicHermite(data, 0.5f));
}

void TestArcLength::testSolve(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);

    float arcLength = spline->arcLength(a, b);

    float calculatedB = ArcLength::solveLength(*spline.get(), a, arcLength);
    QCOMPARE(calculatedB, b);

    //verify that if the desiredLength is longer than the spline, maxT is returned
    float totalLength = spline->totalLength();
    float calculatedOverLength = ArcLength::solveLength(*spline.get(), a, totalLength);
    QCOMPARE(calculatedOverLength, spline->getMaxT());
}


void TestArcLength::testSolveCyclic_data(void)
{
    auto data = TestDataFloat::generateRandomData(10);

    QTest::addColumn<std::shared_ptr<LoopingSpline<Vector2>>>("spline");
    QTest::addColumn<float>("a");
    QTest::addColumn<float>("b");

    auto rowFunction = [=](const char* name, std::shared_ptr<LoopingSpline<Vector2>> spline) {

        //add a row for just part of the spline. we want to make sure a and b fall partway through a segment
        //so we'll explicitly get segment boundaries via spline->getT and lerp between them
        float partialA = lerp(spline->segmentT(1), spline->segmentT(2), 0.75f);
        float partialB = lerp(spline->segmentT(spline->segmentCount() - 3), spline->segmentT(spline->segmentCount() - 2), 0.25f);
        std::string partialName = QString("%1 (Partial)").arg(name).toStdString();
        QTest::newRow(partialName.data()) << spline << partialA << partialB;

        //add a row where a and b are in the same segment, since this is a special case
        size_t testIndex = 3;
        float sameSegmentA = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.2f);
        float sameSegmentB = lerp(spline->segmentT(testIndex), spline->segmentT(testIndex + 1), 0.6f);
        std::string sameSegmentName = QString("%1 (Same)").arg(name).toStdString();
        QTest::newRow(sameSegmentName.data()) << spline << sameSegmentA << sameSegmentB;
    };

    rowFunction("uniformCR", TestDataFloat::createLoopingUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createLoopingCubicHermite(data, 0.5f));
}

void TestArcLength::testSolveCyclic(void)
{
    QFETCH(std::shared_ptr<LoopingSpline<Vector2>>, spline);
    QFETCH(float, a);
    QFETCH(float, b);

    float maxT = spline->getMaxT();
    float desiredLength = spline->arcLength(a, b);
    float totalLength = spline->totalLength();

    //when a and what will eventually be b are both in range, solve and solveCyclic should return the same result
    float calculatedB = ArcLength::solveLength(*spline, a, desiredLength);
    float cyclicB = ArcLength::solveLengthCyclic(*spline, a, desiredLength);
    QCOMPARE(cyclicB, calculatedB);

    //when desiredLength is greater than total length, we should go multiple times around the spline
    float cyclicBCycle1 = ArcLength::solveLengthCyclic(*spline, a, desiredLength + totalLength);
    float cyclicBCycle2 = ArcLength::solveLengthCyclic(*spline, a, desiredLength + totalLength*2);
    QCOMPARE(cyclicBCycle1, calculatedB + maxT);
    QCOMPARE(cyclicBCycle2, calculatedB + 2*maxT);

    //when A isn't in range, the result should respect the unwrapped A, rather than the wrapped A
    float cyclicBCycle3 = ArcLength::solveLengthCyclic(*spline, a + maxT, desiredLength + totalLength*2);
    float cyclicBCycleN = ArcLength::solveLengthCyclic(*spline, a - maxT, desiredLength);
    QCOMPARE(cyclicBCycle3, calculatedB + 3*maxT);
    QCOMPARE(cyclicBCycleN, calculatedB - maxT);

    //test that it works when we give an input value that has to wrap around but do less than one full cycle
    //the result should be on the "next cycle" relative to the input, as opposed to being wrapped around to be less than the input
    // The easiest way to do this is to just reverse the inputs, so that we go from b to a
    float reversedLength = totalLength - desiredLength;
    float calculatedA = ArcLength::solveLengthCyclic(*spline, b, reversedLength);
    float calculatedNegA = ArcLength::solveLengthCyclic(*spline, b - maxT, reversedLength);
    QCOMPARE(calculatedA, a + maxT);
    QCOMPARE(calculatedNegA, a);
}




void TestArcLength::testPartition_data(void)
{
    auto data = TestDataFloat::generateRandomData(10);

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<float>("desiredLength");
    QTest::addColumn<size_t>("expectedPieces");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {
        float totalLength = spline->totalLength();

        //add a row where the desired length is larger than the average segment
        std::string largeName = QString("%1 (Large Pieces)").arg(name).toStdString();
        float largeLength = totalLength / 2.1f;
        QTest::newRow(largeName.data()) << spline << largeLength << size_t(2);

        //add a row where the desired length is low, so there will be many results that begin and end in the same segment, since this is a special case
        std::string smallName = QString("%1 (Small Pieces)").arg(name).toStdString();
        float smallLength = totalLength / 20.5f;
        QTest::newRow(smallName.data()) << spline << smallLength << size_t(20);
    };

    rowFunction("uniformCR", TestDataFloat::createUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createCubicHermite(data, 0.5f));
}

void TestArcLength::testPartition(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(float, desiredLength);
    QFETCH(size_t, expectedPieces);

    std::vector<float> pieces = ArcLength::partition(*spline.get(), desiredLength);

    //verify that we have the correct number of results
    QCOMPARE(pieces.size(), expectedPieces + 1);

    //verify that each piece has the correct arc length
    for(size_t i = 0; i < expectedPieces; i++)
    {
        float pieceLength = spline->arcLength(pieces[i], pieces[i+1]);
        QCOMPARE(pieceLength, desiredLength);
    }
}





void TestArcLength::testPartitionN_data(void)
{
    auto data = TestDataFloat::generateRandomData(10);

    QTest::addColumn<std::shared_ptr<Spline<Vector2>>>("spline");
    QTest::addColumn<size_t>("n");

    auto rowFunction = [=](const char* name, std::shared_ptr<Spline<Vector2>> spline) {
        //add a row where the desired length is larger than the average segment
        std::string largeName = QString("%1 (Large Pieces)").arg(name).toStdString();
        QTest::newRow(largeName.data()) << spline << size_t(3);

        //add a row where the desired length is low, so there will be many results that begin and end in the same segment, since this is a special case
        std::string smallName = QString("%1 (Small Pieces)").arg(name).toStdString();
        QTest::newRow(smallName.data()) << spline << size_t(20);
    };

    rowFunction("uniformCR", TestDataFloat::createUniformCR(data));
    rowFunction("cubicHermiteAlpha", TestDataFloat::createCubicHermite(data, 0.5f));
}

void TestArcLength::testPartitionN(void)
{
    QFETCH(std::shared_ptr<Spline<Vector2>>, spline);
    QFETCH(size_t, n);

    std::vector<float> pieces = ArcLength::partitionN(*spline.get(), n);

    //verify that we have the correct number of results
    QCOMPARE(pieces.size(), n + 1);

    float totalLength = spline->totalLength();

    //verify that each piece has the correct arc length
    for(size_t i = 0; i < n; i++)
    {
        float pieceLength = spline->arcLength(pieces[i], pieces[i+1]);
        QCOMPARE(pieceLength, totalLength/n);
    }
}
