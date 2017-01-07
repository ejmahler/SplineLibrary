#pragma once

#include <QObject>

class TestLinAlg : public QObject
{
    Q_OBJECT
public:
    explicit TestLinAlg(QObject *parent = 0);

signals:

private slots:
    void testTridiagonal_data(void);
    void testTridiagonal(void);

    void testSymmetricTridiagonal_data(void);
    void testSymmetricTridiagonal(void);

    void testCyclicTridiagonal_data(void);
    void testCyclicTridiagonal(void);
};
