#pragma once

#include <QObject>

class TestCalculus : public QObject
{
    Q_OBJECT
public:
    explicit TestCalculus(QObject *parent = 0);

signals:

private slots:
    void testGaussLegendre_data(void);
    void testGaussLegendre(void);
};
