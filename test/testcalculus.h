#pragma once

#include <QObject>

class TestCalculus : public QObject
{
    Q_OBJECT
public:
    explicit TestCalculus(QObject *parent = 0);

signals:

private slots:
    void testAdaptiveSimpsons_data(void);
    void testAdaptiveSimpsons(void);
};
