#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <memory>
#include <vector>

#include <QWidget>

#include "spline_library/spline.h"
#include "spline_library/splineinverter.h"

#include <QVector2D>

namespace Ui {class MainWindow;}

class GraphicsController;
class Graph;
class FpsCalculator;
class SettingsWidget;
class Benchmarker;


class MainWindow : public QWidget
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();

private slots:
	void settingChanged(void);

private: //methods
	
    void keyPressEvent(QKeyEvent *event) override;

    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
	
    void rebuildSpline(std::vector<QVector2D> pointList);
    std::shared_ptr<Spline<QVector2D>> createSpline(
            const std::vector<QVector2D> &pointList,
            const QString &splineType,
            bool looping,
            float alpha,
            bool includeEndpoints,
            int bSplineDegree
            );

	void redraw(void);

	void addVertex(void);
    void deleteVertex(void);

	void createDistanceField(void);

    void runBenchmark(void);


private: //data
	Ui::MainWindow *ui;

    std::shared_ptr<Spline<QVector2D>> mainSpline;
    std::shared_ptr<Spline<QVector2D>> secondarySpline;
    std::shared_ptr<SplineInverter<QVector2D>> splineInverter;

	SettingsWidget *settingsWidget;

	GraphicsController *graphicsController;
    Benchmarker *benchmarker;

	bool leftMousePressed;
	bool rightMousePressed;

	int draggedObject;
	int selectedObject;

    QPoint originalMousePos;
};

template<>
std::array<float, 2> convertPoint<QVector2D, float, 2>(const QVector2D& point);

#endif // MAINWINDOW_H
