#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <memory>
#include <vector>

#include <QWidget>

#include "spline_library/spline.h"
#include "spline_library/splineinverter.h"

class Vector3D;

namespace Ui {class MainWindow;}

class GraphicsController;
class Graph;
class FpsCalculator;
class SettingsWidget;


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
	
	void rebuildSpline(std::vector<Vector3D> pointList);
    std::shared_ptr<Spline<Vector3D>> createSpline(const std::vector<Vector3D> &pointList, const QString &splineType, bool looping, double alpha, bool includeEndpoints);

	void redraw(void);

	void addVertex(void);
	void deleteVertex(void);

	void createDistanceField(void);


private: //data
	Ui::MainWindow *ui;

    std::shared_ptr<Spline<Vector3D>> mainSpline;
    std::shared_ptr<Spline<Vector3D>> secondarySpline;
    std::shared_ptr<SplineInverter<Vector3D>> splineInverter;

	SettingsWidget *settingsWidget;

	GraphicsController *graphicsController;

	bool leftMousePressed;
	bool rightMousePressed;

	int draggedObject;
	int selectedObject;

	QPointF originalMousePos;
};

#endif // MAINWINDOW_H
