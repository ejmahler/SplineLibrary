#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <memory>
#include <vector>

#include <QWidget>

class Spline;
class Vector3D;

namespace Ui {class MainWindow;};

class GraphicsController;
class Graph;
class FpsCalculator;
class SettingsWidget;
class SplineInverter;

class MainWindow : public QWidget
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();

private slots:
	void settingChanged(void);

private: //methods
	
	void keyPressEvent(QKeyEvent *event);

	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseDoubleClickEvent(QMouseEvent *event);
	
	void rebuildSpline(std::vector<Vector3D> pointList);
	void redraw(void);

	void addVertex(void);
	void deleteVertex(void);

	void createDistanceField(void);


private: //data
	Ui::MainWindow *ui;

	std::shared_ptr<Spline> spline;
	std::shared_ptr<SplineInverter> splineInverter;

	SettingsWidget *settingsWidget;

	GraphicsController *graphicsController;

	bool leftMousePressed;
	bool rightMousePressed;

	int draggedObject;
	int selectedObject;

	QPointF originalMousePos;
};

#endif // MAINWINDOW_H
