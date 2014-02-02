#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <ctime>
#include <cmath>

#include <QGLWidget>
#include <QVector>
#include <QKeyEvent>
#include <QFile>
#include <QDir>
#include <QFileDialog>

#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>

#include "graphicscontroller.h"
#include "settingswidget.h"

#include "spline_library/quintic_hermite/quintic_cr_spline.h"
#include "spline_library/quintic_hermite/looping_quintic_cr_spline.h"
#include "spline_library/cubic_hermite/cr_spline.h"
#include "spline_library/cubic_hermite/looping_cr_spline.h"
#include "spline_library/b_spline/cubic_b_spline.h"
#include "spline_library/b_spline/looping_cubic_b_spline.h"
#include "spline_library/natural_spline/natural_spline.h"
#include "spline_library/natural_spline/looping_natural_spline.h"
#include "spline_library/splineinverter.h"

MainWindow::MainWindow(QWidget *parent)
	: QWidget(parent),
	settingsWidget(new SettingsWidget()),
	graphicsController(new GraphicsController(this)),

	leftMousePressed(0),
	rightMousePressed(0),
	draggedObject(-1),
	selectedObject(0)
{
	ui = new Ui::MainWindow();
	ui->setupUi(this);

	ui->containerWidget->layout()->addWidget(graphicsController);

	connect(
		settingsWidget,
		SIGNAL(settingChanged(void)),
		this,
		SLOT(settingChanged(void))
		);

	qsrand(time(0));
	std::vector<Vector3D> points;
	points.push_back(Vector3D(100,100,0));
	points.push_back(Vector3D(400,100,0));
	points.push_back(Vector3D(500,400,0));
	points.push_back(Vector3D(300,600,0));
	points.push_back(Vector3D(300,300,0));
    points.push_back(Vector3D(100,200,0));

	rebuildSpline(points);
}

MainWindow::~MainWindow()
{
	delete ui;
	delete settingsWidget;
}

void MainWindow::settingChanged(void)
{
	//rebuild the spline
	rebuildSpline(spline->getPoints());
}


void MainWindow::keyPressEvent(QKeyEvent *event)
{
	switch(event->key())
	{
	case Qt::Key_S:
		settingsWidget->show();
		settingsWidget->raise();
		break;
    case Qt::Key_G:
		createDistanceField();
		break;
    case Qt::Key_I:
		addVertex();
		break;
    case Qt::Key_D:
		if(spline->getPoints().size() > 3)
			deleteVertex();
		break;
	default:
		event->ignore();
	}
}

void MainWindow::mousePressEvent(QMouseEvent *event)
{
	if(event->button() == Qt::LeftButton)
	{
		leftMousePressed = true;

		//get the mouse position
		QPoint mousePos = graphicsController->mapFromParent(event->pos());

		//see if an object is underneath the mouse
		int objectId = graphicsController->pickVertex(mousePos);

		if(objectId >= 0)
		{
			draggedObject = objectId;
		}
	}
	else if(event->button() == Qt::RightButton)
	{
		rightMousePressed = true;

		//get the mouse position
		QPoint mousePos = graphicsController->mapFromParent(event->pos());

		//see if an object is underneath the mouse
		int objectId = graphicsController->pickVertex(mousePos);

		if(objectId >= 0)
		{
			selectedObject = objectId;
		}
	}
	rebuildSpline(spline->getPoints());
}

void MainWindow::mouseReleaseEvent(QMouseEvent *event)
{
	if(event->button() == Qt::LeftButton)
	{
		leftMousePressed = false;
		draggedObject = -1;
	}
	else if(event->button() == Qt::RightButton)
	{
		rightMousePressed = false;
	}
	rebuildSpline(spline->getPoints());
}

void MainWindow::mouseMoveEvent(QMouseEvent *event)
{
	if(leftMousePressed)
	{
		//get the mouse position
		QPoint mousePos = graphicsController->mapFromParent(event->pos());
		Vector3D realPos = graphicsController->convertPoint(mousePos);

		if(draggedObject >= 0)
		{
			//change the point's position
			std::vector<Vector3D> points = spline->getPoints();
            points[draggedObject] = realPos;

			//rebuild the spline
			rebuildSpline(points);
		}
		else 
		{
			//if we're not currently dragging an object, get the current T to the mouse position
			float closest = splineInverter->findClosestFast(realPos);

			//redraw to highlight this new closest T
			DisplayData d;
            d.showConnectingLines = settingsWidget->getOption("misc_showConnectingLines").toBool();
			d.draggedObject = draggedObject;
			d.selectedObject = selectedObject;

			d.highlightT = true;
			d.highlightedT = closest;

			d.imagePath = settingsWidget->getOption("misc_backgroundImagePath").toString();

			graphicsController->draw(d);
		}
	}
}

void MainWindow::mouseDoubleClickEvent(QMouseEvent *event)
{

}


void MainWindow::rebuildSpline(std::vector<Vector3D> pointList)
{
    QString splineType = settingsWidget->getOption("main_splineType").toString();
    bool isLooping = settingsWidget->getOption("splineType_isLooping").toBool();
    double alpha = settingsWidget->getOption("cubicHermite_alpha").toDouble() / 10;

    if(splineType == "Cubic Catmull-Rom Spline")
    {
        if(isLooping)
        {
            spline = std::shared_ptr<Spline>(new LoopingCRSpline(pointList, alpha));
        }
        else
        {
            spline = std::shared_ptr<Spline>(new CRSpline(pointList, alpha));
        }
    }
    else if(splineType == "Cubic B-Spline")
    {
        if(isLooping)
        {
            spline = std::shared_ptr<Spline>(new LoopingCubicBSpline(pointList));
        }
        else
        {
            spline = std::shared_ptr<Spline>(new CubicBSpline(pointList));
        }
    }
    else if(splineType == "Cubic Natural Spline")
    {
        if(isLooping)
        {
            spline = std::shared_ptr<Spline>(new LoopingNaturalSpline(pointList, alpha));
        }
        else
        {
            bool includeEndpoints = settingsWidget->getOption("naturalSpline_includeEndpoints").toBool();
            spline = std::shared_ptr<Spline>(new NaturalSpline(pointList, includeEndpoints, alpha));
        }
    }
    else
    {
        if(isLooping)
        {
            spline = std::shared_ptr<Spline>(new LoopingQuinticCRSpline(pointList, alpha));
        }
        else
        {
            spline = std::shared_ptr<Spline>(new QuinticCRSpline(pointList, alpha));
        }
    }

	splineInverter = std::shared_ptr<SplineInverter>(new SplineInverter(spline, 100));

	graphicsController->setSpline(spline);

	DisplayData d;
    d.showConnectingLines = settingsWidget->getOption("misc_showConnectingLines").toBool();
	d.draggedObject = draggedObject;
	d.selectedObject = selectedObject;
	d.imagePath = settingsWidget->getOption("misc_backgroundImagePath").toString();

	graphicsController->draw(d);
}

void MainWindow::addVertex(void)
{
    //get the midway point between the chosen vertex and the previous vertex
    //if this is the first vertex, use the next one instead

    std::vector<Vector3D> points = spline->getPoints();
    if(selectedObject == 0)
    {
        int index = 1;
        float currentT = spline->getT(selectedObject);
        float nextT = spline->getT(index);
        Vector3D halfPoint = spline->getPosition((currentT + nextT) * 0.5);

        //insert this new point
        points.insert(points.begin() + index, halfPoint);
        selectedObject = index;
    }
    else
    {
        int index = selectedObject - 1;
        float currentT = spline->getT(selectedObject);
        float nextT = spline->getT(index);
        Vector3D halfPoint = spline->getPosition((currentT + nextT) * 0.5);

        //insert this new point
        points.insert(points.begin() + selectedObject, halfPoint);
    }

	//redraw spline
	rebuildSpline(points);
}

void MainWindow::deleteVertex(void)
{
	std::vector<Vector3D> points = spline->getPoints();

	//remove the object
	points.erase(points.begin() + selectedObject);

    //return the "previous" index. if the current
    if(selectedObject > points.size())
        selectedObject = points.size() - 1;

	//redraw spline
	rebuildSpline(points);
}

void MainWindow::createDistanceField(void)
{
	QString caption = "Save Distance Field Image";

	QDir fileDirectory = QDir::home();
	fileDirectory.cd("Pictures");
	QString dir = fileDirectory.absoluteFilePath("distanceField.png");
	QString filter = "PNG Image (*.png)";

	QString saveFileName = QFileDialog::getSaveFileName(this, caption, dir, filter);

	if(saveFileName.length() > 0)
	{
		graphicsController->createDistanceField(saveFileName);
	}
}
