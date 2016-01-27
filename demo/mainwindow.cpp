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
#include <QProgressDialog>
#include <QMessageBox>
#include <QtConcurrent>
#include <QFutureWatcher>

#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>

#include <QVector2D>

#include "graphicscontroller.h"
#include "settingswidget.h"
#include "benchmarker.h"

#include "spline_library/hermite/quintic/quintic_hermite_spline.h"
#include "spline_library/hermite/quintic/looping_quintic_hermite_spline.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/looping_cubic_hermite_spline.h"
#include "spline_library/basis/uniform_cubic_bspline.h"
#include "spline_library/basis/looping_uniform_cubic_bspline.h"
#include "spline_library/basis/generic_b_spline.h"
#include "spline_library/basis/looping_generic_b_spline.h"
#include "spline_library/natural/natural_spline.h"
#include "spline_library/natural/looping_natural_spline.h"
#include "spline_library/splineinverter.h"

MainWindow::MainWindow(QWidget *parent)
	: QWidget(parent),
    settingsWidget(new SettingsWidget(this)),
    graphicsController(new GraphicsController(this)),
    benchmarker(nullptr),

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
    std::vector<QVector2D> points;
    points.push_back(QVector2D(100,100));
    points.push_back(QVector2D(400,100));
    points.push_back(QVector2D(500,400));
    points.push_back(QVector2D(300,600));
    points.push_back(QVector2D(300,300));
    points.push_back(QVector2D(150,200));
    points.push_back(QVector2D(100,400));

	rebuildSpline(points);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::settingChanged(void)
{
	//rebuild the spline
    rebuildSpline(mainSpline->getOriginalPoints());
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
    case Qt::Key_B:
        runBenchmark();
        break;
    case Qt::Key_I:
		addVertex();
		break;
    case Qt::Key_D:
        if(mainSpline->getOriginalPoints().size() > 3)
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
            rebuildSpline(mainSpline->getOriginalPoints());
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
            rebuildSpline(mainSpline->getOriginalPoints());
		}
	}

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
    rebuildSpline(mainSpline->getOriginalPoints());
}

void MainWindow::mouseMoveEvent(QMouseEvent *event)
{
	if(leftMousePressed)
	{
		//get the mouse position
		QPoint mousePos = graphicsController->mapFromParent(event->pos());
        QVector2D realPos = graphicsController->convertPoint(mousePos);

		if(draggedObject >= 0)
		{
			//change the point's position
            std::vector<QVector2D> points = mainSpline->getOriginalPoints();
            points[draggedObject] = realPos;

			//rebuild the spline
			rebuildSpline(points);
		}
		else 
		{
            //if we're not currently dragging an object, get the closest T to the mouse position
            auto closest = splineInverter->findClosestT(realPos);

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
    Q_UNUSED(event)
}


void MainWindow::rebuildSpline(std::vector<QVector2D> pointList)
{
    QString mainSplineType = settingsWidget->getOption("main_splineType").toString();
    bool mainIsLooping = settingsWidget->getOption("main_isLooping").toBool();
    float mainAlpha = settingsWidget->getOption("main_alpha").toFloat() / 10;

    bool enableSecondary = settingsWidget->getOption("secondary_enable").toBool();
    QString secondarySplineType = settingsWidget->getOption("secondary_splineType").toString();
    bool secondaryIsLooping = settingsWidget->getOption("secondary_isLooping").toBool();
    float secondaryAlpha = settingsWidget->getOption("secondary_alpha").toFloat() / 10;

    bool includeEndpoints = settingsWidget->getOption("naturalSpline_includeEndpoints").toBool();
    int genericBsplineDegree = settingsWidget->getOption("bSpline_degree").toInt();

    mainSpline = createSpline(pointList, mainSplineType, mainIsLooping, mainAlpha, includeEndpoints,genericBsplineDegree);
    graphicsController->setMainSpline(mainSpline);

    if(enableSecondary)
    {
        secondarySpline = createSpline(pointList, secondarySplineType, secondaryIsLooping, secondaryAlpha, includeEndpoints, genericBsplineDegree);
    }
    else
    {
        secondarySpline = nullptr;

    }
    graphicsController->setSecondarySpline(secondarySpline);

    splineInverter = std::make_shared<SplineInverter<QVector2D>>(*mainSpline.get(), 10);

	DisplayData d;
    d.showConnectingLines = settingsWidget->getOption("misc_showConnectingLines").toBool();
	d.draggedObject = draggedObject;
	d.selectedObject = selectedObject;
	d.imagePath = settingsWidget->getOption("misc_backgroundImagePath").toString();

	graphicsController->draw(d);
}

std::shared_ptr<Spline<QVector2D>> MainWindow::createSpline(
        const std::vector<QVector2D> &pointList,
        const QString &splineType,
        bool isLooping,
        float alpha,
        bool includeEndpoints,
        int genericDegree
        )
{
    if(splineType == "Cubic Catmull-Rom Spline")
    {
        if(isLooping)
        {
            return std::make_shared<LoopingCubicHermiteSpline<QVector2D>>(pointList, alpha);
        }
        else
        {
            return std::make_shared<CubicHermiteSpline<QVector2D>>(pointList, alpha);
        }
    }
    else if(splineType == "Cubic B-Spline")
    {
        if(isLooping)
        {
            return std::make_shared<LoopingUniformCubicBSpline<QVector2D>>(pointList);
        }
        else
        {
            return std::make_shared<UniformCubicBSpline<QVector2D>>(pointList);
        }
    }
    else if(splineType == "Generic B-Spline")
    {
        if(isLooping)
        {
            return std::make_shared<LoopingGenericBSpline<QVector2D>>(pointList, genericDegree);
        }
        else
        {
            return std::make_shared<GenericBSpline<QVector2D>>(pointList, genericDegree);
        }
    }
    else if(splineType == "Cubic Natural Spline")
    {
        if(isLooping)
        {
            return std::make_shared<LoopingNaturalSpline<QVector2D>>(pointList, alpha);
        }
        else
        {
            return std::make_shared<NaturalSpline<QVector2D>>(pointList, includeEndpoints, alpha, NaturalSpline<QVector2D>::Natural);
        }
    }
    else if(splineType == "Cubic Natural Spline (\"Not-A-Knot\")")
    {
        if(isLooping)
        {
            return std::make_shared<LoopingNaturalSpline<QVector2D>>(pointList, alpha);
        }
        else
        {
            return std::make_shared<NaturalSpline<QVector2D>>(pointList, includeEndpoints, alpha, NaturalSpline<QVector2D>::NotAKnot);
        }
    }
    else
    {
        if(isLooping)
        {
            return std::make_shared<LoopingQuinticHermiteSpline<QVector2D>>(pointList, alpha);
        }
        else
        {
            return std::make_shared<QuinticHermiteSpline<QVector2D>>(pointList, alpha);
        }
    }
}

void MainWindow::addVertex(void)
{
    //get the midway point between the chosen vertex and the previous vertex
    //if this is the first vertex, use the next one instead

    std::vector<QVector2D> points = mainSpline->getOriginalPoints();
    if(selectedObject == 0)
    {
        int index = 1;
        float currentT = mainSpline->getT(selectedObject);
        float nextT = mainSpline->getT(index);
        QVector2D halfPoint = mainSpline->getPosition((currentT + nextT) * 0.5);

        //insert this new point
        points.insert(points.begin() + index, halfPoint);
        selectedObject = index;
    }
    else
    {
        int index = selectedObject - 1;
        float currentT = mainSpline->getT(selectedObject);
        float nextT = mainSpline->getT(index);
        QVector2D halfPoint = mainSpline->getPosition((currentT + nextT) * 0.5);

        //insert this new point
        points.insert(points.begin() + selectedObject, halfPoint);
    }

	//redraw spline
	rebuildSpline(points);
}

void MainWindow::deleteVertex(void)
{
    std::vector<QVector2D> points = mainSpline->getOriginalPoints();

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
        graphicsController->clearBackground();
        settingsWidget->setOption("misc_backgroundImagePath", saveFileName);
	}
}

void MainWindow::runBenchmark(void)
{
    //if the benchmarker is non-null, there's already one running
    if(benchmarker != nullptr)
        return;

    benchmarker = new Benchmarker(this);

    //use a progress dialog to make sur ethe user knows how long it's going to take
    QProgressDialog dialog(this, Qt::Dialog);
    dialog.setAutoReset(false);
    connect(benchmarker, &Benchmarker::setProgressText, &dialog, &QProgressDialog::setLabelText);
    connect(benchmarker, &Benchmarker::setProgressRange, &dialog, &QProgressDialog::setRange);
    connect(benchmarker, &Benchmarker::setProgressValue, &dialog, &QProgressDialog::setValue);
    connect(&dialog, &QProgressDialog::canceled, benchmarker, &Benchmarker::cancel);

    //begin the task
    auto future = QtConcurrent::run(benchmarker, &Benchmarker::runBenchmark);

    //use a futurewatcher to know when the task is complete. when it is, close the progress dialog
    QFutureWatcher<void> watcher;
    watcher.setFuture(future);
    connect(&watcher, &QFutureWatcher<void>::finished, &dialog, &QProgressDialog::reset);

    //wait for the task to complete
    dialog.exec();

    //if the task was canceled, make sure we wait for the thread to clean up before returning
    if(dialog.wasCanceled())
    {
        future.waitForFinished();
    }
    else
    {
        //build the string that we're going to display to the user
        QString resultText;
        auto result = future.result();
        for(auto it = result.cbegin(); it != result.cend(); it++)
        {
            resultText += QString("%1: %2s\n").arg(it.key(), QString::number(it.value(), 'f', 2));
        }


        //create a messagebox to display the results
        QMessageBox resultBox(this);
        resultBox.setText("Benchmark Results");
        resultBox.setInformativeText(resultText);
        resultBox.setWindowModality(Qt::WindowModal);

        //wait for the user to acknowledge the box
        resultBox.exec();
    }

    //clean up
    delete benchmarker;
    benchmarker = nullptr;
}
