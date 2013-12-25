#include "graphicscontroller.h"

#include <QGLWidget>
#include <QStaticText>

#include <QKeyEvent>

#include <ctime>
#include <algorithm>

#include "spline_library/vector3d.h"
#include "spline_library/splineinverter.h"
#include "spline_library/cubic_hermite/cr_spline.h"
#include "spline_library/cubic_hermite/looping_cr_spline.h"

GraphicsController::GraphicsController(QWidget *parent)
	: QGLWidget(parent), 
	
	spline(nullptr),
	displayControls(false),
	pointRadius(10),
	backgroundImagePath()
{
	setFocusPolicy(Qt::StrongFocus);
}

GraphicsController::~GraphicsController()
{

}

void GraphicsController::setSpline(const std::shared_ptr<Spline> &s)
{
	spline = s;
}

void GraphicsController::draw(const DisplayData &d)
{
	displayData = d;
	update();
}

void GraphicsController::paintEvent(QPaintEvent *event)
{
	//load image if necessary
	if(displayData.imagePath != backgroundImagePath)
	{
		backgroundImagePath = displayData.imagePath;
		backgroundImage = std::shared_ptr<QImage>(new QImage(backgroundImagePath));

		//if loading failed, set the image to null
		if(backgroundImage->isNull())
		{
			backgroundImage = nullptr;
		}
	}


	QPainter painter;

	makeCurrent();
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_LINE_SMOOTH);

	int widgetWidth = width();
	int widgetHeight = height();

	painter.begin(this);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setRenderHint(QPainter::HighQualityAntialiasing);
	painter.setRenderHint(QPainter::SmoothPixmapTransform);

	//paint background
	painter.fillRect(0,0,widgetWidth,widgetHeight,Qt::black);

	painter.save();
	
	if(backgroundImage != nullptr)
	{
		painter.drawImage(0,0,*backgroundImage);
    }
	
    //draw the spline
    double stepSize = 0.01;
    double currentStep = stepSize;
    double limit = spline->getMaxT();
    Vector3D previousPoint = spline->getPosition(0);

    painter.setPen(Qt::blue);

    while(currentStep <= limit)
    {
        auto currentData = spline->getPosition(currentStep);

        painter.drawLine(
            QPointF(previousPoint.x(),previousPoint.y()),
            QPointF(currentData.x(),currentData.y())
            );

        currentStep += stepSize;
        previousPoint = currentData;
    }

    std::vector<Vector3D> points = spline->getPoints();

    //draw control points on top of line
    for(int i = 0; i < points.size(); i++)
	{
		painter.save();

		Vector3D position = points.at(i);

		painter.translate(position.x(),position.y());

		QColor color;
		if(displayData.draggedObject == i)
		{
			color = QColor(qRgb(255,128,64));
		}
		else if(displayData.selectedObject == i)
		{
			color = Qt::cyan;
		}
		else
		{
			color = Qt::white;
		}
		painter.setPen(color);
		painter.setBrush(color);

		float size = pointRadius;
		QRectF pieRect(-size/2,-size/2,size,size);
		painter.drawPie(pieRect,0,5760);

		painter.restore();
	}

	//draw the highlighted point if it exists
	if(displayData.highlightT)
	{
		painter.save();

		Vector3D result = spline->getPosition(displayData.highlightedT);

		painter.translate(result.x(),result.y());

		QColor color = Qt::yellow;
		painter.setPen(color);
		painter.setBrush(color);

		float size = pointRadius;
		QRectF pieRect(-size/2,-size/2,size,size);
		painter.drawPie(pieRect,0,5760);

		painter.restore();
	}

	
	painter.restore();





	//draw container for diagnostic data
	painter.setOpacity(0.75);
	QColor bgColor(qRgb(32,32,32));
	painter.setPen(bgColor);
	painter.setBrush(bgColor);
	//painter.drawRect(0,0,250,120);

	//draw text for diagnostic data
	painter.setOpacity(1);
	painter.setPen(Qt::white);
	//drawDiagnosticText(painter, 5, "Graphics FPS:", QString::number(d.graphicsFps,'f',2));


	//draw container for control data
    int controlBoxWidth = 225;
	int controlBoxHeight = 150;
	if(!displayControls)
		controlBoxHeight = 25;


	painter.setOpacity(0.75);
	painter.setPen(bgColor);
	painter.setBrush(bgColor);
	painter.drawRect(width() - controlBoxWidth, 0, controlBoxWidth, controlBoxHeight);


	//draw text for control data
	painter.setOpacity(1);
	painter.setPen(Qt::white);
	
	if(displayControls)
	{
		drawControlText(painter,5,"[c]", "Hide Controls");
		drawControlText(painter,25,"[insert]", "Add Vertex");
        drawControlText(painter,45,"[delete]", "Delete Vertex");
        drawControlText(painter,65,"[i]", "Generate distance field");
        drawControlText(painter,85,"[s]", "Open Settings");
	}
	else
	{
		drawControlText(painter,5,"[c]", "Show Controls");
	}

	painter.end();
}


QPoint GraphicsController::convertPoint(const Vector3D &point)
{
	return QPoint(int(point.x()), int(point.y()));
}

Vector3D GraphicsController::convertPoint(const QPoint &point)
{
	return Vector3D(point.x(), point.y(),0);
}

void GraphicsController::createDistanceField(const QString &filename)
{
	QImage output(700,700, QImage::Format_ARGB32_Premultiplied);
	QPainter painter(&output);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setRenderHint(QPainter::HighQualityAntialiasing);
	painter.setRenderHint(QPainter::SmoothPixmapTransform);

	qsrand(time(0));

	std::vector<Vector3D> colorList;
	for(int i = 0; i < spline->getPoints().size(); i++)
	{
		colorList.push_back(Vector3D(
			double(qrand()) / RAND_MAX,
			double(qrand()) / RAND_MAX,
			double(qrand()) / RAND_MAX));
	}
    if(spline->isLooping())
        colorSpline = std::shared_ptr<Spline>(new LoopingCRSpline(colorList));
    else
        colorSpline = std::shared_ptr<Spline>(new CRSpline(colorList));

	painter.fillRect(0,0,output.width(),output.height(),Qt::white);

	SplineInverter calc(spline, 40);

	//supersampling amount - 1 is no supersampling
    int supersampling = 4;
	int totalSamples = supersampling * supersampling;
	double base = 1 / (double(supersampling) * 2);
	double step = 1 / double(supersampling);

	//for every pixel in the image, determine the closest T value
	for(int y = 0; y < output.height(); y++)
	{
		for(int x = 0; x < output.width(); x++)
		{
			//use 2x supersampling, with a simple grid
			Vector3D colorVector;
			for(int dy = 0; dy < supersampling; dy++) {
				for(int dx = 0; dx < supersampling; dx++) {
					Vector3D realPoint(
						x + base + dx * step,
						y + base + dy * step, 
						0);
                    double tfast = calc.findClosestPrecise(realPoint);
                    colorVector += getColor(tfast);
				}
			}

			output.setPixel(x,y,
				qRgb(
					colorVector.x() / totalSamples,
					colorVector.y() / totalSamples,
					colorVector.z() / totalSamples
				)
			);
		}
	}

	//draw lines
	painter.setPen(Qt::red);
	std::vector<Vector3D> points = spline->getPoints();

	double stepSize = 1.0 / 100;
	double currentStep = stepSize;
	double limit = spline->getMaxT();
	Vector3D previousPoint = spline->getPosition(0);

	painter.setPen(Qt::red);

	while(currentStep <= limit)
	{
		Vector3D currentPoint = spline->getPosition(currentStep);

		painter.drawLine(
			QPointF(previousPoint.x(),previousPoint.y()),
			QPointF(currentPoint.x(),currentPoint.y())
			);


		currentStep += stepSize;
		previousPoint = currentPoint;
	}

	//draw control points on top of line
	//if this is a loop, don't draw the final point
	int pointLimit;
    if(spline->isLooping())
		pointLimit = points.size() - 1;
	else
		pointLimit = points.size();
	for(int i = 0; i < pointLimit; i++)
	{
		painter.save();

		Vector3D position = points.at(i);

		painter.translate(position.x(),position.y());

		QColor color = Qt::white;
		painter.setPen(color);
		painter.setBrush(color);

		float size = pointRadius/2 + 1;
		QRectF pieRect(-size/2,-size/2,size,size);
		painter.drawPie(pieRect,0,5760);

		painter.restore();
	}

	output.save(filename);
}

int GraphicsController::pickVertex(const QPoint &screenPoint)
{
	std::vector<Vector3D> points = spline->getPoints();

	//convert this point to world coordinates and then just loop through every vertex
	Vector3D convertedPoint = convertPoint(screenPoint);
	for(int i = 0; i < points.size(); i++)
	{
		Vector3D position = points[i];

		if((position - convertedPoint).lengthSquared() < pointRadius * pointRadius)
		{
			return i;
		}
	}

	return -1;
}

void GraphicsController::drawControlText(QPainter &painter, int top, 
	const QString &labelText,const QString &valueText)
{
    int labelLeft = width() - 220 + 5;
	int labelWidth = 40;
	int valueLeft = labelLeft + labelWidth + 10;
	int valueWidth = 150;

	if(!staticText.contains(labelText))
	{
		QStaticText t(labelText);
		t.setPerformanceHint(QStaticText::AggressiveCaching);
		t.setTextFormat(Qt::PlainText);
		t.setTextWidth(labelWidth);
		t.setTextOption(QTextOption(Qt::AlignRight));

		staticText.insert(labelText,t);
	}

	if(!staticText.contains(valueText))
	{
		QStaticText t(valueText);
		t.setPerformanceHint(QStaticText::AggressiveCaching);
		t.setTextFormat(Qt::PlainText);
		t.setTextWidth(valueWidth);
		t.setTextOption(QTextOption(Qt::AlignLeft));

		staticText.insert(valueText,t);
	}

	painter.drawStaticText(labelLeft,top,staticText.value(labelText));
	painter.drawStaticText(valueLeft,top,staticText.value(valueText));
}

void GraphicsController::keyPressEvent(QKeyEvent *event)
{
	switch(event->key())
	{
	case Qt::Key_C:
		displayControls = !displayControls;
		update();
		break;
	default:
		event->ignore();
	}
}

Vector3D GraphicsController::getColor(float t) const
{
	Vector3D value = colorSpline->getPosition(t);

	return Vector3D(
		qBound(0.0,value.x() * 255,255.0),
		qBound(0.0,value.y() * 255,255.0),
		qBound(0.0,value.z() * 255,255.0));
}

