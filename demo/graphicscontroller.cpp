#include "graphicscontroller.h"

#include <QGLWidget>
#include <QStaticText>

#include <QKeyEvent>

#include <QVector2D>
#include <QVector3D>

#include <ctime>
#include <algorithm>

#include "spline_library/splineinverter.h"
#include "spline_library/splinelengthcalculator.h"
#include "spline_library/hermite/cubic/cubic_hermite_spline.h"
#include "spline_library/hermite/cubic/looping_cubic_hermite_spline.h"

GraphicsController::GraphicsController(QWidget *parent)
	: QGLWidget(parent), 
	
    mainSpline(nullptr),
	displayControls(false),
	pointRadius(10),
	backgroundImagePath()
{
	setFocusPolicy(Qt::StrongFocus);
}

GraphicsController::~GraphicsController()
{

}

void GraphicsController::setMainSpline(const std::shared_ptr<Spline<QVector2D>> &s)
{
    mainSpline = s;
}

void GraphicsController::setSecondarySpline(const std::shared_ptr<Spline<QVector2D>> &s)
{
    secondarySpline = s;
}

void GraphicsController::draw(const DisplayData &d)
{
	displayData = d;
	update();
}

void GraphicsController::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event)

	//load image if necessary
	if(displayData.imagePath != backgroundImagePath)
	{
		backgroundImagePath = displayData.imagePath;
        backgroundImage = std::make_shared<QImage>(backgroundImagePath);

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
	
    //draw points
    drawPoints(painter, mainSpline->getOriginalPoints());

	//draw the highlighted point if it exists
	if(displayData.highlightT)
	{
		painter.save();

        QVector2D result = mainSpline->getPosition(displayData.highlightedT);

		painter.translate(result.x(),result.y());

		QColor color = Qt::yellow;
		painter.setPen(color);
		painter.setBrush(color);

		float size = pointRadius;
		QRectF pieRect(-size/2,-size/2,size,size);
		painter.drawPie(pieRect,0,5760);

		painter.restore();
	}

    //draw the spline itself
    drawSpline(painter, mainSpline, Qt::red);
    if(secondarySpline != nullptr)
    {
        drawSpline(painter, secondarySpline, Qt::blue);
    }
	
	painter.restore();

    SplineLengthCalculator<QVector2D> lengthCalc(mainSpline);

	//draw container for diagnostic data
	painter.setOpacity(0.75);
	QColor bgColor(qRgb(32,32,32));
	painter.setPen(bgColor);
    painter.setBrush(bgColor);

    int diagonsticBoxWidth = 200;
    int diagnosticBoxHeight = 100;
    painter.drawRect(0, 0, diagonsticBoxWidth, diagnosticBoxHeight);

	//draw text for diagnostic data
	painter.setOpacity(1);
    painter.setPen(Qt::white);

    drawDiagnosticText(painter, 5, "Spline Length", QString::number(lengthCalc.findLength(0, mainSpline->getMaxT())));

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


QPoint GraphicsController::convertPoint(const QVector2D &point)
{
	return QPoint(int(point.x()), int(point.y()));
}

QVector2D GraphicsController::convertPoint(const QPoint &point)
{
    return QVector2D(point.x(), point.y());
}

void GraphicsController::createDistanceField(const QString &filename)
{
	QImage output(700,700, QImage::Format_ARGB32_Premultiplied);
	QPainter painter(&output);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setRenderHint(QPainter::HighQualityAntialiasing);
	painter.setRenderHint(QPainter::SmoothPixmapTransform);

	qsrand(time(0));

    std::vector<QVector3D> colorList;
    for(size_t i = 0; i < mainSpline->getOriginalPoints().size(); i++)
	{
        colorList.push_back(QVector3D(
            float(qrand()) / RAND_MAX,
            float(qrand()) / RAND_MAX,
            float(qrand()) / RAND_MAX));
	}
    if(mainSpline->isLooping())
        colorSpline = std::make_shared<LoopingCubicHermiteSpline<QVector3D>>(colorList);
    else
        colorSpline = std::make_shared<CubicHermiteSpline<QVector3D>>(colorList);

	painter.fillRect(0,0,output.width(),output.height(),Qt::white);

    SplineInverter<QVector2D> calc(mainSpline, 40);

	//supersampling amount - 1 is no supersampling
    int supersampling = 4;
	int totalSamples = supersampling * supersampling;
    float base = 1 / (float(supersampling) * 2);
    float step = 1 / float(supersampling);

	//for every pixel in the image, determine the closest T value
	for(int y = 0; y < output.height(); y++)
	{
		for(int x = 0; x < output.width(); x++)
		{
			//use 2x supersampling, with a simple grid
            QVector3D colorVector;
			for(int dy = 0; dy < supersampling; dy++) {
				for(int dx = 0; dx < supersampling; dx++) {
                    QVector2D realPoint(
						x + base + dx * step,
                        y + base + dy * step);
                    float tfast = calc.findClosestT(realPoint);
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

    //draw points
    drawPoints(painter, mainSpline->getOriginalPoints());

	//draw lines
    painter.setPen(Qt::red);

    float stepSize = 1.0 / 100;
    float currentStep = stepSize;
    float limit = mainSpline->getMaxT();
    QVector2D previousPoint = mainSpline->getPosition(0);

	painter.setPen(Qt::red);

	while(currentStep <= limit)
	{
        QVector2D currentPoint = mainSpline->getPosition(currentStep);

		painter.drawLine(
			QPointF(previousPoint.x(),previousPoint.y()),
			QPointF(currentPoint.x(),currentPoint.y())
			);


		currentStep += stepSize;
		previousPoint = currentPoint;
    }

	output.save(filename);
}

int GraphicsController::pickVertex(const QPoint &screenPoint)
{
    std::vector<QVector2D> points = mainSpline->getOriginalPoints();

	//convert this point to world coordinates and then just loop through every vertex
    QVector2D convertedPoint = convertPoint(screenPoint);
    for(size_t i = 0; i < points.size(); i++)
	{
        QVector2D position = points[i];

		if((position - convertedPoint).lengthSquared() < pointRadius * pointRadius)
		{
			return i;
		}
	}

	return -1;
}

void GraphicsController::drawDiagnosticText(QPainter &painter, int top,
    const QString &labelText,const QString &valueText)
{
    int labelLeft = 5;
    int labelWidth = 140;
    int valueLeft = labelLeft + labelWidth + 5;
    int valueWidth = 60;

    if(!staticText.contains(labelText))
    {
        QStaticText t(labelText);
        t.setPerformanceHint(QStaticText::AggressiveCaching);
        t.setTextFormat(Qt::PlainText);
        t.setTextWidth(labelWidth);
        t.setTextOption(QTextOption(Qt::AlignRight));

        staticText.insert(labelText,t);
    }

    painter.drawStaticText(labelLeft,top,staticText.value(labelText));

    //don't use static text for the actual numbers, since they change every frame
    painter.drawText(QRect(valueLeft, top, valueWidth, 20), valueText, QTextOption(Qt::AlignRight));
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

QVector3D GraphicsController::getColor(float t) const
{
    QVector3D value = colorSpline->getPosition(t);

    return QVector3D(
        qBound(0.0f,value.x() * 255.0f, 255.0f),
        qBound(0.0f,value.y() * 255.0f, 255.0f),
        qBound(0.0f,value.z() * 255.0f, 255.0f));
}

void GraphicsController::drawSpline(QPainter &painter, const std::shared_ptr<Spline<QVector2D>> &s, const QColor &color)
{
    //draw the spline
    float stepSize = 0.25;
    float currentStep = stepSize;
    float limit = s->getMaxT() + 0.01;

    float thresholdAngle = 0.996;

    painter.setPen(color);

    while(currentStep <= limit)
    {
        drawSplineSegment(painter, s,currentStep - stepSize, currentStep, thresholdAngle);
        currentStep += stepSize;
    }
}

void GraphicsController::drawSplineSegment(QPainter &painter, const std::shared_ptr<Spline<QVector2D>> &s, float beginT, float endT, float thresholdAngle)
{
    auto beginData = s->getCurvature(beginT);
    auto endData = s->getCurvature(endT);

    QVector2D beginNormalizedTangent = beginData.tangent.normalized();
    QVector2D endNormalizedTangent = endData.tangent.normalized();

    //compute the angle between the two tangents
    auto cosAngle = QVector2D::dotProduct(beginNormalizedTangent, endNormalizedTangent);

    //if the angle is too low, subdivide this segment into two segments
    auto minDelta = .001;
    if(cosAngle < thresholdAngle && (endT - beginT) > minDelta)
    {
        //we dont want the exact middle, give a bias to the end whose curvature rejected against the tangent is longer
        //if one side's curvature rejection is larger, it means that side is turning faster, so we can reduce the number of segments
        //by making the segment division closer to that side
        //not completely necessary or practical here but it shows off a potential use of the tangent and curvature
        auto beginProjection = QVector2D::dotProduct(beginNormalizedTangent, beginData.curvature);
        auto endProjection = QVector2D::dotProduct(endNormalizedTangent, endData.curvature);

        QVector2D beginRejection = beginData.curvature - beginNormalizedTangent * beginProjection;
        QVector2D endRejection = endData.curvature - endNormalizedTangent * endProjection;

        auto beginRejectionLength = beginRejection.length();
        auto endRejectionLength = endRejection.length();

        auto lengthSum = (beginRejectionLength + endRejectionLength) * 0.05;

        auto beginPercent = (beginRejectionLength + lengthSum) / (lengthSum * 2 + beginRejectionLength + endRejectionLength);

        auto middle = beginT + (endT - beginT) * (1 - beginPercent);
        drawSplineSegment(painter, s, beginT, middle, thresholdAngle);
        drawSplineSegment(painter, s, middle, endT, thresholdAngle);
    }
    else
    {
        painter.drawLine(
            QPointF(beginData.position.x(),beginData.position.y()),
            QPointF(endData.position.x(),endData.position.y())
            );
    }
}

void GraphicsController::drawPoints(QPainter &painter, const std::vector<QVector2D> &points)
{
    int size = points.size();

    //draw straight lines connecting each control point
    painter.setPen(qRgb(32,32,32));

    if(displayData.showConnectingLines)
    {
        for(size_t i = 0; i < points.size() - 1; i++)
        {
            painter.drawLine(
                        QPointF(points.at(i).x(),points.at(i).y()),
                        QPointF(points.at(i + 1).x(),points.at(i + 1).y())
                        );
        }

        if(mainSpline->isLooping())
        {
            painter.drawLine(
                        QPointF(points.at(0).x(),points.at(0).y()),
                        QPointF(points.at(points.size() - 1).x(),points.at(points.size() - 1).y())
                        );
        }
    }

    //draw control points on top of lines
    for(int i = 0; i < size; i++)
    {
        painter.save();

        QVector2D position = points.at(i);

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
}

