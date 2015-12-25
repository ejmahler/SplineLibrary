#ifndef GRAPHICSCONTROLLER_H
#define GRAPHICSCONTROLLER_H

#include <memory>

#include <QGLWidget>

#include "spline_library/spline.h"

class QVector2D;

struct DisplayData
{
    bool showConnectingLines;
	int selectedObject;
	int draggedObject;
	bool highlightT;
    float highlightedT;
	QString imagePath;

	DisplayData(void)
        :showConnectingLines(false), selectedObject(0), draggedObject(-1), highlightT(false), highlightedT(0),
         imagePath()
	{}
};

class Graph;
class QPainter;

class GraphicsController : public QGLWidget
{
	Q_OBJECT

public:
	GraphicsController(QWidget *parent);
	~GraphicsController();

    void setMainSpline(const std::shared_ptr<Spline<QVector2D>> &s);
    void setSecondarySpline(const std::shared_ptr<Spline<QVector2D>> &s);

	void draw(const DisplayData &d);

    QPoint convertPoint(const QVector2D &point);
    QVector2D convertPoint(const QPoint &point);

    void clearBackground(void);

	void createDistanceField(const QString &filename);

	int pickVertex(const QPoint &screenPoint);

protected:
     void paintEvent(QPaintEvent *event) override;

private:

     void drawSpline(QPainter &painter, const std::shared_ptr<Spline<QVector2D>> &s, const QColor &color);
     void drawSplineSegment(
             QPainter &painter,
             const std::shared_ptr<Spline<QVector2D>> &s,
             float beginT,
             float endT);

     void drawSplineDerivative(QPainter &painter, const std::shared_ptr<Spline<QVector2D>> &s, const QColor &color);
     void drawSplineSegmentDerivative(
             QPainter &painter,
             const std::shared_ptr<Spline<QVector2D>> &s,
             float beginT,
             float endT);

    void drawPoints(QPainter &painter, const std::vector<QVector2D> &points);

    void drawDiagnosticText(QPainter &painter, int top,
            const QString &labelText,const QString &valueText);
	void drawControlText(QPainter &painter, int top, 
		const QString &labelText,const QString &valueText);

    void keyPressEvent(QKeyEvent *event) override;

    QVector3D getColor(float t) const;

private:

	QMap<QString,QStaticText> staticText;

    std::shared_ptr<Spline<QVector2D>> mainSpline;
    std::shared_ptr<Spline<QVector2D>> secondarySpline;

	DisplayData displayData;

	bool displayControls;

    int pointRadius;

	std::shared_ptr<QImage> backgroundImage;
	QString backgroundImagePath;

    std::shared_ptr<Spline<QVector3D>> colorSpline;
};

#endif // GRAPHICSCONTROLLER_H
