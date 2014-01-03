#ifndef GRAPHICSCONTROLLER_H
#define GRAPHICSCONTROLLER_H

#include <memory>

#include <QGLWidget>

class Vector3D;

class Spline;

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

class GraphicsController : public QGLWidget
{
	Q_OBJECT

public:
	GraphicsController(QWidget *parent);
	~GraphicsController();

    void setSpline(const std::shared_ptr<Spline> &s);
	void draw(const DisplayData &d);

	QPoint convertPoint(const Vector3D &point);
	Vector3D convertPoint(const QPoint &point);

	void createDistanceField(const QString &filename);

	int pickVertex(const QPoint &screenPoint);

protected:
	 void paintEvent(QPaintEvent *event);

private:

	void drawControlText(QPainter &painter, int top, 
		const QString &labelText,const QString &valueText);

	void keyPressEvent(QKeyEvent *event);

    Vector3D getColor(float t) const;

private:

	QMap<QString,QStaticText> staticText;
	std::shared_ptr<Spline> spline;

	DisplayData displayData;

	bool displayControls;

	int pointRadius;

	std::shared_ptr<QImage> backgroundImage;
	QString backgroundImagePath;

	std::shared_ptr<Spline> colorSpline;
};

#endif // GRAPHICSCONTROLLER_H
