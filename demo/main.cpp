#include "mainwindow.h"
#include <QApplication>

#include <QGLFormat>

int main(int argc, char *argv[])
{
	//set up antialiasing for opengl rendering
	QGLFormat glf = QGLFormat::defaultFormat(); 
	glf.setSampleBuffers(true); 
	glf.setSamples(8); 
	QGLFormat::setDefaultFormat(glf);

	QApplication a(argc, argv);
	MainWindow w;
	w.show();
	return a.exec();
}
