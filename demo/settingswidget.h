#ifndef SETTINGSWIDGET_H
#define SETTINGSWIDGET_H

#include <QDialog>
#include <QMap>

namespace Ui {class SettingsWidget;}

class SettingsWidget : public QDialog
{
	Q_OBJECT

public:
	SettingsWidget(QWidget *parent = 0);
	~SettingsWidget();

    QVariant getOption(const QString &option) const;
    void setOption(const QString &option, const QVariant &value);

signals:
	void settingChanged(void);

private slots:
	void on_lineEdit_changed(const QString &text);
	void on_checkBox_changed(bool checked);
	void on_spinBox_changed(int value);
	void on_doubleSpinBox_changed(double value);
	void on_slider_changed(int value);

	void on_buttonChooseFile_clicked(void);

private:
	QString getSettingCategory(const QString &objectName);
	QString getSettingName(const QString &objectName);
	void loadSettings(void);
	void setSignals(void);

private:
	Ui::SettingsWidget *ui;

	QMap<QString,QVariant> options;
};

#endif // SETTINGSWIDGET_H
