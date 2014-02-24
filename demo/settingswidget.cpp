#include "settingswidget.h"
#include "ui_settingswidget.h"

#include <QFileDialog>
#include <QRadioButton>
#include <QCheckBox>
#include <QSpinBox>

#include "settings.h"

SettingsWidget::SettingsWidget(QWidget *parent)
    : QDialog(parent)
{
	ui = new Ui::SettingsWidget();
	ui->setupUi(this);

	loadSettings();
	setSignals();
}

SettingsWidget::~SettingsWidget()
{
	delete ui;
}

QVariant SettingsWidget::getOption(const QString &option) const
{
	return options.value(option);
}



void SettingsWidget::on_lineEdit_changed(const QString &text)
{
	QString name = sender()->objectName();
	options[name] = text;

	Settings::saveSetting(getSettingCategory(name),getSettingName(name),text);
	emit settingChanged();
}

void SettingsWidget::on_checkBox_changed(bool checked)
{
	QString name = sender()->objectName();
	options[name] = checked;

	Settings::saveSetting(getSettingCategory(name),getSettingName(name),checked);
	emit settingChanged();
}

void SettingsWidget::on_slider_changed(int value)
{
	QString name = sender()->objectName();
	options[name] = value;

	Settings::saveSetting(getSettingCategory(name),getSettingName(name), value);
	emit settingChanged();
}

void SettingsWidget::on_spinBox_changed(int value)
{
	QString name = sender()->objectName();
	options[name] = value;

	Settings::saveSetting(getSettingCategory(name),getSettingName(name), value);
	emit settingChanged();
}

void SettingsWidget::on_doubleSpinBox_changed(double value)
{
	QString name = sender()->objectName();
	options[name] = value;

	Settings::saveSetting(getSettingCategory(name),getSettingName(name), value);
	emit settingChanged();
}


void SettingsWidget::on_buttonChooseFile_clicked(void)
{
	QString caption = "Choose Background Image";
	QString filter = "Images (*.png *.jpg)";

	QString fileName = QFileDialog::getOpenFileName(this, caption, QString(), filter);
	if(fileName.size() > 0)
	{
		ui->misc_backgroundImagePath->setText(fileName);
	}
}


QString SettingsWidget::getSettingCategory(const QString &objectName)
{
	return objectName.split("_")[0];
}
QString SettingsWidget::getSettingName(const QString &objectName)
{
	return objectName.split("_")[1];
}


void SettingsWidget::loadSettings(void)
{
    //go through all line edits
    for(QLineEdit *x: findChildren<QLineEdit*>())
    {
        QString category = getSettingCategory(x->objectName());
        QString name = getSettingName(x->objectName());

        if(category != "qt" && Settings::hasSetting(category,name))
        {
            x->setText(Settings::getSetting(category,name).toString());
            options[x->objectName()] = Settings::getSetting(category,name);
        }
        else
        {
            options[x->objectName()] = x->text();
        }
    }

    //go through all combo boxes
    for(QComboBox *x: findChildren<QComboBox*>())
    {
        QString category = getSettingCategory(x->objectName());
        QString name = getSettingName(x->objectName());

        if(category != "qt" && Settings::hasSetting(category,name))
        {
            x->setCurrentText(Settings::getSetting(category,name).toString());
            options[x->objectName()] = Settings::getSetting(category,name);
        }
        else
        {
            options[x->objectName()] = x->currentText();
        }
    }

	//go through all checkboxes
    for(QRadioButton *x: findChildren<QRadioButton*>())
	{
		QString category = getSettingCategory(x->objectName());
		QString name = getSettingName(x->objectName());

		if(Settings::hasSetting(category,name))
		{
			x->setChecked(Settings::getSetting(category,name).toBool());
			options[x->objectName()] = Settings::getSetting(category,name);
		}
		else
		{
			options[x->objectName()] = x->isChecked();
		}
	}

	//go through all radio buttons
    for(QCheckBox *x: findChildren<QCheckBox*>())
	{
		QString category = getSettingCategory(x->objectName());
		QString name = getSettingName(x->objectName());

		if(Settings::hasSetting(category,name))
		{
			x->setChecked(Settings::getSetting(category,name).toBool());
			options[x->objectName()] = Settings::getSetting(category,name);
		}
		else
		{
			options[x->objectName()] = x->isChecked();
		}
	}

	//go through all spin boxes
    for(QSlider *x: findChildren<QSlider*>())
	{
		QString category = getSettingCategory(x->objectName());
		QString name = getSettingName(x->objectName());

		if(Settings::hasSetting(category,name))
		{
			x->setValue(Settings::getSetting(category,name).toInt());
			options[x->objectName()] = Settings::getSetting(category,name);
		}
		else
		{
			options[x->objectName()] = x->value();
		}
	}

	//go through all spin boxes
    for(QSpinBox *x: findChildren<QSpinBox*>())
	{
		QString category = getSettingCategory(x->objectName());
		QString name = getSettingName(x->objectName());

		if(Settings::hasSetting(category,name))
		{
			x->setValue(Settings::getSetting(category,name).toInt());
			options[x->objectName()] = Settings::getSetting(category,name);
		}
		else
		{
			options[x->objectName()] = x->value();
		}
	}

	//go through all double spin boxes
    for(QDoubleSpinBox *x: findChildren<QDoubleSpinBox*>())
	{
		QString category = getSettingCategory(x->objectName());
		QString name = getSettingName(x->objectName());

		if(Settings::hasSetting(category,name))
		{
			x->setValue(Settings::getSetting(category,name).toDouble());
			options[x->objectName()] = Settings::getSetting(category,name);
		}
		else
		{
			options[x->objectName()] = x->value();
		}
	}
}

void SettingsWidget::setSignals(void)
{
    //connect all line edits
    for(QLineEdit *x: findChildren<QLineEdit*>())
    {
        connect(
            x,
            SIGNAL(textChanged(const QString&)),
            this,
            SLOT(on_lineEdit_changed(const QString&))
            );
    }

    //connect all combo boxes
    for(QComboBox *x: findChildren<QComboBox*>())
    {
        connect(
            x,
            SIGNAL(currentTextChanged(const QString&)),
            this,
            SLOT(on_lineEdit_changed(const QString&))
            );
    }

	//connect all radio buttons
    for(QRadioButton *x: findChildren<QRadioButton*>())
	{
		connect(
			x,
			SIGNAL(toggled(bool)),
			this,
			SLOT(on_checkBox_changed(bool))
			);
	}

	//connect all checkboxes
    for(QCheckBox *x: findChildren<QCheckBox*>())
	{
		connect(
			x,
			SIGNAL(toggled(bool)),
			this,
			SLOT(on_checkBox_changed(bool))
			);
	}

	//connect all sliders
    for(QSlider *x: findChildren<QSlider*>())
	{
		connect(
			x,
			SIGNAL(valueChanged(int)),
			this,
			SLOT(on_slider_changed(int))
			);
	}

	//connect all spin boxes
    for(QDoubleSpinBox *x: findChildren<QDoubleSpinBox*>())
	{
		connect(
			x,
			SIGNAL(valueChanged(double)),
			this,
			SLOT(on_doubleSpinBox_changed(double))
			);
	}

	//connect all double spin boxes
    for(QSpinBox *x: findChildren<QSpinBox*>())
	{
		connect(
			x,
			SIGNAL(valueChanged(double)),
			this,
			SLOT(on_spinBox_changed(double))
			);
	}
}

