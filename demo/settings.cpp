#include "settings.h"

QSettings Settings::settings(QSettings::IniFormat, QSettings::UserScope, "Splines", "Splines");

Settings::Settings(QObject *parent)
	: QObject(parent)
{

}

Settings::~Settings()
{

}

bool Settings::hasSetting(const QString& category, const QString &name)
{
	return settings.contains(getKey(category,name));
}

QVariant Settings::getSetting(const QString& category, const QString &name)
{
	return settings.value(getKey(category,name), 0);
}

void Settings::saveSetting(const QString& category, const QString &name, const QVariant &value)
{
	settings.setValue(getKey(category,name),value);
}

void Settings::clearSetting(const QString& category, const QString &name)
{
	settings.remove(getKey(category,name));
}



QString Settings::getKey(const QString& category, const QString &name)
{
	return QString("%1/%2").arg(category,name);
}
