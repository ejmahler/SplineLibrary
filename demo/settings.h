#ifndef SETTINGS_H
#define SETTINGS_H

#include <QObject>
#include <QSettings>

class Settings : public QObject
{
	Q_OBJECT

public:
	~Settings();

	static bool hasSetting(const QString& category, const QString &name);
	static QVariant getSetting(const QString& category, const QString &name);
	static void saveSetting(const QString& category, const QString &name, const QVariant &value);
	static void clearSetting(const QString& category, const QString &name);

private:
	static inline QString getKey(const QString& category, const QString &name);

	Settings(QObject *parent);

	static QSettings settings;
};

#endif // SETTINGS_H
