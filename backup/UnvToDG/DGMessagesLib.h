#ifndef DGMESSAGESLIB_H_INCLUDED
#define DGMESSAGESLIB_H_INCLUDED
#include <QString>
namespace message
{
	/*Function create header of DG's data files*/
    QString headerFile();

	/*Function display undefined keyWord error*/
    QString undfKeyW(QString keyW, QString location);

	/*Function display opening file error*/
    QString opFError(QString fileName, QString location);

	/*Function writes logFile to report error to user*/
    void writeLog(QString location, QString caseName, QString str);

	/*Functions gets time data from system*/
    QString getTime();

	/*Function create header of files p T U*/
    QString headerpTU(QString file);
}
#endif // DGMESSAGESLIB_H_INCLUDED
