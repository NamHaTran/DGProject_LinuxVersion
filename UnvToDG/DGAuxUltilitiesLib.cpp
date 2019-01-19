#include "DGAuxUltilitiesLib.h"
#include <QString>
#include <QDir>

namespace auxUlti
{
    QString workingdir()
	{
        QString workingDirectory(QDir::currentPath());
        return workingDirectory;
	}
}
