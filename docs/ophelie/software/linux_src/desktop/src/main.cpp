#include "MainWindow.h"

#include <QApplication>
#include <QSurfaceFormat>

#include <QVTKOpenGLNativeWidget.h>

int main(int argc, char *argv[]) {
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
  QApplication app(argc, argv);
  QApplication::setApplicationName(QStringLiteral("OPHELIE Desktop"));
  QApplication::setOrganizationName(QStringLiteral("OPHELIE"));

  MainWindow w;
  w.show();
  return app.exec();
}
