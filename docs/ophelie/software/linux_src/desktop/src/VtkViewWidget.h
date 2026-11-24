#pragma once

#include <QWidget>
#include <QStringList>
#include <memory>

class QVTKOpenGLNativeWidget;
class vtkRenderer;
class vtkGenericOpenGLRenderWindow;
class vtkActor;

class VtkViewWidget : public QWidget {
  Q_OBJECT
public:
  explicit VtkViewWidget(QWidget *parent = nullptr);

  void clearActors();
  void loadStl(const QString &path, const QString &name, double r, double g, double b, double opacity = 1.0);
  void resetCamera();
  void setVisibleRole(const QString &name, bool visible);

private:
  QVTKOpenGLNativeWidget *vtk_widget_ = nullptr;
  vtkGenericOpenGLRenderWindow *window_ = nullptr;
  vtkRenderer *renderer_ = nullptr;
  struct NamedActor {
    QString name;
    vtkActor *actor = nullptr;
  };
  QList<NamedActor> actors_;
};
