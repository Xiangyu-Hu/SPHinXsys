#include "VtkViewWidget.h"

#include <QVBoxLayout>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>

VtkViewWidget::VtkViewWidget(QWidget *parent) : QWidget(parent) {
  auto *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  vtk_widget_ = new QVTKOpenGLNativeWidget(this);
  layout->addWidget(vtk_widget_);

  window_ = vtkGenericOpenGLRenderWindow::New();
  vtk_widget_->setRenderWindow(window_);
  renderer_ = vtkRenderer::New();
  renderer_->SetBackground(0.12, 0.14, 0.18);
  window_->AddRenderer(renderer_);
}

void VtkViewWidget::clearActors() {
  for (auto &a : actors_) {
    if (a.actor) {
      renderer_->RemoveActor(a.actor);
      a.actor->Delete();
    }
  }
  actors_.clear();
  window_->Render();
}

void VtkViewWidget::loadStl(const QString &path, const QString &name, double r, double g, double b, double opacity) {
  if (path.isEmpty()) return;
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  if (path.endsWith(QStringLiteral(".vtp"), Qt::CaseInsensitive)) {
    auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(path.toUtf8().constData());
    reader->Update();
    mapper->SetInputConnection(reader->GetOutputPort());
  } else {
    auto reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(path.toUtf8().constData());
    reader->Update();
    mapper->SetInputConnection(reader->GetOutputPort());
  }
  auto *actor = vtkActor::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(r, g, b);
  actor->GetProperty()->SetOpacity(opacity);
  renderer_->AddActor(actor);
  actors_.push_back({name, actor});
  window_->Render();
}

void VtkViewWidget::resetCamera() {
  renderer_->ResetCamera();
  window_->Render();
}

void VtkViewWidget::setVisibleRole(const QString &name, bool visible) {
  for (auto &a : actors_) {
    if (a.name == name && a.actor) a.actor->SetVisibility(visible ? 1 : 0);
  }
  window_->Render();
}
