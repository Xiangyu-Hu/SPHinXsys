#pragma once

#include "ProjectModel.h"
#include "SolverProcess.h"

#include <QMainWindow>

class QTreeWidget;
class QTreeWidgetItem;
class QTextEdit;
class QLabel;
class QCheckBox;
class QDoubleSpinBox;
class QLineEdit;
class QTabWidget;
class VtkViewWidget;

class MainWindow : public QMainWindow {
  Q_OBJECT
public:
  explicit MainWindow(QWidget *parent = nullptr);

private slots:
  void openProject();
  void openExampleNc();
  void openExampleStir();
  void saveProject();
  void saveProjectAs();
  void validateProject();
  void prepareParticles();
  void runSolver();
  void stopSolver();
  void onStirringToggled(bool on);
  void browseResultsRoot();
  void clearResultsRoot();
  void replaceSelectedStl();
  void openRunDirectory();
  void onStageChanged(SolverProcess::Stage stage);
  void onEventLine(const QString &line);
  void onLogLine(const QString &line);
  void onSolverFinished(int code, SolverProcess::Stage stage);

private:
  void loadProjectPath(const QString &path);
  void refreshUiFromProject();
  void reloadAssemblyView();
  void applyParamsToProject();
  QString defaultSolverPath() const;
  QString projectPreparedDir() const;
  QString newRunDir() const;
  void setBusy(bool busy);

  ProjectModel project_;
  SolverProcess solver_;

  QLabel *status_mode_ = nullptr;
  QLabel *status_stage_ = nullptr;
  QTreeWidget *tree_ = nullptr;
  VtkViewWidget *view_ = nullptr;
  QTextEdit *log_ = nullptr;
  QTextEdit *scope_ = nullptr;
  QCheckBox *stir_box_ = nullptr;
  QDoubleSpinBox *melt_spin_ = nullptr;
  QDoubleSpinBox *dp_spin_ = nullptr;
  QDoubleSpinBox *end_spin_ = nullptr;
  QDoubleSpinBox *rpm_spin_ = nullptr;
  QLineEdit *results_root_edit_ = nullptr;
  QTabWidget *center_tabs_ = nullptr;
};
