#include "MainWindow.h"
#include "VtkViewWidget.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QDateTime>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTextEdit>
#include <QToolBar>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QVBoxLayout>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
  setWindowTitle(QStringLiteral("OPHELIE 冷坩埚桌面（Linux 开发版）"));
  resize(1400, 900);
  solver_.setSolverPath(defaultSolverPath());

  auto *tb = addToolBar(QStringLiteral("main"));
  tb->addAction(QStringLiteral("打开项目"), this, &MainWindow::openProject);
  tb->addAction(QStringLiteral("示例:自然对流"), this, &MainWindow::openExampleNc);
  tb->addAction(QStringLiteral("示例:搅拌"), this, &MainWindow::openExampleStir);
  tb->addAction(QStringLiteral("保存"), this, &MainWindow::saveProject);
  tb->addAction(QStringLiteral("另存项目…"), this, &MainWindow::saveProjectAs);
  tb->addSeparator();
  tb->addAction(QStringLiteral("检查"), this, &MainWindow::validateProject);
  tb->addAction(QStringLiteral("生成粒子"), this, &MainWindow::prepareParticles);
  tb->addAction(QStringLiteral("运行"), this, &MainWindow::runSolver);
  tb->addAction(QStringLiteral("停止"), this, &MainWindow::stopSolver);
  tb->addSeparator();
  tb->addAction(QStringLiteral("打开结果目录…"), this, &MainWindow::openRunDirectory);

  auto *split = new QSplitter(this);
  setCentralWidget(split);

  auto *left = new QWidget;
  auto *left_layout = new QVBoxLayout(left);
  tree_ = new QTreeWidget;
  tree_->setHeaderLabels({QStringLiteral("装配 / 角色"), QStringLiteral("状态")});
  left_layout->addWidget(tree_, 2);
  auto *stl_row = new QHBoxLayout;
  auto *replace_btn = new QPushButton(QStringLiteral("浏览替换选中 STL…"));
  stl_row->addWidget(replace_btn);
  left_layout->addLayout(stl_row);
  connect(replace_btn, &QPushButton::clicked, this, &MainWindow::replaceSelectedStl);

  auto *params = new QGroupBox(QStringLiteral("物理设置"));
  auto *form = new QFormLayout(params);
  stir_box_ = new QCheckBox(QStringLiteral("启用机械搅拌"));
  melt_spin_ = new QDoubleSpinBox;
  melt_spin_->setRange(0.01, 1.0);
  melt_spin_->setDecimals(4);
  melt_spin_->setSuffix(QStringLiteral(" m"));
  dp_spin_ = new QDoubleSpinBox;
  dp_spin_->setRange(0.001, 0.05);
  dp_spin_->setDecimals(4);
  dp_spin_->setSuffix(QStringLiteral(" m"));
  end_spin_ = new QDoubleSpinBox;
  end_spin_->setRange(0.001, 100.0);
  end_spin_->setDecimals(4);
  end_spin_->setSuffix(QStringLiteral(" s"));
  rpm_spin_ = new QDoubleSpinBox;
  rpm_spin_->setRange(0.0, 1000.0);
  rpm_spin_->setDecimals(2);
  rpm_spin_->setSuffix(QStringLiteral(" rpm"));
  form->addRow(stir_box_);
  form->addRow(QStringLiteral("液深"), melt_spin_);
  form->addRow(QStringLiteral("粒子间距"), dp_spin_);
  form->addRow(QStringLiteral("结束时间"), end_spin_);
  form->addRow(QStringLiteral("转速"), rpm_spin_);
  left_layout->addWidget(params, 1);

  auto *paths_box = new QGroupBox(QStringLiteral("目录（安装包只读；项目/结果可写）"));
  auto *paths_form = new QFormLayout(paths_box);
  results_root_edit_ = new QLineEdit;
  results_root_edit_->setPlaceholderText(QStringLiteral("默认: <项目>/runs （可改到任意可写目录）"));
  auto *path_btns = new QHBoxLayout;
  auto *browse_out = new QPushButton(QStringLiteral("选择结果目录…"));
  auto *clear_out = new QPushButton(QStringLiteral("恢复默认"));
  path_btns->addWidget(browse_out);
  path_btns->addWidget(clear_out);
  paths_form->addRow(QStringLiteral("结果输出根目录"), results_root_edit_);
  paths_form->addRow(path_btns);
  left_layout->addWidget(paths_box, 0);
  connect(browse_out, &QPushButton::clicked, this, &MainWindow::browseResultsRoot);
  connect(clear_out, &QPushButton::clicked, this, &MainWindow::clearResultsRoot);

  scope_ = new QTextEdit;
  scope_->setReadOnly(true);
  left_layout->addWidget(scope_, 1);
  split->addWidget(left);

  auto *center = new QWidget;
  auto *center_layout = new QVBoxLayout(center);
  center_tabs_ = new QTabWidget;
  view_ = new VtkViewWidget;
  center_tabs_->addTab(view_, QStringLiteral("几何 / 粒子 / 结果"));
  center_layout->addWidget(center_tabs_, 3);
  log_ = new QTextEdit;
  log_->setReadOnly(true);
  center_layout->addWidget(log_, 1);
  split->addWidget(center);
  split->setStretchFactor(0, 1);
  split->setStretchFactor(1, 3);

  status_mode_ = new QLabel(QStringLiteral("模式: —"));
  status_stage_ = new QLabel(QStringLiteral("状态: idle"));
  statusBar()->addWidget(status_mode_);
  statusBar()->addPermanentWidget(status_stage_);

  connect(stir_box_, &QCheckBox::toggled, this, &MainWindow::onStirringToggled);
  connect(&solver_, &SolverProcess::stageChanged, this, &MainWindow::onStageChanged);
  connect(&solver_, &SolverProcess::eventLine, this, &MainWindow::onEventLine);
  connect(&solver_, &SolverProcess::logLine, this, &MainWindow::onLogLine);
  connect(&solver_, &SolverProcess::finished, this, &MainWindow::onSolverFinished);
}

static QString linuxSrcRoot() {
  const QByteArray env = qgetenv("OPHELIE_LINUX_SRC");
  if (!env.isEmpty()) return QFileInfo(QString::fromUtf8(env)).absoluteFilePath();
  // Typical: <linux_src>/desktop/build/ophelie_desktop
  QDir d(QCoreApplication::applicationDirPath());
  const QString cand = d.absoluteFilePath(QStringLiteral("../.."));
  if (QFileInfo::exists(cand + QStringLiteral("/bin/coldcrucible-solver")) &&
      QFileInfo::exists(cand + QStringLiteral("/examples/natural-convection/case.json"))) {
    return QFileInfo(cand).absoluteFilePath();
  }
  // Fallback: cwd-relative docs/ophelie/software/linux_src (when launched from repo root)
  const QString from_cwd = QDir::current().absoluteFilePath(QStringLiteral("docs/ophelie/software/linux_src"));
  if (QFileInfo::exists(from_cwd + QStringLiteral("/bin/coldcrucible-solver"))) {
    return QFileInfo(from_cwd).absoluteFilePath();
  }
  return {};
}

QString MainWindow::defaultSolverPath() const {
  const QByteArray env = qgetenv("OPHELIE_SOLVER");
  if (!env.isEmpty()) return QString::fromUtf8(env);
  const QString root = linuxSrcRoot();
  if (!root.isEmpty()) {
    const QString p = root + QStringLiteral("/bin/coldcrucible-solver");
    if (QFileInfo::exists(p)) return QFileInfo(p).absoluteFilePath();
  }
  return {};
}

void MainWindow::openProject() {
  const QString path = QFileDialog::getOpenFileName(this, QStringLiteral("打开 case.json"), {},
                                                    QStringLiteral("Case (case.json)"));
  if (!path.isEmpty()) loadProjectPath(path);
}

void MainWindow::openExampleNc() {
  const QString root = linuxSrcRoot();
  if (root.isEmpty()) {
    QMessageBox::warning(this, QStringLiteral("找不到示例"),
                         QStringLiteral("请设置 OPHELIE_LINUX_SRC，或从仓库根目录启动桌面程序。"));
    return;
  }
  loadProjectPath(root + QStringLiteral("/examples/natural-convection/case.json"));
}

void MainWindow::openExampleStir() {
  const QString root = linuxSrcRoot();
  if (root.isEmpty()) {
    QMessageBox::warning(this, QStringLiteral("找不到示例"),
                         QStringLiteral("请设置 OPHELIE_LINUX_SRC，或从仓库根目录启动桌面程序。"));
    return;
  }
  loadProjectPath(root + QStringLiteral("/examples/mechanical-stirring/case.json"));
}

void MainWindow::loadProjectPath(const QString &path) {
  QString err;
  if (!project_.load(path, &err)) {
    QMessageBox::warning(this, QStringLiteral("打开失败"), err);
    return;
  }
  refreshUiFromProject();
  reloadAssemblyView();
  onLogLine(QStringLiteral("已打开: %1").arg(path));
}

void MainWindow::refreshUiFromProject() {
  stir_box_->blockSignals(true);
  stir_box_->setChecked(project_.stirringEnabled());
  stir_box_->blockSignals(false);
  melt_spin_->setValue(project_.meltHeightM());
  dp_spin_->setValue(project_.dpM());
  end_spin_->setValue(project_.endTimeS());
  rpm_spin_->setValue(project_.rpm());
  rpm_spin_->setEnabled(project_.stirringEnabled());
  results_root_edit_->setText(project_.resultsRoot());
  status_mode_->setText(project_.stirringEnabled() ? QStringLiteral("模式: 机械搅拌 / elastic_driven")
                                                   : QStringLiteral("模式: 自然对流"));
  scope_->setPlainText(project_.acceptanceSummary() + QStringLiteral("\n\n目录约定:\n") +
                       QStringLiteral("- 安装目录: 程序/求解器（只读）\n") +
                       QStringLiteral("- 项目目录: case.json / assets / prepared\n") +
                       QStringLiteral("- 结果目录: 默认可在项目下，也可指定到任意可写路径"));
  setWindowTitle(QStringLiteral("OPHELIE — %1").arg(project_.projectName()));

  tree_->clear();
  for (const auto &role : project_.assemblyRoles()) {
    auto *item = new QTreeWidgetItem(tree_);
    item->setText(0, role);
    const QString stl = project_.roleStl(role);
    item->setText(1, QFileInfo::exists(stl) ? QStringLiteral("OK") : QStringLiteral("缺失"));
    item->setToolTip(1, stl);
    item->setData(0, Qt::UserRole, role);
  }
}

void MainWindow::reloadAssemblyView() {
  view_->clearActors();
  const struct {
    const char *role;
    double r, g, b, a;
  } colors[] = {
      {"crucible", 0.75, 0.75, 0.8, 0.35},
      {"coil", 0.9, 0.55, 0.15, 0.9},
      {"cooling_pipes", 0.2, 0.6, 0.9, 0.5},
      {"thermocouple", 0.2, 0.9, 0.3, 1.0},
      {"paddle", 0.85, 0.85, 0.2, 1.0},
      {"paddle_solver", 0.95, 0.3, 0.3, 0.7},
  };
  for (const auto &c : colors) {
    const QString path = project_.roleStl(QString::fromUtf8(c.role));
    if (path.isEmpty()) continue;
    if (!project_.stirringEnabled() && (QString::fromUtf8(c.role).startsWith(QStringLiteral("paddle")))) {
      // preview-only possible later; still show paddle if present but dimmer
      view_->loadStl(path, QString::fromUtf8(c.role), c.r, c.g, c.b, 0.25);
    } else {
      view_->loadStl(path, QString::fromUtf8(c.role), c.r, c.g, c.b, c.a);
    }
  }
  view_->resetCamera();
}

void MainWindow::applyParamsToProject() {
  project_.setStirringEnabled(stir_box_->isChecked());
  project_.setMeltHeightM(melt_spin_->value());
  project_.setDpM(dp_spin_->value());
  project_.setEndTimeS(end_spin_->value());
  if (stir_box_->isChecked()) project_.setRpm(rpm_spin_->value());
  project_.setResultsRoot(results_root_edit_->text().trimmed());
}

void MainWindow::saveProject() {
  applyParamsToProject();
  if (!results_root_edit_->text().trimmed().isEmpty()) {
    project_.setResultsRoot(results_root_edit_->text().trimmed());
  }
  QString err;
  if (!project_.save(&err)) {
    QMessageBox::warning(this, QStringLiteral("保存失败"), err);
    return;
  }
  onLogLine(QStringLiteral("已保存 case.json"));
  refreshUiFromProject();
}

void MainWindow::saveProjectAs() {
  const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("另存项目为 case.json"),
                                                    QDir::homePath() + QStringLiteral("/ophelie_project/case.json"),
                                                    QStringLiteral("Case (case.json)"));
  if (path.isEmpty()) return;
  applyParamsToProject();
  // Ensure destination directory exists; copy current JSON there then reload.
  QDir().mkpath(QFileInfo(path).absolutePath());
  QFile out(path);
  if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    QMessageBox::warning(this, QStringLiteral("另存失败"), path);
    return;
  }
  out.write(QJsonDocument(project_.root()).toJson(QJsonDocument::Indented));
  out.close();
  loadProjectPath(path);
  onLogLine(QStringLiteral("已另存并打开: %1").arg(path));
}

void MainWindow::browseResultsRoot() {
  const QString dir = QFileDialog::getExistingDirectory(this, QStringLiteral("选择结果输出根目录"),
                                                        results_root_edit_->text().isEmpty() ? QDir::homePath()
                                                                                             : results_root_edit_->text());
  if (dir.isEmpty()) return;
  results_root_edit_->setText(dir);
  project_.setResultsRoot(dir);
  onLogLine(QStringLiteral("结果将写到: %1/<run-id>/").arg(dir));
}

void MainWindow::clearResultsRoot() {
  results_root_edit_->clear();
  project_.setResultsRoot(QString());
  onLogLine(QStringLiteral("结果目录恢复为默认 <项目>/runs/"));
}

void MainWindow::openRunDirectory() {
  QString start = results_root_edit_->text().trimmed();
  if (start.isEmpty() && !project_.projectDir().isEmpty()) {
    start = QDir(project_.projectDir()).filePath(QStringLiteral("runs"));
  }
  if (start.isEmpty()) start = QDir::homePath();
  const QString path = QFileDialog::getExistingDirectory(this, QStringLiteral("打开结果/run 目录（含 result-index.json）"), start);
  if (path.isEmpty()) return;
  const QString index_path = QDir(path).filePath(QStringLiteral("result-index.json"));
  if (!QFileInfo::exists(index_path)) {
    QMessageBox::warning(this, QStringLiteral("无效结果目录"),
                         QStringLiteral("未找到 result-index.json:\n%1").arg(index_path));
    return;
  }
  QFile f(index_path);
  if (!f.open(QIODevice::ReadOnly)) return;
  const auto doc = QJsonDocument::fromJson(f.readAll());
  const auto o = doc.object();
  onLogLine(QStringLiteral("打开结果: %1").arg(path));
  onLogLine(QStringLiteral("  status=%1 physics=%2 passed=%3")
                .arg(o.value(QStringLiteral("status")).toString())
                .arg(o.value(QStringLiteral("physics")).toString())
                .arg(o.value(QStringLiteral("summary_passed")).toInt(-1)));
  // Load available VTP frames into viewer (geometry/results tab).
  view_->clearActors();
  const auto frames = o.value(QStringLiteral("frames")).toArray();
  int loaded = 0;
  for (const auto &v : frames) {
    const QString rel = v.toString();
    if (!rel.endsWith(QStringLiteral(".vtp"))) continue;
    // Prefer a mid/last GlassBody frame if many; load up to a few bodies' last frame.
    if (!rel.contains(QStringLiteral("GlassBody")) && !rel.contains(QStringLiteral("Rotor")) &&
        !rel.contains(QStringLiteral("Wall"))) {
      continue;
    }
  }
  // Load last frame per body from index.bodies if present.
  const auto bodies = o.value(QStringLiteral("bodies")).toObject();
  double colors[][3] = {{0.3, 0.7, 0.95}, {0.95, 0.55, 0.2}, {0.6, 0.6, 0.6}, {0.2, 0.9, 0.4}};
  int ci = 0;
  for (auto it = bodies.begin(); it != bodies.end(); ++it) {
    const auto files = it.value().toObject().value(QStringLiteral("files")).toArray();
    if (files.isEmpty()) continue;
    const QString rel = files.last().toString();
    const QString abs = QDir(path).filePath(rel);
    if (!QFileInfo::exists(abs)) continue;
    const auto *c = colors[ci % 4];
    view_->loadStl(abs, it.key(), c[0], c[1], c[2], 0.85);
    ++ci;
    ++loaded;
  }
  if (loaded == 0) {
    // fallback: any vtp under frames/
    const QDir frames_dir(QDir(path).filePath(QStringLiteral("frames")));
    const auto list = frames_dir.entryList(QStringList() << QStringLiteral("*.vtp"), QDir::Files, QDir::Name);
    for (const auto &name : list) {
      if (!name.contains(QStringLiteral("GlassBody"))) continue;
      view_->loadStl(frames_dir.filePath(name), name, 0.3, 0.7, 0.95, 0.85);
      ++loaded;
      break;
    }
  }
  view_->resetCamera();
  onLogLine(QStringLiteral("  已加载 %1 个结果体到视窗（需图形会话查看）").arg(loaded));
  if (!o.value(QStringLiteral("claims_scope")).toArray().isEmpty()) {
    onLogLine(QStringLiteral("  注意: 请遵守 result-index 中的 claims_scope，勿夸大物理结论"));
  }
}

void MainWindow::replaceSelectedStl() {
  if (project_.casePath().isEmpty()) {
    QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("请先打开项目"));
    return;
  }
  auto *item = tree_->currentItem();
  if (!item) {
    QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("请先在左侧选中一个装配角色"));
    return;
  }
  const QString role = item->data(0, Qt::UserRole).toString();
  const QString path = QFileDialog::getOpenFileName(this, QStringLiteral("选择 STL（浏览本机硬盘）"), QDir::homePath(),
                                                    QStringLiteral("STL (*.stl *.STL);;All (*)"));
  if (path.isEmpty()) return;
  QString err;
  if (!project_.setRoleStl(role, path, &err)) {
    QMessageBox::warning(this, QStringLiteral("导入失败"), err);
    return;
  }
  project_.save(nullptr);
  refreshUiFromProject();
  reloadAssemblyView();
  onLogLine(QStringLiteral("已导入 %1 <= %2（已复制到项目 assets/）").arg(role, path));
}

void MainWindow::onStirringToggled(bool on) {
  rpm_spin_->setEnabled(on);
  applyParamsToProject();
  refreshUiFromProject();
  reloadAssemblyView();
  if (on && project_.roleStl(QStringLiteral("paddle")).isEmpty()) {
    QMessageBox::warning(this, QStringLiteral("缺桨"), QStringLiteral("启用搅拌但缺少桨 STL；将阻止 prepare/run。"));
  }
}

void MainWindow::validateProject() {
  if (project_.casePath().isEmpty()) return;
  applyParamsToProject();
  project_.save(nullptr);
  setBusy(true);
  solver_.validate(project_.casePath());
}

void MainWindow::prepareParticles() {
  if (project_.casePath().isEmpty()) return;
  applyParamsToProject();
  project_.save(nullptr);
  QDir().mkpath(projectPreparedDir());
  setBusy(true);
  solver_.prepare(project_.casePath(), projectPreparedDir());
}

void MainWindow::runSolver() {
  if (project_.casePath().isEmpty()) return;
  applyParamsToProject();
  project_.save(nullptr);
  const QString run = newRunDir();
  QDir().mkpath(run);
  setBusy(true);
  solver_.run(project_.casePath(), projectPreparedDir(), run);
}

void MainWindow::stopSolver() { solver_.requestStop(); }

QString MainWindow::projectPreparedDir() const {
  return QDir(project_.projectDir()).filePath(QStringLiteral("prepared"));
}

QString MainWindow::newRunDir() const {
  const QString id = QDateTime::currentDateTimeUtc().toString(QStringLiteral("yyyyMMdd_HHmmss"));
  const QString root = project_.resultsRoot().trimmed().isEmpty()
                           ? QDir(project_.projectDir()).filePath(QStringLiteral("runs"))
                           : project_.resultsRoot().trimmed();
  return QDir(root).filePath(id);
}

void MainWindow::setBusy(bool busy) {
  // toolbar actions remain clickable for stop; parameter edits locked while busy
  melt_spin_->setEnabled(!busy);
  dp_spin_->setEnabled(!busy);
  end_spin_->setEnabled(!busy);
  stir_box_->setEnabled(!busy);
  rpm_spin_->setEnabled(!busy && stir_box_->isChecked());
}

void MainWindow::onStageChanged(SolverProcess::Stage stage) {
  static const char *names[] = {"idle", "validating", "preparing", "running", "stopping",
                                "completed", "cancelled", "failed"};
  const int i = static_cast<int>(stage);
  status_stage_->setText(QStringLiteral("状态: %1").arg(QString::fromUtf8(names[i])));
}

void MainWindow::onEventLine(const QString &line) {
  onLogLine(QStringLiteral("[event] %1").arg(line));
  const auto doc = QJsonDocument::fromJson(line.toUtf8());
  if (!doc.isObject()) return;
  const auto o = doc.object();
  if (o.value(QStringLiteral("event")).toString() == QStringLiteral("error")) {
    statusBar()->showMessage(o.value(QStringLiteral("message")).toString(), 8000);
  }
}

void MainWindow::onLogLine(const QString &line) {
  log_->append(line);
}

void MainWindow::onSolverFinished(int code, SolverProcess::Stage stage) {
  setBusy(false);
  onLogLine(QStringLiteral("求解器结束 exit=%1 stage=%2").arg(code).arg(static_cast<int>(stage)));
  if (code != 0) {
    QMessageBox::warning(this, QStringLiteral("任务失败"),
                         QStringLiteral("exit=%1。详见日志；勿将失败标为通过。").arg(code));
  }
}
