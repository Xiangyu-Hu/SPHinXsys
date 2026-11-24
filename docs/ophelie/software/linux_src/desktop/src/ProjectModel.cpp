#include "ProjectModel.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

bool ProjectModel::load(const QString &casePath, QString *error) {
  QFile f(casePath);
  if (!f.open(QIODevice::ReadOnly)) {
    if (error) *error = QStringLiteral("无法打开: %1").arg(casePath);
    return false;
  }
  const auto doc = QJsonDocument::fromJson(f.readAll());
  if (!doc.isObject()) {
    if (error) *error = QStringLiteral("case.json 不是对象");
    return false;
  }
  root_ = doc.object();
  case_path_ = casePath;
  project_dir_ = QFileInfo(casePath).absolutePath();
  return true;
}

bool ProjectModel::save(QString *error) const {
  if (case_path_.isEmpty()) {
    if (error) *error = QStringLiteral("无项目路径");
    return false;
  }
  QFile f(case_path_);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
    if (error) *error = QStringLiteral("无法写入: %1").arg(case_path_);
    return false;
  }
  f.write(QJsonDocument(root_).toJson(QJsonDocument::Indented));
  return true;
}

QString ProjectModel::projectName() const {
  return root_.value(QStringLiteral("project_name")).toString(QStringLiteral("(未命名)"));
}

bool ProjectModel::stirringEnabled() const {
  return root_.value(QStringLiteral("stirring")).toObject().value(QStringLiteral("enabled")).toBool(false);
}

void ProjectModel::setStirringEnabled(bool on) {
  auto stirring = root_.value(QStringLiteral("stirring")).toObject();
  stirring.insert(QStringLiteral("enabled"), on);
  if (on) {
    stirring.insert(QStringLiteral("mechanics_mode"), QStringLiteral("elastic_driven"));
  } else {
    stirring.insert(QStringLiteral("mechanics_mode"), QStringLiteral("none"));
  }
  root_.insert(QStringLiteral("stirring"), stirring);
}

QString ProjectModel::mechanicsMode() const {
  return root_.value(QStringLiteral("stirring")).toObject().value(QStringLiteral("mechanics_mode")).toString(QStringLiteral("none"));
}

double ProjectModel::meltHeightM() const {
  return root_.value(QStringLiteral("geometry")).toObject().value(QStringLiteral("melt_height_m")).toDouble(0.18);
}

void ProjectModel::setMeltHeightM(double v) {
  auto g = root_.value(QStringLiteral("geometry")).toObject();
  g.insert(QStringLiteral("melt_height_m"), v);
  root_.insert(QStringLiteral("geometry"), g);
}

double ProjectModel::dpM() const {
  return root_.value(QStringLiteral("numerics")).toObject().value(QStringLiteral("dp_m")).toDouble(0.01);
}

void ProjectModel::setDpM(double v) {
  auto n = root_.value(QStringLiteral("numerics")).toObject();
  n.insert(QStringLiteral("dp_m"), v);
  root_.insert(QStringLiteral("numerics"), n);
}

double ProjectModel::endTimeS() const {
  return root_.value(QStringLiteral("run")).toObject().value(QStringLiteral("end_time_s")).toDouble(0.1);
}

void ProjectModel::setEndTimeS(double v) {
  auto r = root_.value(QStringLiteral("run")).toObject();
  r.insert(QStringLiteral("end_time_s"), v);
  root_.insert(QStringLiteral("run"), r);
}

double ProjectModel::rpm() const {
  return root_.value(QStringLiteral("stirring")).toObject().value(QStringLiteral("rpm")).toDouble(0.0);
}

void ProjectModel::setRpm(double v) {
  auto s = root_.value(QStringLiteral("stirring")).toObject();
  s.insert(QStringLiteral("rpm"), v);
  root_.insert(QStringLiteral("stirring"), s);
}

QString ProjectModel::roleStl(const QString &role) const {
  const auto rel = roleStlRelative(role);
  if (rel.isEmpty()) return {};
  const QFileInfo fi(rel);
  if (fi.isAbsolute()) return fi.absoluteFilePath();
  return QFileInfo(project_dir_ + QLatin1Char('/') + rel).absoluteFilePath();
}

QString ProjectModel::roleStlRelative(const QString &role) const {
  const auto roles = root_.value(QStringLiteral("assembly")).toObject().value(QStringLiteral("roles")).toObject();
  return roles.value(role).toString();
}

bool ProjectModel::setRoleStl(const QString &role, const QString &absoluteOrRelativePath, QString *error) {
  if (project_dir_.isEmpty()) {
    if (error) *error = QStringLiteral("请先打开或保存项目");
    return false;
  }
  QFileInfo src(absoluteOrRelativePath);
  if (!src.exists() || !src.isFile()) {
    if (error) *error = QStringLiteral("STL 不存在: %1").arg(absoluteOrRelativePath);
    return false;
  }
  // Prefer copying into project assets/ so project stays portable; keep original name.
  const QString assets = project_dir_ + QStringLiteral("/assets");
  QDir().mkpath(assets);
  const QString dest = assets + QLatin1Char('/') + src.fileName();
  if (QFileInfo(dest).absoluteFilePath() != src.absoluteFilePath()) {
    QFile::remove(dest);
    if (!QFile::copy(src.absoluteFilePath(), dest)) {
      if (error) *error = QStringLiteral("无法复制到项目 assets/: %1").arg(dest);
      return false;
    }
  }
  auto assembly = root_.value(QStringLiteral("assembly")).toObject();
  auto roles = assembly.value(QStringLiteral("roles")).toObject();
  roles.insert(role, QStringLiteral("assets/") + src.fileName());
  assembly.insert(QStringLiteral("roles"), roles);
  root_.insert(QStringLiteral("assembly"), assembly);
  return true;
}

QStringList ProjectModel::assemblyRoles() const {
  const auto roles = root_.value(QStringLiteral("assembly")).toObject().value(QStringLiteral("roles")).toObject();
  return roles.keys();
}

QString ProjectModel::resultsRoot() const {
  return root_.value(QStringLiteral("paths")).toObject().value(QStringLiteral("results_root")).toString();
}

void ProjectModel::setResultsRoot(const QString &absolutePath) {
  auto paths = root_.value(QStringLiteral("paths")).toObject();
  if (absolutePath.trimmed().isEmpty()) {
    paths.remove(QStringLiteral("results_root"));
  } else {
    paths.insert(QStringLiteral("results_root"), absolutePath.trimmed());
  }
  if (paths.isEmpty()) {
    root_.remove(QStringLiteral("paths"));
  } else {
    root_.insert(QStringLiteral("paths"), paths);
  }
}

QString ProjectModel::purpose() const {
  return root_.value(QStringLiteral("purpose")).toString();
}

QString ProjectModel::acceptanceSummary() const {
  const auto scope = root_.value(QStringLiteral("acceptance_scope")).toObject();
  QStringList lines;
  lines << QStringLiteral("用途: %1").arg(purpose());
  lines << QStringLiteral("模式: %1").arg(stirringEnabled() ? QStringLiteral("机械搅拌") : QStringLiteral("自然对流"));
  if (stirringEnabled()) {
    lines << QStringLiteral("力学: %1（禁止刚体动画冒充）").arg(mechanicsMode());
  }
  const auto forbids = scope.value(QStringLiteral("forbids")).toArray();
  if (!forbids.isEmpty()) {
    lines << QStringLiteral("禁止夸大:");
    for (const auto &v : forbids) lines << QStringLiteral("  - %1").arg(v.toString());
  }
  return lines.join(QLatin1Char('\n'));
}
