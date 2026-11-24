#pragma once

#include <QJsonObject>
#include <QString>
#include <optional>

class ProjectModel {
public:
  bool load(const QString &casePath, QString *error);
  bool save(QString *error) const;

  QString casePath() const { return case_path_; }
  QString projectDir() const { return project_dir_; }
  QString projectName() const;
  QJsonObject root() const { return root_; }

  bool stirringEnabled() const;
  void setStirringEnabled(bool on);
  QString mechanicsMode() const;

  double meltHeightM() const;
  void setMeltHeightM(double v);
  double dpM() const;
  void setDpM(double v);
  double endTimeS() const;
  void setEndTimeS(double v);
  double rpm() const;
  void setRpm(double v);

  QString roleStl(const QString &role) const;
  QString roleStlRelative(const QString &role) const;
  bool setRoleStl(const QString &role, const QString &absoluteOrRelativePath, QString *error);
  QStringList assemblyRoles() const;

  QString purpose() const;
  QString acceptanceSummary() const;

  /** User-writable results root; empty means <project>/runs */
  QString resultsRoot() const;
  void setResultsRoot(const QString &absolutePath);

private:
  QString case_path_;
  QString project_dir_;
  QJsonObject root_;
};
