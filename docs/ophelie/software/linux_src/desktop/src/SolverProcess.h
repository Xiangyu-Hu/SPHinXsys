#pragma once

#include <QObject>
#include <QProcess>
#include <QString>

class SolverProcess : public QObject {
  Q_OBJECT
public:
  enum class Stage { Idle, Validating, Preparing, Running, Stopping, Completed, Cancelled, Failed };

  explicit SolverProcess(QObject *parent = nullptr);

  void setSolverPath(const QString &path);
  QString solverPath() const { return solver_path_; }

  Stage stage() const { return stage_; }
  bool isBusy() const;

  void validate(const QString &casePath);
  void prepare(const QString &casePath, const QString &outDir);
  void run(const QString &casePath, const QString &preparedDir, const QString &runDir);
  void requestStop();

signals:
  void stageChanged(SolverProcess::Stage stage);
  void eventLine(const QString &jsonLine);
  void logLine(const QString &line);
  void finished(int exitCode, SolverProcess::Stage finalStage);

private:
  void startCommand(const QStringList &args, Stage stage);
  void setStage(Stage s);
  void onReadyStdout();
  void onReadyStderr();
  void onFinished(int code, QProcess::ExitStatus st);

  QString solver_path_;
  QProcess process_;
  Stage stage_ = Stage::Idle;
  QByteArray stdout_buf_;
};
