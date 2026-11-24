#include "SolverProcess.h"

SolverProcess::SolverProcess(QObject *parent) : QObject(parent) {
  connect(&process_, &QProcess::readyReadStandardOutput, this, &SolverProcess::onReadyStdout);
  connect(&process_, &QProcess::readyReadStandardError, this, &SolverProcess::onReadyStderr);
  connect(&process_, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished), this,
          &SolverProcess::onFinished);
}

void SolverProcess::setSolverPath(const QString &path) { solver_path_ = path; }

bool SolverProcess::isBusy() const {
  return stage_ == Stage::Validating || stage_ == Stage::Preparing || stage_ == Stage::Running ||
         stage_ == Stage::Stopping;
}

void SolverProcess::setStage(Stage s) {
  if (stage_ == s) return;
  stage_ = s;
  emit stageChanged(s);
}

void SolverProcess::startCommand(const QStringList &args, Stage stage) {
  if (isBusy()) {
    emit logLine(QStringLiteral("已有任务在运行，拒绝启动新任务"));
    return;
  }
  if (solver_path_.isEmpty()) {
    emit logLine(QStringLiteral("未设置 coldcrucible-solver 路径"));
    setStage(Stage::Failed);
    emit finished(2, Stage::Failed);
    return;
  }
  stdout_buf_.clear();
  setStage(stage);
  process_.setProgram(solver_path_);
  process_.setArguments(args);
  process_.setProcessChannelMode(QProcess::SeparateChannels);
  process_.start();
  if (!process_.waitForStarted(5000)) {
    emit logLine(QStringLiteral("进程启动失败: %1").arg(process_.errorString()));
    setStage(Stage::Failed);
    emit finished(2, Stage::Failed);
  }
}

void SolverProcess::validate(const QString &casePath) {
  startCommand({QStringLiteral("validate"), QStringLiteral("--config"), casePath}, Stage::Validating);
}

void SolverProcess::prepare(const QString &casePath, const QString &outDir) {
  startCommand({QStringLiteral("prepare"), QStringLiteral("--config"), casePath, QStringLiteral("--output"), outDir},
               Stage::Preparing);
}

void SolverProcess::run(const QString &casePath, const QString &preparedDir, const QString &runDir) {
  startCommand({QStringLiteral("run"), QStringLiteral("--config"), casePath, QStringLiteral("--prepared"), preparedDir,
                QStringLiteral("--output"), runDir},
               Stage::Running);
}

void SolverProcess::requestStop() {
  if (!isBusy()) return;
  setStage(Stage::Stopping);
  emit logLine(QStringLiteral("请求停止求解进程…"));
  process_.terminate();
  if (!process_.waitForFinished(8000)) {
    emit logLine(QStringLiteral("强制结束求解进程"));
    process_.kill();
  }
}

void SolverProcess::onReadyStdout() {
  stdout_buf_.append(process_.readAllStandardOutput());
  while (true) {
    const int nl = stdout_buf_.indexOf('\n');
    if (nl < 0) break;
    const QByteArray line = stdout_buf_.left(nl).trimmed();
    stdout_buf_.remove(0, nl + 1);
    if (line.isEmpty()) continue;
    emit eventLine(QString::fromUtf8(line));
  }
}

void SolverProcess::onReadyStderr() {
  const QByteArray raw = process_.readAllStandardError();
  for (const auto &line : QString::fromUtf8(raw).split(QLatin1Char('\n'))) {
    if (!line.trimmed().isEmpty()) emit logLine(line);
  }
}

void SolverProcess::onFinished(int code, QProcess::ExitStatus st) {
  Q_UNUSED(st);
  // flush remaining stdout
  onReadyStdout();
  Stage final = Stage::Failed;
  if (stage_ == Stage::Stopping) {
    final = Stage::Cancelled;
  } else if (code == 0) {
    final = Stage::Completed;
  }
  setStage(final);
  emit finished(code, final);
}
