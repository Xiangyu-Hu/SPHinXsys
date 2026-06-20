#ifndef ELECTROMAGNETIC_OPHELIE_PROGRESS_H
#define ELECTROMAGNETIC_OPHELIE_PROGRESS_H

#include "io_environment.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

class OphelieProgressLogger
{
  public:
    explicit OphelieProgressLogger(const std::string &stage_name)
        : stage_name_(stage_name), start_(std::chrono::steady_clock::now())
    {
        std::cout << "[ophelie] >>> " << stage_name_ << " started" << std::endl;
        std::cout.flush();
    }

    void log(const std::string &message) const
    {
        std::cout << "[ophelie] " << stage_name_ << " | " << message << " | elapsed_s=" << elapsedSeconds() << std::endl;
        std::cout.flush();
    }

    void finish(const std::string &message = "done") const
    {
        std::cout << "[ophelie] <<< " << stage_name_ << " " << message << " | total_s=" << elapsedSeconds() << std::endl;
        std::cout.flush();
    }

  private:
    double elapsedSeconds() const
    {
        const auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - start_).count();
    }

    std::string stage_name_;
    std::chrono::steady_clock::time_point start_;
};

inline void logOphelieRunContext()
{
    const std::string cwd = std::filesystem::current_path().string();
    std::cout << "[ophelie] cwd=" << cwd << " | VTP/reload under **./output** (relative to cwd)." << std::endl;
    std::cout << "[ophelie] If ./output already exists, SphinxSys moves it to ./output_backup at startup." << std::endl;
    std::cout << "[ophelie] Typical: cd ~/SPHinXsysSYCL/build then run; relax VTPs may land in output_backup after EM."
              << std::endl;
}

/** Log absolute path of a written artifact (VTP, reload, etc.). */
inline void logOphelieOutputArtifact(const std::string &relative_path)
{
    const std::filesystem::path absolute_path = std::filesystem::absolute(relative_path);
    std::cout << "[ophelie] wrote " << relative_path << " -> " << absolute_path.string() << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PROGRESS_H
