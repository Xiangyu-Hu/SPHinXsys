#ifndef ELECTROMAGNETIC_OPHELIE_PROGRESS_H
#define ELECTROMAGNETIC_OPHELIE_PROGRESS_H

#include <chrono>
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
    std::cout << "[ophelie] VTP/reload files go to ./output under your **current working directory**." << std::endl;
    std::cout << "[ophelie] Recommended: cd ~/SPHinXsysSYCL/build && ./tests/.../bin/test_3d_ophelie_team7 ..." << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PROGRESS_H
