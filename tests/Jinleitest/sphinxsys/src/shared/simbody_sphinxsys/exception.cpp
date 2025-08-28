#include "exception.h"

namespace SPH
{
//=================================================================================================//
Exception::Exception(const std::string &aMsg, const std::string &aFile, int aLine) : exception()
{
    setNull();

    setMessage(aMsg);
    file_ = aFile;
    line_ = aLine;
}
//=============================================================================================//
namespace
{
std::string findFileName(const std::string &filepath)
{
    std::string::size_type pos = filepath.find_last_of("/\\");
    if (pos + 1 >= filepath.size())
        pos = 0;
    return filepath.substr(pos + 1);
}
} // namespace
//=============================================================================================//
Exception::Exception(const std::string &file, size_t line,
                     const std::string &func)
{
    addMessage("\tThrown at " + findFileName(file) + ":" +
               std::to_string(line) + " in " + func + "().");
}
//=============================================================================================//
Exception::Exception(const std::string &file, size_t line,
                     const std::string &func, const std::string &msg)
    : Exception{file, line, func}
{
    addMessage(msg);
}
//=============================================================================================//
void Exception::addMessage(const std::string &msg)
{
    if (msg_.length() == 0)
        msg_ = msg;
    else
        msg_ = msg + "\n" + msg_;
}
//=============================================================================================//
const char *Exception::what() const noexcept
{
    return getMessage();
}
//=============================================================================================//
void Exception::setNull()
{
    setMessage("");
    line_ = -1;
}
//=============================================================================================//
void Exception::setMessage(const std::string &aMsg)
{
    msg_ = aMsg;
}
//=============================================================================================//
const char *Exception::getMessage() const
{
    return (msg_.c_str());
}
//=============================================================================================//
void Exception::print(std::ostream &aOut) const
{
    aOut << "\nException:\n";

    if (file_.size() > 0)
    {
        aOut << "  file= " << file_ << '\n';
    }

    if (line_ >= 0)
    {
        aOut << "  line= " << line_ << '\n';
    }
    aOut << '\n'
         << std::endl;
}
//=============================================================================================//
} // namespace SPH