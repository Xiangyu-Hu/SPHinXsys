#ifndef LOG_H
#define LOG_H
#include<tuple>
#include<map>
#include<cstdio>
#include<cstdarg>
namespace gpu{
  namespace Logging{
    using ElementType = std::tuple<FILE*, const char*, const char*>;
    enum LogLevel{CRITICAL=0, ERROR, EXCEPTION, WARNING, MESSAGE, STDERR, STDOUT,
		  DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4, DEBUG5, DEBUG6, DEBUG7};

#ifdef MAXLOGLEVEL
    constexpr int maxLogLevel = MAXLOGLEVEL;
#else
    constexpr int maxLogLevel = 6;
#endif

    ElementType getLogLevelInfo(int level){
      static const std::map<int, ElementType> printMap{
	{CRITICAL , std::make_tuple(stderr, "\e[101m[CRITICAL] ", "\e[0m\n")},
	{ERROR , std::make_tuple(stderr, "\e[91m[ERROR] \e[0m", "\n")},
	{EXCEPTION , std::make_tuple(stderr, "\e[1m\e[91m[EXCEPTION] \e[0m", "\n")},
	{WARNING , std::make_tuple(stderr, "\e[93m[WARNING] \e[0m", "\n")},
	{MESSAGE , std::make_tuple(stderr, "\e[92m[MESSAGE] \e[0m", "\n")},
	{STDERR , std::make_tuple(stderr, " ", "\n")},
	{STDOUT , std::make_tuple(stdout, " ", "\n")},
	{DEBUG , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG1 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG2 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG3 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG4 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG5 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG6 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")},
	{DEBUG7 , std::make_tuple(stderr, "\e[96m[ DEBUG ] \e[0m", "\n")}
      };
      return  printMap.at(level);
    }

    template<int level>
    static inline void log(char const *fmt, ...){
      if(level<=maxLogLevel){
	const auto currentLevelInfo = getLogLevelInfo(level);
	auto stream = std::get<0>(currentLevelInfo);
	auto prefix = std::get<1>(currentLevelInfo);
	auto suffix = std::get<2>(currentLevelInfo);
	va_list args;
	va_start(args, fmt);
	fprintf(stream, "%s", prefix);
	vfprintf(stream, fmt, args);
	fprintf(stream, "%s", suffix);
	va_end(args);
      }
    }

    template<int level>
    static inline void log(const std::string &msg){
      log<level>("%s", msg.c_str());
    }

  }
}
#endif
