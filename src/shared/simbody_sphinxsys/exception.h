/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	exception.h
 * @details A toolkit for exception functionality.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef EXCEPTION_SIMBODY_H
#define EXCEPTION_SIMBODY_H

#include "base_data_type_package.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

namespace SPH
{
/**
 * @class Exception
 * @brief A class for basic exception functionality.
 * @details Exception class manages the concatenation of error messages from all the
 * 			derived classes. When creating new exceptions, remember to call addMessage()
 * 			as shown above if the exception class does have any error message.
 */
class Exception : public std::exception
{
  private:
    std::string msg_;  /**< A user set message for the exception. */
    std::string file_; /**< File in which the error occurred. */
    int line_;         /**< Line number at which the error occurred. */
  public:
    /**
     * @brief Default Constructor
     * @details This constructor is for backward compatibility. Use the constructor
     * 			taking file, line, func.
     */
    Exception(const std::string &aMsg = "", const std::string &aFile = "", int aLine = -1);

    /**
     * @brief Call this constructor from derived classes to add file, line and
     * 			function information to the error message. Use this when throwing
     * 			Derived classes. Use OPENSIM_THROW_<> macros at throw sites.
     */
    Exception(const std::string &file, size_t line, const std::string &func);
    /**
     * @brief Use this when you want to throw an Exception and also provide a message.
     */
    Exception(const std::string &file, size_t line,
              const std::string &func, const std::string &msg);
    /**
     * @ brief Destructor.
     */
    virtual ~Exception() throw() {}

  protected:
    /** Add to the error message that will be returned for the exception.     */
    void addMessage(const std::string &msg);

  private:
    void setNull();

  public:
    void setMessage(const std::string &aMsg);
    const char *getMessage() const;

    virtual void print(std::ostream &aOut) const;

    const char *what() const noexcept override;
};

class InvalidArgument : public Exception
{
  public:
    InvalidArgument(const std::string &file,
                    size_t line,
                    const std::string &func,
                    const std::string &msg = "") : Exception(file, line, func)
    {
        std::string mesg = "Invalid Argument. " + msg;

        addMessage(mesg);
    }
};

class InvalidCall : public Exception
{
  public:
    InvalidCall(const std::string &file,
                size_t line,
                const std::string &func,
                const std::string &msg = "") : Exception(file, line, func)
    {
        std::string mesg = "Invalid Call. " + msg;

        addMessage(mesg);
    }
};

class InvalidTemplateArgument : public Exception
{
  public:
    InvalidTemplateArgument(const std::string &file,
                            size_t line,
                            const std::string &func,
                            const std::string &msg) : Exception(file, line, func)
    {
        std::string mesg = "Invalid Template argument. " + msg;

        addMessage(mesg);
    }
};

class IndexOutOfRange : public Exception
{
  public:
    IndexOutOfRange(const std::string &file,
                    size_t line,
                    const std::string &func,
                    size_t index,
                    size_t min,
                    size_t max) : Exception(file, line, func)
    {
        std::string msg = "min = " + std::to_string(min);
        msg += " max = " + std::to_string(max);
        msg += " index = " + std::to_string(index);

        addMessage(msg);
    }
};

class KeyNotFound : public Exception
{
  public:
    KeyNotFound(const std::string &file,
                size_t line,
                const std::string &func,
                const std::string &key) : Exception(file, line, func)
    {
        std::string msg = "Key '" + key + "' not found.";

        addMessage(msg);
    }
};
} // namespace SPH
#endif // EXCEPTION_SIMBODY_H