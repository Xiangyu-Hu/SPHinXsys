/**
 * @file 	exception.h
 * @details A toolkit for exception functionality. 
 * @author	Chi Zhang and Xiangyu Hu.
 */
#pragma once

#include "base_data_package.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <cassert>

namespace SPH {
	/**
 	 * @class Exception 
 	 * @brief A class for basic exception functionality.
 	 * @details Exception class manages the concatenation of error messages from all the 
 	 * 			derived classes. When creating new exceptions, remember to call addMessage()
 	 * 			as shown above if the exception class does have any error message.
 	 */
	class Exception : public std::exception {
	private:
    	std::string msg_;	/**< A user set message for the exception. */ 
    	std::string file_;	/**< File in which the error occurred. */
    	int line_;			/**< Line number at which the error occurred. */
public:
    /** 
     * @brief Default Constructor 
     * @details This constructor is for backward compatibility. Use the constructor
     * 			taking file, line, func.                                                  
     */
    Exception(const std::string &aMsg="", const std::string &aFile="", int aLine=-1);

    /** 
     * @brief Call this constructor from derived classes to add file, line and 
     * 			function information to the error message. Use this when throwing
     * 			Derived classes. Use OPENSIM_THROW_<> macros at throw sites.              
     */
    Exception(const std::string& file, size_t line, const std::string& func);
    /** 
     * @brief Use this when you want to throw an Exception and also provide a message.
     */ 
    Exception(const std::string& file, size_t line,
              const std::string& func, const std::string& msg);
    /**
     * @ brief Destructor.
     */
    virtual ~Exception() throw() {}

	protected:
    	/** Add to the error message that will be returned for the exception.     */
    	void addMessage(const std::string& msg);

	private:
    	void setNull();

	public:
    	void setMessage(const std::string &aMsg);
    	const char* getMessage() const;

    	virtual void print(std::ostream &aOut) const;

    	const char* what() const noexcept override;
	};

	class InvalidArgument : public Exception {
	public:
    	InvalidArgument(const std::string& file,
                    size_t line,
                    const std::string& func,
                    const std::string& msg = "") :
        	Exception(file, line, func) {
        	std::string mesg = "Invalid Argument. " + msg;

        	addMessage(mesg);
    	}
	};

	class InvalidCall : public Exception {
	public:
    	InvalidCall(const std::string& file,
                size_t line,
                const std::string& func,
                const std::string& msg = "") :
        	Exception(file, line, func) {
        	std::string mesg = "Invalid Call. " + msg;

        	addMessage(mesg);
    	}
	};

	class InvalidTemplateArgument : public Exception {
	public:
    	InvalidTemplateArgument(const std::string& file,
                            size_t line,
                            const std::string& func,
                            const std::string& msg) :
        	Exception(file, line, func) {
        	std::string mesg = "Invalid Template argument. " + msg;

        	addMessage(mesg);
    	}
	};

	class IndexOutOfRange : public Exception {
	public:
    	IndexOutOfRange(const std::string& file,
                    size_t line,
                    const std::string& func,
                    size_t index,
                    size_t min, 
                    size_t max) :
        	Exception(file, line, func) {
        	std::string msg = "min = " + std::to_string(min);
        	msg += " max = " + std::to_string(max);
        	msg += " index = " + std::to_string(index);

        	addMessage(msg);
    	}
	};

	class KeyNotFound : public Exception {
	public:
    	KeyNotFound(const std::string& file,
                size_t line,
                const std::string& func,
                const std::string& key) :
        	Exception(file, line, func) {
        	std::string msg = "Key '" + key + "' not found.";

        	addMessage(msg);
    	}
	};
}/** end of namespace */