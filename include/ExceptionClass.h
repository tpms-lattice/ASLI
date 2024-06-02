/* ==========================================================================
 *  OmegaException class. Taken from P. Muldoon "Exceptionally Bad: The
 *  Misuse of Exceptions in C++ & How to Do Better" presented at CppCon 2023
 * ==========================================================================*/

#ifndef OMEGA_EXCEPTION_H
#define OMEGA_EXCEPTION_H

#include <cassert>
#include <future>
#include <iostream>
#include <memory>
#include <string>
#include <variant>
#include <vector>
#include <stdexcept>
#include <source_location>
#include <map>
//#include <stacktrace>

template<typename DATA_T>
class OmegaException {
  public:
  OmegaException(std::string str, DATA_T data, const std::source_location& loc =
    std::source_location::current())//, std::stacktrace trace = std::stacktrace::current()
    : err_str_{std::move(str)}, user_data_{std::move(data)}, location_{loc}{}//, backtrace_{trace}

  std::string& what() { return err_str_; }
  const std::string& what() const noexcept { return err_str_; }
  const std::source_location& where() const noexcept { return location_; }
  //const std::stacktrace& stack() const noexcept { return backtrace_; }
  DATA_T& data(){ return user_data_;}
  const DATA_T& data() const noexcept { return user_data_; }

  private:
  std::string err_str_;
  DATA_T user_data_;
  const std::source_location location_;
  //const std::stacktrace backtrace_;
};

inline std::ostream& operator << (std::ostream& os, const std::source_location& location) {
    os << location.file_name() << "("
       << location.line() << ":"
       << location.column() << "), function `"
       << location.function_name() << "`";
    return os;
}

//inline std::ostream& operator << (std::ostream& os, const std::stacktrace& backtrace) {
//    for(auto iter = backtrace.begin(); iter != (backtrace.end()-3); ++iter) {
//        os << iter->source_file() << "(" << iter->source_line() 
//        << ") : " << iter->description() << "\n";
//    }
//    return os;
//}

using ExceptionError = OmegaException<void*>;

#endif