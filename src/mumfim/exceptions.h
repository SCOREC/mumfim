#ifndef MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#define MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#include <stdexcept>
namespace mumfim
{
  class mumfim_error : public std::logic_error
  {
    public:
    explicit mumfim_error(const std::string & err) : std::logic_error{err} {}
    explicit mumfim_error(const char * err) : std::logic_error{err} {}
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
