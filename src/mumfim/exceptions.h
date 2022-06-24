#ifndef MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#define MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#include <stdexcept>
#include <petsc.h>
namespace mumfim
{
  class mumfim_error : public std::runtime_error
  {
    public:
    explicit mumfim_error(const std::string & err) : std::runtime_error{err} {}
    explicit mumfim_error(const char * err) : std::runtime_error{err} {}
  };
  class petsc_error : public std::exception
  {
    public:
    explicit petsc_error(PetscErrorCode error_code)
    {
      const char * text;
      char * specific;
      PetscErrorMessage(error_code, &text, &specific);
      std::stringstream error_message_stream;
      error_message_stream << "Petsc Error " << error_code << "\n"
                    << text << "\n"
                    << specific;
      error_message_ = error_message_stream.str();
    }
    [[ nodiscard]]
    const char * what() const noexcept override {
      return error_message_.c_str();
    }
    std::string error_message_;
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
