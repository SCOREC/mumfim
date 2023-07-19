#ifndef MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#define MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
#include <petsc.h>

#include <stdexcept>

namespace mumfim
{
  struct mumfim_error : std::runtime_error
  {
    using std::runtime_error::runtime_error;
  };

  struct material_error : mumfim_error
  {
    using mumfim_error::mumfim_error;
  };

  class petsc_error : public mumfim_error
  {
    public:
    explicit petsc_error(PetscErrorCode error_code)
        : mumfim_error("Petsc Error")
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

    [[nodiscard]] const char * what() const noexcept override
    {
      return error_message_.c_str();
    }

    std::string error_message_;
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_EXCEPTIONS_H
