#ifndef ENUMS_H
#define ENUMS_H

#include <iostream>

enum class error_io
{
  success,
  error_read,
  error_memory,
  incorrect_arguments,
  division_by_zero,
  unsupported_mu,
  unsupported_C_rho,
  unsupported_programm_mode,
  big_relative_error,
  not_stable,
  COUNT
};

enum class equation_of_state_type_t
{
  linear,
  exponential,
  COUNT
};

enum class programm_mode_t
{
  debugging_test,
  COUNT
};

inline void print_error_io (error_io error)
{
  switch (error)
    {
    case error_io::success:
      break;
    case error_io::error_read:
      printf ("ERROR: Error read\n");
      break;
    case error_io::error_memory:
      printf ("ERROR: Error memory\n");
      break;
    case error_io::incorrect_arguments:
      printf ("ERROR: Incorrect arguments\n");
      break;
    case error_io::division_by_zero:
      printf ("ERROR: Division by zero\n");
      break;
    case error_io::unsupported_mu:
      printf ("WARNING: Unsupported mu\n");
      break;
    case error_io::unsupported_C_rho:
      printf ("WARNING: Unsupported C_rho\n");
      break;
    case error_io::unsupported_programm_mode:
      printf ("ERROR: Unsupported programm mode\n");
      break;
    case error_io::big_relative_error:
      printf ("WARNING: Big relative error\n");
      break;
    case error_io::not_stable:
      printf ("WARNING: Not stable\n");
      break;
    default:
      printf ("ERROR: Unknown error\n");
      break;
    }
}

#endif // ENUMS_H
