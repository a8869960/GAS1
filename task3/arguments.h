#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "math.h"
#include <cassert>
#include <iostream>

#include "enums.h"
#include "functions.h"
#include "initial_conditions.h"

#define EPS 1e-16

struct arguments_t
{
  const int task = 5;

  double mu = 1e-3;
  double C_rho = 0.;
  double tau = 1.;
  double h = 1.;

  unsigned int k;
  double eps;

  bool DEBUG;

  pressure_function_t p;
  initial_conditions_t functions;

  const char *error_message = " <mu> <C_rho> <tau> <h> <1 or 2> <p_type> <is_debug> <eps>";

  void set_all_functions_for_debugging_test ();

  error_io is_mu_correct () const;
  error_io is_C_rho_correct () const;
  bool is_correct () const;
  error_io read (const int argc, char *argv[]);

};

#endif // ARGUMENTS_H
