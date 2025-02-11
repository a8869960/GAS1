#include <iostream>
#include <fenv.h>

#include "arguments.h"
#include "grid.h"
#include "solver.h"

int main (int argc, char *argv[])
{
  feenableexcept (FE_ALL_EXCEPT ^ FE_INEXACT);
  arguments_t args;

  // Parse command line
  const error_io parse_res = args.read (argc, argv);
  if (parse_res != error_io::success)
    {
      print_error_io (parse_res);
      printf ("%s%s\n", argv[0], args.error_message);
      return -1;
    }

  const grid_t grid (args.tau, args.h);
  solver_t solver (grid, args);
  const error_io status = solver.solve_systems ();
  print_error_io (status);


//  printf ("%s: task = %d, x = %.3lf, X = %.3lf, t = %.3lf, T = %.3lf, tau = %10.3e, h = %10.3e, mu = %.3lf, C_rho = %.3lf, (r_u, r_rho) = (%10.7e, %10.7e)\n",
//          argv[0], args.task, grid.get_x (), grid.get_X (), grid.get_t (), grid.get_T (), args.tau, args.h, args.mu, args.C_rho,
//          solver.compute_error_u (), solver.compute_error_rho ());

  printf ("%s: tau = %.0e, h = %.0e, C_rho = %.0lf, mu = %.4lf, (r_u, r_rho) = (%10.3e, %10.3e)\n",
          argv[0], args.tau, args.h, args.C_rho, args.mu,
          solver.compute_error_u (), solver.compute_error_rho ());

  return 0;
}
