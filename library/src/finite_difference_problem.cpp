#include <finite_difference_problem.h>

namespace mjrfd
{
    FiniteDifferenceProblem::FiniteDifferenceProblem(const unsigned n_var,
                                                     const unsigned n_dof_per_var,
                                                     const double dt) :
        Problem(n_var, n_dof_per_var, dt)
    {
    }

    FiniteDifferenceProblem::~FiniteDifferenceProblem()
    {
    }
}
