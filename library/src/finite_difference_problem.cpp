#include <finite_difference_problem.h>

namespace mjrfd
{
    FiniteDifferenceProblem::FiniteDifferenceProblem(const unsigned n_dof,
                                                     const double dt) :
        Problem(n_dof, dt)
    {
    }

    FiniteDifferenceProblem::~FiniteDifferenceProblem()
    {
    }

    const double FiniteDifferenceProblem::stencil_1_central(const unsigned t,
                                                            const unsigned i) const
    {
        return -0.5*u(t,i-1) + 0.5*u(t,i+1);
    }

    const double FiniteDifferenceProblem::stencil_1_forward(const unsigned t,
                                                            const unsigned i) const
    {
        return -1.5*u(t,i) + 2.0*u(t,i+1) - 0.5*u(t,i+2);
    }

    const double FiniteDifferenceProblem::stencil_1_backward(const unsigned t,
                                                             const unsigned i) const
    {
        return 0.5*u(t,i-2) - 2.0*u(t,i-1) + 1.5*u(t,i);
    }

    const double FiniteDifferenceProblem::stencil_2_central(const unsigned t,
                                                            const unsigned i) const
    {
        return 1.0*u(t,i-1) - 2.0*u(t,i) + 1.0*u(t,i+1);
    }
}
