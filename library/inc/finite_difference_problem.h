#pragma once

#include <problem.h>

namespace mjrfd
{
    class FiniteDifferenceProblem : public Problem
    {
    public:
        FiniteDifferenceProblem(const unsigned n_dof, const double dt = 0.0);
        virtual ~FiniteDifferenceProblem();

        const double stencil_1_central(const unsigned t, const unsigned i) const;
        const double stencil_1_forward(const unsigned t, const unsigned i) const;
        const double stencil_1_backward(const unsigned t, const unsigned i) const;
        const double stencil_2_central(const unsigned t, const unsigned i) const;
    };
}
