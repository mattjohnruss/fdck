#pragma once

#include <problem.h>

namespace mjrfd
{
    class FiniteDifferenceProblem : public Problem
    {
    public:
        FiniteDifferenceProblem(const unsigned n_var,
                                const unsigned n_dof_per_var,
                                const double dt = 0.0);

        virtual ~FiniteDifferenceProblem();

        template<class T>
        const double stencil_1_central(const unsigned t,
                                       const T variable,
                                       const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            return -0.5*u(t,v,i-1) + 0.5*u(t,v,i+1);
        }

        template<class T>
        const double stencil_1_forward(const unsigned t,
                                       const T variable,
                                       const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            return -1.5*u(t,v,i) + 2.0*u(t,v,i+1) - 0.5*u(t,v,i+2);
        }

        template<class T>
        const double stencil_1_backward(const unsigned t,
                                        const T variable,
                                        const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            return 0.5*u(t,v,i-2) - 2.0*u(t,v,i-1) + 1.5*u(t,v,i);
        }

        template<class T>
        const double stencil_2_central(const unsigned t,
                                       const T variable,
                                       const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            return 1.0*u(t,v,i-1) - 2.0*u(t,v,i) + 1.0*u(t,v,i+1);
        }


    };
}
