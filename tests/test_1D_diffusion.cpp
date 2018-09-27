#include <advection_diffusion_reaction_problem.h>

#include <catch2/catch.hpp>

#include <iomanip>

enum Variable
{
    c = 0,
};

using namespace mjrfd;

struct AdvectionDiffusionParams
{
    double Pe;
    double a;
    double b;
};

#define PRINT_VAR(var, width) \
    #var" = " << std::setw(width) << var

class AdvectionDiffusionTestProblem : public AdvectionDiffusionReactionProblem
{
public:
    AdvectionDiffusionTestProblem(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(1, n_node)
    {
        enable_bc(Boundary::Left,  { c });
        enable_bc(Boundary::Right, { c });

        enable_spatial_terms({ c });

        set_variable_names({ "c" });

        Max_residual = 1.0e-14;

        u(c, 0) = 1.0;

        for(unsigned i = 1; i < n_node_; ++i)
        {
            //u(c, i) = 1.0 - this->x(i);
            u(c, i) = 0.0;
        }
    }

    ~AdvectionDiffusionTestProblem()
    {
    }

    double absolute_error_with_exact_solution(const unsigned i) const
    {
        const double x = this->x(i);
        std::vector<double> exact_u(1);
        exact_solution(time(), x, exact_u);

        return std::abs(u(c,i) - exact_u[0]);
    }

    double relative_error_with_exact_solution(const unsigned i) const
    {
        const double x = this->x(i);
        std::vector<double> exact_u(1);
        exact_solution(time(), x, exact_u);

        const double u = this->u(c,i);

        //std::cout << std::setprecision(16) << "std::abs(u) = " << std::setw(22) << std::abs(u) << ", std::abs(exact_u[0]) = " << std::setw(22) << std::abs(exact_u[0]) << '\n';
        //std::cout << std::setprecision(16) << PRINT_VAR(std::abs(u), 22) << ", " << PRINT_VAR(std::abs(exact_u[0]), 22) << '\n';
        if(std::abs(u) <= 1.0e-16 && std::abs(exact_u[0]) <= 1.0e-16)
        {
            return 0.0;
        }

        return std::abs((u - exact_u[0])/exact_u[0]);
    }

    AdvectionDiffusionParams p;

private:
    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c] = 1.0; a2[c] = 0.0; a3[c] = 1.0;
        }
        else if(b == Boundary::Right)
        {
            a1[c] = 1.0; a2[c] = 0.0; a3[c] = 0.0;
        }
    }

    void get_d(const unsigned,
               const unsigned,
               std::vector<double> &d) const override
    {
        d[c] = 1.0;
    }

    void get_v(const unsigned,
               const unsigned,
               double &v) const override
    {
        v = p.Pe;
    }

    void get_r(const std::vector<double> &,
               std::vector<double> &r) const override
    {
        r[c] = 0.0;
    }

    void exact_solution(const double,
                        const double x,
                        std::vector<double> &sol) const override
    {
        if(std::abs(p.Pe) <= 1.0e-16)
            sol[c] = 1.0 - x;
        else
            sol[c] = (std::exp(p.Pe*x) - std::exp(p.Pe))/(1.0 - std::exp(p.Pe));
    }
};

namespace mjrfd
{
    TEST_CASE( "Simple diffusion with Dirichlet BCs", "[diffusion]" )
    {
        const unsigned n_node = 10001;

        AdvectionDiffusionTestProblem problem(n_node);
        problem.p.Pe = 0.2;

        problem.steady_solve();

        for(unsigned i = 0; i < n_node; ++i)
        {
            CHECK( problem.absolute_error_with_exact_solution(i) <= 1.0e-8 );
        }
    }
}