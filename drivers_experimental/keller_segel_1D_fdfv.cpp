#include "problem.h"
#include "stencil.h"
#include "log.h"
#include "config.h"
#include "utilities.h"

#include <tuple>
#include <fstream>
#include <random>
#include <functional>
#include <optional>

#include <fenv.h>

using namespace fdck;

enum Variable
{
    rho_bar = 0,
    c = 1,
};

enum class KellerSegelICs
{
    constant,
    perturbed,
    gaussian,
    exact,
};

struct KellerSegelParameters
{
    bool imposed_c_profile = false;
    double rho_bar_init = 1.0;
    double c_init = 1.0;
    KellerSegelICs initial_conditions = KellerSegelICs::constant;
    double D = 1.0;
    double r = 1.0;
    double chi = 0.0;
    double gamma_rho = 1.0;
    double gamma_c = 1.0;
    unsigned n_interior_cell_1d = 3;
    double L = 1.0;
    double dx = L/n_interior_cell_1d;
};

//namespace exact_solution
//{
    //double m;
    //double omega;
    //const double L = 1;

    //using std::sin;
    //using std::cos;
    //using std::exp;
    //using std::pow;

    //double exact_u(double x, double t, const KellerSegelParameters &)
    //{
        //const double pi = std::acos(-1.0);
        //return pow(sin(pi*omega*t), 2.0)*exp(-m*(x - 0.5)*(x - 0.5));
    //}

    //// NOTE this is actually rho itself, not the cell-averaged rho_bar. Do we
    //// need to cell-average this and rho_force for meaningful comparison with
    //// numerics?
    //double exact_rho(double x, double t, const KellerSegelParameters &p)
    //{
        //return exact_u(x, t, p);
    //}

    //double exact_c(double x, double t, const KellerSegelParameters &p)
    //{
        //return exact_u(x, t, p);
    //}

    //// NOTE above comment on exact_rho applies here as well
    //double force_rho(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //const double term_1 =
            //exp(-2.0*m*pow(x - 0.5, 2.0))*
            //(
             //-2.0*m*(m*pow(1.0 - 2.0*x, 2.0) - 1.0)*p.chi*pow(sin(pi*t*omega), 4.0)
            //);

        //const double term_2 =
            //exp(-m*pow(x - 0.5, 2.0))*sin(pi*t*omega)*
            //(
             //2.0*pi*omega*cos(pi*t*omega) +
             //m*(m*pow(1.0 - 2.0*x, 2.0) - 2.0)*sin(pi*t*omega)
            //);

        //return term_1 + term_2;
    //}

    //double force_c(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //return
            //exp(-m*pow(x - 0.5, 2.0))*sin(pi*t*omega)*
            //(
             //2.0*pi*omega*cos(pi*t*omega) +
             //sin(pi*t*omega)*(m*(m*pow(1.0 - 2.0*x, 2.0) - 2.0) + p.gamma_c - p.gamma_rho)
            //);
    //}
//}

//namespace exact_solution
//{
    //double omega;
    //const double L = 1;

    //using std::sin;
    //using std::cos;
    //using std::exp;
    //using std::pow;

    //// NOTE this is the cell-averaged rho_bar. Let's see if this is the right
    //// thing to do...
    //double exact_rho_bar(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);
        //return pow(sin(pi*omega*t), 2.0)*(1.0/12.0*p.dx*(12*x*(1.0 - x) - p.dx*p.dx));
    //}

    //double exact_c(double x, double t, const KellerSegelParameters &)
    //{
        //const double pi = std::acos(-1.0);
        //return pow(sin(pi*omega*t), 2.0)*x*(1.0 - x);
    //}

    //// NOTE above comment on exact_rho applies here as well
    //double force_rho_bar(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //return
            //1.0/6.0*p.dx*sin(pi*omega*t)*
            //(
             //- pi*(12.0*x*(x - 1.0) + p.dx*p.dx)*omega*cos(pi*omega*t)
             //+ 12.0*sin(pi*omega*t)
             //+ (6.0 + 36.0*x*(x - 1.0) + p.dx*p.dx)*p.chi*pow(sin(pi*omega*t), 3.0)
            //);
    //}

    //double force_c(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //return
            //1.0/12.0*pow(sin(pi*omega*t), 2.0)*
            //(
             //24.0
             //- 12.0*(x - 1.0)*x*p.gamma_c
             //+ 12.0*(x - 1.0)*x*p.gamma_rho*p.dx
             //+ p.gamma_rho*pow(p.dx, 3.0)
            //)
            //- pi*(x - 1.0)*x*omega*sin(2.0*pi*omega*t);
    //}
//}

//namespace exact_solution
//{
    //double omega;
    //const double L = 1;

    //using std::sin;
    //using std::cos;
    //using std::exp;
    //using std::pow;

    //// NOTE this is the cell-averaged rho_bar. Let's see if this is the right
    //// thing to do...
    //double exact_rho_bar(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);
        //return 1.0/240.0*p.dx*pow(sin(pi*t*omega), 2.0)*
            //(
             //240.0*x*x*(x - 1.0)*(x - 1.0)
             //+ 20.0*(1.0 + 6.0*x*(x - 1.0))*pow(p.dx, 2.0)
             //+ 3.0*pow(p.dx, 4.0)
            //);
    //}

    //double exact_c(double x, double t, const KellerSegelParameters &)
    //{
        //const double pi = std::acos(-1.0);
        //return sin(pi*omega*t)*x*x*(1.0 - x)*(1.0 - x);
    //}

    //// NOTE above comment on exact_rho applies here as well
    //double force_rho_bar(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //// TODO clean up this mess from Mathematica
        //return
            //(p.dx*sin(omega*pi*t)*(omega*pi*cos(omega*pi*t)*
            //(20*(1 + 6*(-1 + x)*x)*pow(p.dx,2) + 3*pow(p.dx,4) + 
             //240*pow(-1 + x,2)*pow(x,2)) + 
            //p.chi*(20*(1 + 6*(-1 + x)*x*(3 + 10*(-1 + x)*x))*pow(p.dx,2) + 
                   //3*(1 + 6*(-1 + x)*x)*pow(p.dx,4) + 
                   //240*(3 + 14*(-1 + x)*x)*pow(-1 + x,2)*pow(x,2))*
            //pow(sin(omega*pi*t),3) - 
            //120*(2 + 12*(-1 + x)*x + pow(p.dx,2))*sin(omega*pi*t)))/120.;
    //}

    //double force_c(double x, double t, const KellerSegelParameters &p)
    //{
        //const double pi = std::acos(-1.0);

        //return
            //((240*(-2 + (-1 + x)*x*(-12 + p.gamma_c*(-1 + x)*x)) - 20*p.gamma_rho*(1 + 6*(-1 + x)*x)*pow(p.dx,3) - 3*p.gamma_rho*pow(p.dx,5) - 240*p.dx*p.gamma_rho*pow(-1 + x,2)*pow(x,2))*
             //pow(sin(omega*pi*t),2))/240. + omega*pi*pow(-1 + x,2)*pow(x,2)*sin(2*omega*pi*t);
    //}
//}

namespace exact_solution
{
    double omega;
    const double L = 1;

    using std::sin;
    using std::cos;
    using std::exp;
    using std::pow;

    double exact_rho(double x, double t, const KellerSegelParameters &p)
    {
        const double pi = std::acos(-1.0);
        return sin(pi*omega*t)*x*x*(1.0 - x)*(1.0 - x);
    }

    double exact_c(double x, double t, const KellerSegelParameters &)
    {
        const double pi = std::acos(-1.0);
        return sin(pi*omega*t)*x*x*(1.0 - x)*(1.0 - x);
    }

    double force_rho(double x, double t, const KellerSegelParameters &p)
    {
        const double pi = std::acos(-1.0);

        return
            -((2 + 12*(-1 + x)*x + pow(p.dx,2))*pow(sin(omega*pi*t),2)) + 
            (p.chi*pow(p.dx,-1)*((-1 - p.dx + 2*x)*pow(2 + p.dx - 2*x,3)*pow(-p.dx + 2*x,3) + 
                             (-1 + p.dx + 2*x)*pow(-2 + p.dx + 2*x,3)*pow(p.dx + 2*x,3))*
             pow(sin(omega*pi*t),4))/32. + 
            (omega*pi*(20*(1 + 6*(-1 + x)*x)*pow(p.dx,2) + 3*pow(p.dx,4) + 
                       240*pow(-1 + x,2)*pow(x,2))*sin(2*omega*pi*t))/240.;
    }

    double force_c(double x, double t, const KellerSegelParameters &p)
    {
        const double pi = std::acos(-1.0);

        return
            sin(omega*pi*t)*
            (
             2*omega*pi*cos(omega*pi*t)*pow(-1 + x,2)*pow(x,2)
             + (-2 + (-1 + x)*x*(-12 + (p.gamma_c - p.gamma_rho)*(-1 + x)*x))*sin(omega*pi*t)
            );
    }
}

// Hybrid FVFD method for Keller-Segel based on Chertock et al (2017)
class KellerSegelProblem1D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // Square grid of cells, with one layer of ghost cells added outside each
    // boundary.
    // We need no aux dofs but we need 3 history values for SSP R-K timestepping.
    // TODO check if dx is correctly defined
    //KellerSegelProblem1D(const unsigned n_interior_cell_1d, double L) :
    KellerSegelProblem1D(KellerSegelParameters p_in) :
        Problem(2, p_in.n_interior_cell_1d + 2, 0, 3),
        p{p_in},
        force_rho{},
        force_c{},
        exact_rho{},
        exact_c{}
    {
        p.dx = p.L/static_cast<double>(p.n_interior_cell_1d);

        FDCK_DEBUG("p.n_interior_cell_1d = {}, p.dx = {}", p.n_interior_cell_1d, p.dx);

        FDCK_INFO("Creating grid with {} interior cells ({} including ghost cells).", p.n_interior_cell_1d, p.n_interior_cell_1d+2);
        Max_residual = 1.0e-4;
        // TODO investigate why setting this to 1 causes (possibly) no
        // iterations... at the very least the logging is confusing
        Max_newton_iterations = 2;
    }

    ~KellerSegelProblem1D() override
    {
    }

    void output(std::ostream &out) const override
    {
        double time = this->time();

        if(exact_rho && exact_c)
        {
            out << "t x rho_bar c exact_rho exact_c\n";
            for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
            {
                double x = this->x(i);

                out << time << ' '
                    << x << ' '
                    << u(rho_bar, index(i)) << ' '
                    << u(c, index(i)) << ' '
                    << (*exact_rho)(x, time, p) << ' '
                    << (*exact_c)(x, time, p) << '\n';
            }
        }
        else
        {
            out << "t x rho_bar c\n";
            for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
            {
                double x = this->x(i);

                out << time << ' '
                    << x << ' '
                    << u(rho_bar, index(i)) << ' '
                    << u(c, index(i)) << '\n';
            }
        }
    }

    void set_initial_conditions()
    {
        FDCK_TRACE("set_initial_conditions()");

        clear_solution();

        switch(p.initial_conditions)
        {
            case KellerSegelICs::constant:
                for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
                {
                    u(rho_bar, index(i)) = p.rho_bar_init;
                    u(c, index(i)) = p.c_init;
                }
                break;
            case KellerSegelICs::perturbed:
                {
                    std::random_device rd{};
                    std::mt19937 gen{rd()};
                    std::normal_distribution<double> dist{0.0, 0.001};

                    for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
                    {
                        u(rho_bar, index(i)) = p.rho_bar_init + dist(gen);
                        u(c, index(i)) = p.c_init + dist(gen);
                    }
                }
                break;
            case KellerSegelICs::gaussian:
                for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
                {
                    const double x_c = this->x(i) - 0.5;
                    u(rho_bar, index(i)) = 1000.0*std::exp(-100.0*x_c*x_c);
                    u(c, index(i)) = 500.0*std::exp(-50.0*x_c*x_c);
                }
                break;
            case KellerSegelICs::exact:
                {
                    assert(exact_rho && exact_c);
                    const double time = this->time();
                    for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
                    {
                        const double x = this->x(i);
                        u(rho_bar, index(i)) = (*exact_rho)(x, time, p);
                        u(c, index(i)) = (*exact_c)(x, time, p);
                    }
                }
                break;
            default:
                FDCK_FATAL("Invalid initial conditions");
                std::exit(1);
        }

        if(p.imposed_c_profile == true)
        {
            for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
            {
                u(c, index(i)) = c_value(x(i));
            }
        }
    }

    double x(unsigned i) const
    {
        FDCK_DEBUG("i = {}, p.dx = {}, x = {}",
                    i,
                    p.dx,
                    (static_cast<double>(i) - 0.5)*p.dx);

        return (static_cast<double>(i) - 0.5)*p.dx;
    }

    void ssp_rk_3_timestep(double dt)
    {
        // backup the current time for working around the time += dt update
        const double current_time = time();

        // take a forward Euler step to give w^(1)
        unsteady_solve(dt);

        // reset time before next solve
        time() = current_time;

        // take another forward Euler step to give the [ ... ] brackets in w^(2)
        unsteady_solve(dt);

        // calculate w^(2)
        for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
        {
            u(rho_bar, index(i)) =
                0.75*u(2, rho_bar, index(i)) + 0.25*u(rho_bar, index(i));
        }

        // reset time before next solve
        time() = current_time;

        // take another forward Euler step to give the [ ... ] brackets in w^(3)
        unsteady_solve(dt);

        // calculate w^(3)
        for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
        {
            u(rho_bar, index(i)) =
                (1.0/3.0)*u(3, rho_bar, index(i)) + (2.0/3.0)*u(rho_bar, index(i));
        }

        // update the time manually
        time() = current_time + dt;
    }

    KellerSegelParameters p;

    typedef std::function<double(double, double, const KellerSegelParameters &)> exact_solution_function;

    typedef std::function<double(double, double, const KellerSegelParameters &)> force_function;

    std::optional<force_function> force_rho;
    std::optional<force_function> force_c;

    std::optional<exact_solution_function> exact_rho;
    std::optional<exact_solution_function> exact_c;

private:
    enum class Face
    {
        East,
        West,
    };

    // TODO double check all p.dx/p.dx/dt denominators and multiply through as in
    // other problems

    // Returns the node number corresponding to the grid location i
    // Indexing scheme is as follows.
    // Interior cells: i = 1, ..., p.n_interior_cell_1d,
    // Ghost cells: i = 0, i = p.n_interior_cell_1d + 1
    unsigned index(unsigned i) const
    {
        FDCK_TRACE("index(i = {})", i);

        assert(i <= p.n_interior_cell_1d + 1);

        // The grid location and node numbering coincide in 1d
        return i;
    }
    
    // Equation (2.3)(a)
    // Calculates F_{i+1/2}
    double f_p(int i) const
    {
        FDCK_TRACE("f_p(i = {})", i);
        double u_p = u_p_at_midpoint(i);
        double drho_dx_p = drho_dx_p_at_midpoint(i);
        return p.chi*rho_point_value_i_p(i)*u_p - p.D*drho_dx_p;
    }

    // Calculates F_{i-1/2}
    double f_m(unsigned i) const
    {
        FDCK_TRACE("f_m(i = {})", i);
        // for the next two, we can call the plus (_p) functions with args i-1,
        // j-1 as it gives the same result as a manual impl of the minus functions
        //return f_p(i - 1, j);
        double u_m = u_p_at_midpoint(i-1);
        double drho_dx_m = drho_dx_p_at_midpoint(i-1);
        // rho_point_value_i_m must be written explicitly because the upwinding
        // should be in the opposite direction than the plus case
        return p.chi*rho_point_value_i_m(i)*u_m - p.D*drho_dx_m;
    }

    // Equation (2.4)(a)
    double drho_dx_p_at_midpoint(unsigned i) const
    {
        FDCK_TRACE("drho_dx_p_at_midpoint(i = {})", i);
        double drho_dx_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            drho_dx_p += w*u(1, rho_bar, index(i+k))/p.dx;
        }
        return drho_dx_p;
    }

    // Equation (2.4)(c)
    double u_p_at_midpoint(unsigned i) const
    {
        FDCK_TRACE("u_p_at_midpoint(i = {})", i);
        double u_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            u_p += w*u(1, c, index(i+k))/p.dx;
        }
        return u_p;
    }

    // Equation (2.5)(a)
    double rho_point_value_i_p(unsigned i) const
    {
        FDCK_TRACE("rho_point_value_i_p(i = {})", i);
        double u_p = u_p_at_midpoint(i);

        if(u_p > 0.0)
        {
            //FDCK_DEBUG("i = {}, u_p > 0.0 - upwinding", i);
            return rho_point_value_at_face(i, Face::East);
        }
        else
        {
            //FDCK_DEBUG("i = {}, u_p <= 0.0 - NOT upwinding", i);
            return rho_point_value_at_face(i+1, Face::West);
        }
    }

    double rho_point_value_i_m(unsigned i) const
    {
        FDCK_TRACE("rho_point_value_i_m(i = {})", i);
        double u_m = u_p_at_midpoint(i-1);

        if(u_m > 0.0)
        {
            //FDCK_DEBUG("i = {}, u_m > 0.0 - upwinding", i);
            return rho_point_value_at_face(i-1, Face::East);
        }
        else
        {
            //FDCK_DEBUG("i = {}, u_m <= 0.0 - NOT upwinding", i);
            return rho_point_value_at_face(i, Face::West);
        }
    }

    // Equation (2.7)
    // We have shifted some of the indices here but it should be equivalent to
    // the paper 
    double rho_point_value_at_face(unsigned i, Face f) const
    {
        FDCK_TRACE("rho_point_value_at_face(i = {})", i);
        double result = u(1, rho_bar, index(i));

        switch(f)
        {
            case Face::East:
                result += 0.5*p.dx*drho_dx(i);
                break;
            case Face::West:
                result -= 0.5*p.dx*drho_dx(i);
                break;
        }

        return result;
    }

    // Equation (2.8)(a)
    // TODO What exactly does the +- mean in the conditions here? At the moment
    // we are assuming it means that both the + and - quantities must be >= 0
    // separately
    double drho_dx(unsigned i) const
    {
        FDCK_TRACE("drho_dx(i = {})", i);
        double drho_bar_central = 0.0;
        for(auto [k, w] : stencil::central_1::weights)
        {
            drho_bar_central += w*u(1, rho_bar, index(i+k));
        }

        // Central stencil has weights +-0.5, so here we only need another
        // +-0.5 factor to give the +-0.25 required - see paper
        const double test_p = u(1, rho_bar, index(i)) + 0.5*drho_bar_central;
        const double test_m = u(1, rho_bar, index(i)) - 0.5*drho_bar_central;

        if(test_p >= 0.0 && test_m >= 0.0)
        {
            //FDCK_DEBUG("Both tests positive - using central difference");
            return drho_bar_central/p.dx;
        }
        else
        {
            //FDCK_DEBUG("One or more test failed - using minmod");
            double drho_bar_forward = 0.0;
            for(auto [k, w] : stencil::first_order::forward_1::weights)
            {
                drho_bar_forward += w*u(1, rho_bar, index(i+k));
            }
            double drho_bar_backward = 0.0;
            for(auto [k, w] : stencil::first_order::backward_1::weights)
            {
                drho_bar_backward += w*u(1, rho_bar, index(i+k));
            }
            return utilities::minmod<double>({ 2.0*drho_bar_forward/p.dx,
                                               drho_bar_central/p.dx,
                                               2.0*drho_bar_backward/p.dx });
        }
    }

    double fc_p(unsigned i) const
    {
        FDCK_TRACE("fc_p(i = {})", i);
        double u_p = u_p_at_midpoint(i);
        return -u_p;
    }

    double fc_m(unsigned i) const
    {
        FDCK_TRACE("fc_m(i = {})", i);
        double u_m = u_p_at_midpoint(i-1);
        return -u_m;
    }

    double extrapolate_onto_boundary_p(Variable var, unsigned i) const
    {
        const double quadratic_extrapolation =
            + 3.0*u(1, var, index(i+1))
            - 3.0*u(1, var, index(i+2))
            + 1.0*u(1, var, index(i+3));

        if (quadratic_extrapolation >= 0.0)
        {
            return quadratic_extrapolation;
        }
        else
        {
            return u(1, var, index(i+1));
        }
    }

    double extrapolate_onto_boundary_m(Variable var, unsigned i) const
    {
        const double quadratic_extrapolation =
            + 3.0*u(1, var, index(i-1))
            - 3.0*u(1, var, index(i-2))
            + 1.0*u(1, var, index(i-3));

        if (quadratic_extrapolation >= 0.0)
        {
            return quadratic_extrapolation;
        }
        else
        {
            return u(1, var, index(i-1));
        }
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        // Set the values in the ghost cells using quadratic interpolation
        // of the previous values in the neighbouring interior cells
        {
            // Left boundary
            unsigned i = 0;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}) [left]", i);

            unsigned idx = rho_bar*n_dof_per_var_ + index(i);
            residual(idx) += u(rho_bar, index(i)) - extrapolate_onto_boundary_p(rho_bar, i);
                //u(rho_bar, index(i)) - (u(1, rho_bar, index(i+1)));
                //u(rho_bar, index(i)) - (3.0*u(1, rho_bar, index(i+1)) - 3.0*u(1, rho_bar, index(i+2)) + 1.0*u(1, rho_bar, index(i+3)));

            idx = c*n_dof_per_var_ + index(i);
            residual(idx) += u(c, index(i)) - extrapolate_onto_boundary_p(c, i);
                //u(c, index(i)) - (u(1, c, index(i+1)));
                //u(c, index(i)) - (3.0*u(1, c, index(i+1)) - 3.0*u(1, c, index(i+2)) + 1.0*u(1, c, index(i+3)));

            // Right boundary
            i = p.n_interior_cell_1d+1;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}) [right]", i);

            idx = rho_bar*n_dof_per_var_ + index(i);
            residual(idx) += u(rho_bar, index(i)) - extrapolate_onto_boundary_m(rho_bar, i);
                //u(rho_bar, index(i)) - (u(1, rho_bar, index(i-1)));
                //u(rho_bar, index(i)) - (3.0*u(1, rho_bar, index(i-1)) - 3.0*u(1, rho_bar, index(i-2)) + 1.0*u(1, rho_bar, index(i-3)));

            idx = c*n_dof_per_var_ + index(i);
            residual(idx) += u(c, index(i)) - extrapolate_onto_boundary_m(c, i);
                //u(c, index(i)) - (u(1, c, index(i-1)));
                //u(c, index(i)) - (3.0*u(1, c, index(i-1)) - 3.0*u(1, c, index(i-2)) + 1.0*u(1, c, index(i-3)));
        }

        // Loop over the interior cells
        for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
        {
            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}) [interior]", i);

            // rho_bar
            unsigned idx = rho_bar*n_dof_per_var_ + index(i);

            // Time derivative
            residual(idx) +=
                1.0/dt_*(u(rho_bar, index(i)) - u(1, rho_bar, index(i)));

            // Storage for the fluxes
            double local_f_m = 0.0;
            double local_f_p = 0.0;

            // Unless we're in a boundary cell, calculate the flux. If
            // we're in a boundary cell, the flux at the appropriate
            // boundary is left as 0.0 to implement the no-flux BCs
            if(i != 1)
                local_f_m = f_m(i);
            if(i != p.n_interior_cell_1d)
                local_f_p = f_p(i);

            residual(idx) -= -(local_f_p - local_f_m)/p.dx;

            // Logistic growth term
            residual(idx) -=
                p.r*u(1, rho_bar, index(i))*(1.0 - u(1, rho_bar, index(i)));

            // Forcing (for testing exact solutions)
            // TODO does the time() here work with the SSP RK timestepper,
            // which performs multiple forward Euler steps...?
            if(force_rho)
            {
                residual(idx) -= (*force_rho)(x(i), time(), p);
            }

            // c
            idx = c*n_dof_per_var_ + index(i);

            if(p.imposed_c_profile == true)
            {
                residual(idx) += u(c, index(i)) - c_value(x(i));
            }
            else
            {
                // Time derivative
                residual(idx) +=
                    1.0/dt_*(u(c, index(i)) - u(1, c, index(i)));

                // TODO double check this for 1d!

                // Laplacian in x-direction

                // Test using a flux-based Laplacian, like for the rho_bar equation.
                // Equivalent to the standard second-order stencil in the bulk,
                // but allows us to set the flux to zero at the borders
                double local_fc_p = 0.0;
                double local_fc_m = 0.0;

                if(i != 1)
                    local_fc_m = fc_m(i);
                if(i != p.n_interior_cell_1d)
                    local_fc_p = fc_p(i);
                
                residual(idx) -= -(local_fc_p - local_fc_m)/p.dx;

                // Standard second-order Laplacian stencil at cell centres
                //for(auto [k, w] : stencil::central_2::weights)
                //{
                    //residual(idx) -= w*u(1, c, index(i+k))/(p.dx*p.dx);
                //}

                // Using the cell-centred Laplacian stencil, which is the same
                // as regular 3-point stencil in the interior. At
                // boundaries/corners we use the modified stencil from the
                // notes by Long that automatically imposes Neumann BCs (double
                // check the definitions in Long)
                //if(i != 1 && i != p.n_interior_cell_1d)
                //{
                    //for(auto [k, w] : stencil::central_2::weights)
                    //{
                        //residual(idx) -= w*u(1, c, index(i+k))/(p.dx*p.dx);
                    //}
                //}
                //else if(i == 1)
                //{
                    //for(auto [k, w] : stencil::first_order::forward_1::weights)
                    //{
                        //residual(idx) -= w*u(1, c, index(i+k))/(p.dx*p.dx);
                    //}
                //}
                //else if(i == p.n_interior_cell_1d)
                //{
                    //for(auto [k, w] : stencil::first_order::backward_1::weights)
                    //{
                        //residual(idx) -= w*u(1, c, index(i+k))/(p.dx*p.dx);
                    //}
                //}

                // Reaction terms
                residual(idx) -=
                    - p.gamma_c*u(1, c, index(i))
                    + p.gamma_rho*u(1, rho_bar, index(i));

                // Forcing (for testing exact solutions)
                // TODO does the time() here work with the SSP RK timestepper,
                // which performs multiple forward Euler steps...?
                if(force_c)
                {
                    residual(idx) -= (*force_c)(x(i), time(), p);
                }
            }
        }
    }

    void calculate_jacobian(std::vector<Triplet> &triplet_list) const override
    {
        // Jacobian entries for the ghost cells are all ones down the diagonal
        {
            // Left boundary
            unsigned i = 0;

            unsigned idx = rho_bar*n_dof_per_var_ + index(i);
            triplet_list.emplace_back(idx, idx, 1.0);

            idx = c*n_dof_per_var_ + index(i);
            triplet_list.emplace_back(idx, idx, 1.0);

            // Right boundary
            i = p.n_interior_cell_1d+1;

            idx = rho_bar*n_dof_per_var_ + index(i);
            triplet_list.emplace_back(idx, idx, 1.0);

            idx = c*n_dof_per_var_ + index(i);
            triplet_list.emplace_back(idx, idx, 1.0);
        }

        // Loop over the interior cells
        for(unsigned i = 1; i <= p.n_interior_cell_1d; ++i)
        {
            // rho_bar
            unsigned idx = rho_bar*n_dof_per_var_ + index(i);
            triplet_list.emplace_back(idx, idx, 1.0/dt_);

            // c
            idx = c*n_dof_per_var_ + index(i);
            if(p.imposed_c_profile == true)
            {
                triplet_list.emplace_back(idx, idx, 1.0);
            }
            else
            {
                triplet_list.emplace_back(idx, idx, 1.0/dt_);
            }
        }
    }

    double c_value(double x) const
    {
        //const double s = x*p.dx;
        (void) x;
        return 1.0;
        //static const double pi = std::acos(-1.0);
        //return 1.0 - 0.5*(1.0 + std::cos(pi*s));
    }
};

int main_exact_solution_test(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    Config cf;
    cf.parse_command_line(argc, argv);

    log::set_level(cf.get_or<std::string>("log_level", "info"));

    KellerSegelParameters p;

    p.initial_conditions = KellerSegelICs::exact;
    p.D = 1.0;
    p.r = 0.0;
    p.chi = cf.get_or("chi", 1.0);
    p.gamma_rho = cf.get_or("gamma_rho", 1.0);
    p.gamma_c = cf.get_or("gamma_c", 1.0);
    p.L = exact_solution::L;
    p.n_interior_cell_1d = cf.get_or("n", 3u);

    const double t_max = cf.get_or("t_max", 1.0);
    const double dt = cf.get_or("dt", 0.1);
    const unsigned output_interval = cf.get_or("output_interval", 1u);

    KellerSegelProblem1D problem(p);

    exact_solution::omega = cf.get_or("omega", 1.0);

    problem.exact_rho = exact_solution::exact_rho;
    problem.exact_c = exact_solution::exact_c;

    problem.force_rho = exact_solution::force_rho;
    problem.force_c = exact_solution::force_c;

    cf.print_all();

    problem.set_initial_conditions();

    char filename[200];
    std::ofstream outfile;

    std::sprintf(filename, "output_%05u.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() < t_max)
    {
        problem.ssp_rk_3_timestep(dt);
        //problem.unsteady_solve(dt);

        if(i % output_interval == 0)
        {
            std::sprintf(filename, "output_%05u.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }

    return 0;
}

int main_normal(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    Config cf;
    cf.parse_command_line(argc, argv);

    log::set_level(cf.get_or<std::string>("log_level", "info"));

    KellerSegelParameters p;

    p.initial_conditions = KellerSegelICs::perturbed;
    p.imposed_c_profile = cf.get_or("imposed_c_profile", false);
    p.rho_bar_init = cf.get_or("rho_bar_init", 1.0);
    p.c_init = cf.get_or("c_init", 1.0);
    p.D = cf.get_or("D", 1.0);
    p.r = cf.get_or("r", 0.0);
    p.chi = cf.get_or("chi", 1.0);
    p.gamma_rho = cf.get_or("gamma_rho", 1.0);
    p.gamma_c = cf.get_or("gamma_c", 1.0);
    p.n_interior_cell_1d = cf.get_or("n", 3u);
    p.L = cf.get_or("L", 1.0);

    const double t_max = cf.get_or("t_max", 1.0);
    const double dt = cf.get_or("dt", 0.1);
    const unsigned output_interval = cf.get_or("output_interval", 1u);

    cf.print_all();

    KellerSegelProblem1D problem(p);

    problem.set_initial_conditions();

    char filename[200];
    std::ofstream outfile;

    std::sprintf(filename, "output_%05u.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() < t_max)
    {
        problem.ssp_rk_3_timestep(dt);
        //problem.unsteady_solve(dt);

        if(i % output_interval == 0)
        {
            std::sprintf(filename, "output_%05u.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }

    return 0;
}

int main(int argc, char **argv)
{
    return main_exact_solution_test(argc, argv);
    //return main_normal(argc, argv);
}
