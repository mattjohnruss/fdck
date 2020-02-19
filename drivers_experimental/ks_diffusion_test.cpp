#include "problem.h"
#include "stencil.h"
#include "log.h"
#include "config.h"
#include "utilities.h"

#include <tuple>
#include <fstream>

#include <fenv.h>

using namespace fdck;

enum Variable
{
    rho_bar = 0,
};

// Hybrid FVFD method for Keller-Segel based on Chertock et al (2017)
class KellerSegelProblem2D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // Square grid of cells, with one layer of ghost cells added outside each
    // boundary.
    // We need no aux dofs but we need 3 history values for SSP R-K timestepping.
    // TODO check if dx, dy are correctly defined
    KellerSegelProblem2D(const unsigned n_interior_cell_1d) :
        Problem(1, n_interior_cell_1d*n_interior_cell_1d + 4*n_interior_cell_1d, 0, 3),
        dx_(1.0/(n_interior_cell_1d)),
        dy_(1.0/(n_interior_cell_1d)),
        n_interior_cell_1d_(n_interior_cell_1d)
    {
        FDCK_INFO("Creating grid with {} interior cells in each direction ({} including ghost cells).", n_interior_cell_1d, n_interior_cell_1d_+2);
        Max_residual = 1.0e-14;
        // TODO investigate why setting this to 1 causes (possibly) no
        // iterations... at the very least the logging is confusing
        Max_newton_iterations = 2;
    }

    ~KellerSegelProblem2D() override
    {
    }

    void output(std::ostream &out) const override
    {
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            out << "t x y rho\\\\_bar\n";
            double x = this->x(i);

            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                double y = this->y(j);

                out << time() << ' '
                    << x << ' '
                    << y << ' '
                    << u(rho_bar, index_2d(i, j)) << '\n';
            }

            if(i != n_interior_cell_1d_)
            {
                out << '\n';
            }
        }
    }

    void set_initial_conditions()
    {
        FDCK_TRACE("set_initial_conditions()");

        clear_solution();

        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                double x_c = this->x(i) - 0.5;
                double y_c = this->y(j) - 0.5;

                u(rho_bar, index_2d(i, j)) = 10.0*std::exp(-100.0*(x_c*x_c + y_c*y_c));
            }
        }
    }

    double x(unsigned i) const
    {
        return (static_cast<double>(i) - 0.5)*dx_;
    }

    double y(unsigned j) const
    {
        return (static_cast<double>(j) - 0.5)*dy_;
    }

private:
    enum class Face
    {
        North,
        East,
        South,
        West,
    };

    double dx_;
    double dy_;
    unsigned n_interior_cell_1d_;

    unsigned index_2d(unsigned i, unsigned j) const
    {
        FDCK_TRACE("index_2d(i = {}, j = {})", i, j);

        if(i >= 1 && i <= n_interior_cell_1d_ && j >= 1 && j <= n_interior_cell_1d_)
        {
            return (i-1) + (j-1)*n_interior_cell_1d_;
        }
        else if(i == 0 && j != 0 && j != n_interior_cell_1d_+1)
        {
            return index_ghost(j, Face::West);
        }
        else if(i == n_interior_cell_1d_+1 && j != 0 && j != n_interior_cell_1d_+1)
        {
            return index_ghost(j, Face::East);
        }
        else if(j == 0 && i != 0 && i != n_interior_cell_1d_+1)
        {
            return index_ghost(i, Face::South);
        }
        else if(j == n_interior_cell_1d_+1 && i != 0 && i != n_interior_cell_1d_+1)
        {
            return index_ghost(i, Face::North);
        }
        else
        {
            FDCK_FATAL("Cell index (i = {}, j = {}) invalid; exiting", i, j);
            std::abort();
        }
    }
    
    // Returns the node number for the i-th ghost cell on face f of the grid
    unsigned index_ghost(unsigned i, Face f) const
    {
        FDCK_TRACE("index_ghost(i = {}, f = {})", i, static_cast<unsigned>(f));
        const unsigned n_interior_cell = n_interior_cell_1d_*n_interior_cell_1d_;

        switch(f)
        {
            case Face::North:
                return n_interior_cell + (i-1);
            case Face::East:
                return n_interior_cell + n_interior_cell_1d_ + (i-1);
            case Face::South:
                return n_interior_cell + 2*n_interior_cell_1d_ + (i-1);
            case Face::West:
                return n_interior_cell + 3*n_interior_cell_1d_ + (i-1);
        }
    }

    // Equation (2.3)(a)
    // Calculates F_{i+1/2, j}
    double f_p(int i, int j) const
    {
        FDCK_TRACE("f_p(i = {}, j = {})", i, j);
        double drho_dx_p = drho_dx_p_at_midpoint(i, j);
        return -drho_dx_p;
    }

    // Calculates F_{i-1/2, j}
    double f_m(unsigned i, unsigned j) const
    {
        FDCK_TRACE("f_m(i = {}, j = {})", i, j);
        double drho_dx_m = drho_dx_p_at_midpoint(i-1, j);
        return -drho_dx_m;
    }

    // Equation (2.3)(b)
    // Calculates G_{i, j+1/2}
    double g_p(int i, int j) const
    {
        FDCK_TRACE("g_p(i = {}, j = {})", i, j);
        double drho_dy_p = drho_dy_p_at_midpoint(i, j);
        return -drho_dy_p;
    }

    // Calculates G_{i, j-1/2}
    double g_m(unsigned i, unsigned j) const
    {
        FDCK_TRACE("g_m(i = {}, j = {})", i, j);
        double drho_dy_m = drho_dy_p_at_midpoint(i, j-1);
        return -drho_dy_m;
    }

    // Equation (2.4)(a)
    double drho_dx_p_at_midpoint(unsigned i, unsigned j) const
    {
        FDCK_TRACE("drho_dy_p_at_midpoint(i = {}, j = {})", i, j);
        double drho_dx_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            drho_dx_p += w*u(1, rho_bar, index_2d(i+k, j))/dx_;
        }
        return drho_dx_p;
    }

    // Equation (2.4)(b)
    double drho_dy_p_at_midpoint(unsigned i, unsigned j) const
    {
        FDCK_TRACE("drho_dy_p_at_midpoint(i = {}, j = {})", i, j);
        double drho_dy_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            drho_dy_p += w*u(1, rho_bar, index_2d(i, j+k))/dy_;
        }
        return drho_dy_p;
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        for(unsigned b = 1; b <= n_interior_cell_1d_; ++b)
        {
            // Set the values in the ghost cells using quadratic interpolation
            // of the previous values in the neighbouring interior cells

            // Top boundary
            unsigned j = n_interior_cell_1d_+1;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}, j = {}) [top]", b, j);

            unsigned index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            residual(index) +=
                //u(rho_bar, index_2d(b, j)) - (3.0*u(1, rho_bar, index_2d(b, j-1)) - 3.0*u(1, rho_bar, index_2d(b, j-2)) + 1.0*u(1, rho_bar, index_2d(b, j-3)));
                u(rho_bar, index_2d(b, j));

            // Right boundary
            unsigned i = n_interior_cell_1d_+1;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}, j = {}) [right]", i, b);

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            residual(index) +=
                u(rho_bar, index_2d(i, b));

            // Bottom boundary
            j = 0;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}, j = {}) [bottom]", b, j);

            index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            residual(index) +=
                u(rho_bar, index_2d(b, j));

            // Left boundary
            i = 0;

            FDCK_TRACE("");
            FDCK_TRACE("calculate_residual(i = {}, j = {}) [left]", i, b);

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            residual(index) +=
                u(rho_bar, index_2d(i, b));
        }

        // Loop over the interior cells
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                FDCK_TRACE("");
                FDCK_TRACE("calculate_residual(i = {}, j = {}) [interior]", i, j);

                // rho_bar
                unsigned index = rho_bar*n_dof_per_var_ + index_2d(i, j);

                // Time derivative
                residual(index) +=
                    1.0/dt_*(u(rho_bar, index_2d(i, j)) - u(1, rho_bar, index_2d(i, j)));

                // Storage for the fluxes
                double local_f_m = 0.0;
                double local_f_p = 0.0;
                double local_g_m = 0.0;
                double local_g_p = 0.0;

                // Unless we're in a boundary cell, calculate the flux. If
                // we're in a boundary cell, the flux at the appropriate
                // boundary is left as 0.0 to implement the no-flux BCs
                if(i != 1)
                    local_f_m = f_m(i, j);
                if(i != n_interior_cell_1d_)
                    local_f_p = f_p(i, j);
                if(j != 1)
                    local_g_m = g_m(i, j);
                if(j != n_interior_cell_1d_)
                    local_g_p = g_p(i, j);

                residual(index) += (local_f_p - local_f_m)/dx_;
                residual(index) += (local_g_p - local_g_m)/dy_;
            }
        }
    }

    void calculate_jacobian(std::vector<Triplet> &triplet_list) const override
    {
        // Jacobian entries for the ghost cells are all ones down the diagonal
        for(unsigned b = 1; b <= n_interior_cell_1d_; ++b)
        {
            // Top boundary
            unsigned j = n_interior_cell_1d_+1;

            // rho_bar
            unsigned index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            triplet_list.emplace_back(index, index, 1.0);

            // Right boundary
            unsigned i = n_interior_cell_1d_+1;

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            triplet_list.emplace_back(index, index, 1.0);

            // Bottom boundary
            j = 0;

            index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            triplet_list.emplace_back(index, index, 1.0);

            // Left boundary
            i = 0;

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            triplet_list.emplace_back(index, index, 1.0);
        }

        // Loop over the interior cells
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                // rho_bar
                unsigned index = rho_bar*n_dof_per_var_ + index_2d(i, j);
                triplet_list.emplace_back(index, index, 1.0/dt_);
            }
        }
    }
};

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    Config cf;
    cf.parse_command_line(argc, argv);

    log::set_level("debug");

    unsigned n_interior_cell_1d = cf.get_or("n", 3u);
    double dt = cf.get_or("dt", 0.1);
    double t_max = cf.get_or("t_max", 1.0);

    KellerSegelProblem2D problem(n_interior_cell_1d);

    problem.set_initial_conditions();
    problem.disable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    std::sprintf(filename, "output_%05u.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() < t_max)
    {
        //problem.ssp_rk_3_timestep(dt);
        problem.unsteady_solve(dt);

        std::sprintf(filename, "output_%05u.csv", i);
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        ++i;
    }

    return 0;
}
