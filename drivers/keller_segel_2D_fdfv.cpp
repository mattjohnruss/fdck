#include "problem.h"
#include "stencil.h"
#include "log.h"
#include "config.h"
#include "utilities.h"

#include <tuple>
#include <fstream>

using namespace mjrfd;

enum Variable
{
    rho_bar = 0,
    c = 1,
};

struct KellerSegelParameters
{
    double chi = 1.0;
    double gamma_rho = 1.0;
    double gamma_c = 1.0;
};

// Hybrid FVFD method for Keller-Segel based on Chertock et al (2017)
class KellerSegelProblem2D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // Square grid of cells, with one layer of ghost cells added outside each
    // boundary.
    // We need no aux dofs but we need 3 history values for SSP R-K timestepping.
    KellerSegelProblem2D(const unsigned n_interior_cell_1d) :
        Problem(2, n_interior_cell_1d*n_interior_cell_1d + 4*n_interior_cell_1d, 0, 3),
        dx_(1.0/(n_interior_cell_1d-1)),
        dy_(1.0/(n_interior_cell_1d-1)),
        n_interior_cell_1d_(n_interior_cell_1d)
    {
        MJRFD_INFO("Creating grid with {} interior cells in each direction ({} including ghost cells).", n_interior_cell_1d, n_interior_cell_1d_+2);
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
            out << "t x y rho_bar c\n";
            double x = this->x(i);

            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                double y = this->y(j);

                out << time() << ' '
                    << x << ' '
                    << y << ' '
                    << u(rho_bar, index_2d(i, j)) << ' '
                    << u(c, index_2d(i, j)) << '\n';
            }

            if(i != n_interior_cell_1d_)
            {
                out << '\n';
            }
        }
    }

    void set_initial_conditions()
    {
        MJRFD_TRACE("set_initial_conditions()");

        clear_solution();

        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                double x_c = this->x(i) - 0.5;
                double y_c = this->y(j) - 0.5;

                u(rho_bar, index_2d(i, j)) = 1000.0*std::exp(-100.0*(x_c*x_c + y_c*y_c));
                u(c, index_2d(i, j)) = 500.0*std::exp(-50.0*(x_c*x_c + y_c*y_c));
            }
        }
    }

    double x(unsigned i) const
    {
        return static_cast<double>(i-1)/static_cast<double>(n_interior_cell_1d_-1);
    }

    double y(unsigned j) const
    {
        return static_cast<double>(j-1)/static_cast<double>(n_interior_cell_1d_-1);
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
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                u(rho_bar, index_2d(i, j)) =
                    0.75*u(2, rho_bar, index_2d(i, j)) + 0.25*u(rho_bar, index_2d(i, j));
            }
        }

        // reset time before next solve
        time() = current_time;

        // take another forward Euler step to give the [ ... ] brackets in w^(3)
        unsteady_solve(dt);

        // calculate w^(3)
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                u(rho_bar, index_2d(i, j)) =
                    (1.0/3.0)*u(3, rho_bar, index_2d(i, j)) + (2.0/3.0)*u(rho_bar, index_2d(i, j));
            }
        }

        // update the time manually
        time() = current_time + dt;
    }

    KellerSegelParameters p;

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

    // TODO double check all dx_/dx_/dt denominators and multiply through as in
    // other problems

    // Returns the node number corresponding to the 2d grid location (i, j)
    // Indexing scheme is as follows.
    // Interior cells: i, j = 1, ..., n_interior_cell_1d_,
    // Ghost cells:
    //     North: j = n_interior_cell_1d_ + 1, i = 1, ..., n_interior_cell_1d_,
    //     East:  i = n_interior_cell_1d_ + 1, j = 1, ..., n_interior_cell_1d_,
    //     South: j = 0,                       i = 1, ..., n_interior_cell_1d_,
    //     North: j = 0,                       j = 1, ..., n_interior_cell_1d_.
    // Node ordering: interior cells are first (row-major), then ghost cells in
    // the order above.
    unsigned index_2d(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("index_2d(i = {}, j = {})", i, j);

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
            MJRFD_FATAL("Cell index (i = {}, j = {}) invalid; exiting", i, j);
            std::abort();
        }
    }
    
    // Returns the node number for the i-th ghost cell on face f of the grid
    unsigned index_ghost(unsigned i, Face f) const
    {
        MJRFD_TRACE("index_ghost(i = {}, f = {})", i, static_cast<unsigned>(f));
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
    double f_p(int i, int j) const
    {
        MJRFD_TRACE("f_p(i = {}, j = {})", i, j);
        double u_p = u_p_at_midpoint(i, j);
        double drho_dx_p = drho_dx_p_at_midpoint(i, j);
        return p.chi*rho_point_value_i_p(i, j)*u_p - drho_dx_p;
    }

    double f_m(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("f_m(i = {}, j = {})", i, j);
        // for the next two, we can call the plus (_p) functions with args i-1,
        // j-1 as it gives the same result as a manual impl of the minus functions
        double u_m = u_p_at_midpoint(i-1, j-1);
        double drho_dx_m = drho_dx_p_at_midpoint(i-1, j-1);
        // rho_point_value_i_m must be written explicitly because the upwinding
        // should be in the opposite direction than the plus case
        return p.chi*rho_point_value_i_m(i, j)*u_m - drho_dx_m;
    }

    // Equation (2.3)(b)
    double g_p(int i, int j) const
    {
        MJRFD_TRACE("g_p(i = {}, j = {})", i, j);
        double v_p = v_p_at_midpoint(i, j);
        double drho_dy_p = drho_dy_p_at_midpoint(i, j);
        return p.chi*rho_point_value_j_p(i, j)*v_p - drho_dy_p;
    }

    double g_m(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("g_m(i = {}, j = {})", i, j);
        // the comments in f_m apply equally to g_m
        double v_m = v_p_at_midpoint(i-1, j-1);
        double drho_dy_m = drho_dy_p_at_midpoint(i-1, j-1);
        return p.chi*rho_point_value_j_m(i, j)*v_m - drho_dy_m;
    }

    // Equation (2.4)(a)
    double drho_dx_p_at_midpoint(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("drho_dy_p_at_midpoint(i = {}, j = {})", i, j);
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
        MJRFD_TRACE("drho_dy_p_at_midpoint(i = {}, j = {})", i, j);
        double drho_dy_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            drho_dy_p += w*u(1, rho_bar, index_2d(i, j+k))/dy_;
        }
        return drho_dy_p;
    }

    // Equation (2.4)(c)
    double u_p_at_midpoint(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("u_p_at_midpoint(i = {}, j = {})", i, j);
        double u_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            u_p += w*u(1, c, index_2d(i+k, j))/dx_;
        }
        return u_p;
    }

    // Equation (2.4)(d)
    double v_p_at_midpoint(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("v_p_at_midpoint(i = {}, j = {})", i, j);
        double v_p = 0.0;
        for(auto [k, w] : stencil::first_order::forward_1::weights)
        {
            v_p += w*u(1, c, index_2d(i, j+k))/dy_;
        }
        return v_p;
    }

    // Equation (2.5)(a)
    double rho_point_value_i_p(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("rho_point_value_i_p(i = {}, j = {})", i, j);
        double u_p = u_p_at_midpoint(i, j);

        if(u_p > 0.0)
        {
            return rho_point_value_p_at_face(i, j, Face::East);
        }
        else
        {
            return rho_point_value_p_at_face(i+1, j, Face::West);
        }
    }

    double rho_point_value_i_m(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("rho_point_value_i_m(i = {}, j = {})", i, j);
        double u_m = u_p_at_midpoint(i-1, j-1);

        if(u_m < 0.0)
        {
            return rho_point_value_p_at_face(i, j, Face::East);
        }
        else
        {
            return rho_point_value_p_at_face(i-1, j, Face::West);
        }
    }

    // Equation (2.5)(b)
    double rho_point_value_j_p(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("rho_point_value_j_p(i = {}, j = {})", i, j);
        double v_p = v_p_at_midpoint(i, j);

        if(v_p > 0.0)
        {
            return rho_point_value_p_at_face(i, j, Face::North);
        }
        else
        {
            return rho_point_value_p_at_face(i, j+1, Face::South);
        }
    }

    double rho_point_value_j_m(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("rho_point_value_j_m(i = {}, j = {})", i, j);
        double v_m = v_p_at_midpoint(i-1, j-1);

        if(v_m < 0.0)
        {
            return rho_point_value_p_at_face(i, j, Face::North);
        }
        else
        {
            return rho_point_value_p_at_face(i, j-1, Face::South);
        }
    }

    // Equation (2.7)
    // We have shifted some of the indices here but it should be equivalent to
    // the paper 
    double rho_point_value_p_at_face(unsigned i, unsigned j, Face f) const
    {
        MJRFD_TRACE("rho_point_value_p_at_face(i = {}, j = {})", i, j);
        double result = u(1, rho_bar, index_2d(i, j));

        switch(f)
        {
            case Face::North:
                result += 0.5*dy_*drho_dy(i, j);
                break;
            case Face::South:
                result -= 0.5*dy_*drho_dy(i, j);
                break;
            case Face::East:
                result += 0.5*dx_*drho_dx(i, j);
                break;
            case Face::West:
                result -= 0.5*dx_*drho_dx(i, j);
                break;
        }

        return result;
    }

    // Equation (2.8)(a)
    // TODO What exactly does the +- mean in the conditions here? At the moment
    // we are assuming it means that both the + and - quantities must be >= 0
    // separately
    double drho_dx(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("drho_dx(i = {}, j = {})", i, j);
        double drho_bar_central = 0.0;
        for(auto [k, w] : stencil::central_1::weights)
        {
            drho_bar_central += w*u(1, rho_bar, index_2d(i+k, j));
        }

        const double test_p = u(1, rho_bar, index_2d(i, j)) + 0.5*drho_bar_central;
        const double test_m = u(1, rho_bar, index_2d(i, j)) - 0.5*drho_bar_central;

        if(test_p >= 0.0 && test_m >= 0.0)
        {
            return drho_bar_central/dx_;
        }
        else
        {
            double drho_bar_forward = 0.0;
            for(auto [k, w] : stencil::first_order::forward_1::weights)
            {
                drho_bar_forward += w*u(1, rho_bar, index_2d(i+k, j));
            }
            double drho_bar_backward = 0.0;
            for(auto [k, w] : stencil::first_order::backward_1::weights)
            {
                drho_bar_backward += w*u(1, rho_bar, index_2d(i+k, j));
            }
            return utilities::minmod<double>({ 2.0*drho_bar_forward/dx_,
                                               drho_bar_central/dx_,
                                               2.0*drho_bar_backward/dx_ });
        }
    }

    // Equation (2.8)(b)
    // TODO Same comment here as for (2.8)(a) regarding +-
    double drho_dy(unsigned i, unsigned j) const
    {
        MJRFD_TRACE("drho_dy(i = {}, j = {})", i, j);
        double drho_bar_central = 0.0;
        for(auto [k, w] : stencil::central_1::weights)
        {
            drho_bar_central += w*u(1, rho_bar, index_2d(i, j+k));
        }

        const double test_p = u(1, rho_bar, index_2d(i, j)) + 0.5*drho_bar_central;
        const double test_m = u(1, rho_bar, index_2d(i, j)) - 0.5*drho_bar_central;

        if(test_p >= 0.0 && test_m >= 0.0)
        {
            return drho_bar_central/dy_;
        }
        else
        {
            double drho_bar_forward = 0.0;
            for(auto [k, w] : stencil::first_order::forward_1::weights)
            {
                drho_bar_forward += w*u(1, rho_bar, index_2d(i, j+k));
            }
            double drho_bar_backward = 0.0;
            for(auto [k, w] : stencil::first_order::backward_1::weights)
            {
                drho_bar_backward += w*u(1, rho_bar, index_2d(i, j+k));
            }
            return utilities::minmod<double>({ 2.0*drho_bar_forward/dy_,
                                               drho_bar_central/dy_,
                                               2.0*drho_bar_backward/dy_ });
        }
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        for(unsigned b = 1; b <= n_interior_cell_1d_; ++b)
        {
            // Set the values in the ghost cells using quadratic interpolation
            // of the previous values in the neighbouring interior cells

            // Top boundary
            unsigned j = n_interior_cell_1d_+1;

            MJRFD_TRACE("");
            MJRFD_TRACE("calculate_residual(i = {}, j = {}) [top]", b, j);

            // rho_bar
            unsigned index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            residual(index) +=
                u(rho_bar, index_2d(b, j)) - (3.0*u(1, rho_bar, index_2d(b, j-1)) - 3.0*u(1, rho_bar, index_2d(b, j-2)) + 1.0*u(1, rho_bar, index_2d(b, j-3)));

            // c
            index = c*n_dof_per_var_ + index_2d(b, j);
            residual(index) += u(c, index_2d(b, j)) - u(1, c, index_2d(b, j-1));

            // Right boundary
            unsigned i = n_interior_cell_1d_+1;

            MJRFD_TRACE("");
            MJRFD_TRACE("calculate_residual(i = {}, j = {}) [right]", i, b);

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            residual(index) +=
                u(rho_bar, index_2d(i, b)) - (3.0*u(1, rho_bar, index_2d(i-1, b)) - 3.0*u(1, rho_bar, index_2d(i-2, b)) + 1.0*u(1, rho_bar, index_2d(i-3, b)));

            index = c*n_dof_per_var_ + index_2d(i, b);
            residual(index) += u(c, index_2d(i, b)) - u(1, c, index_2d(i-1, b));

            // Bottom boundary
            j = 0;

            MJRFD_TRACE("");
            MJRFD_TRACE("calculate_residual(i = {}, j = {}) [bottom]", b, j);

            index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            residual(index) +=
                u(rho_bar, index_2d(b, j)) - (3.0*u(1, rho_bar, index_2d(b, j+1)) - 3.0*u(1, rho_bar, index_2d(b, j+2)) + 1.0*u(1, rho_bar, index_2d(b, j+3)));

            index = c*n_dof_per_var_ + index_2d(b, j);
            residual(index) += u(c, index_2d(b, j)) - u(1, c, index_2d(b, j+1));

            // Left boundary
            i = 0;

            MJRFD_TRACE("");
            MJRFD_TRACE("calculate_residual(i = {}, j = {}) [left]", i, b);

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            residual(index) +=
                u(rho_bar, index_2d(i, b)) - (3.0*u(1, rho_bar, index_2d(i+1, b)) - 3.0*u(1, rho_bar, index_2d(i+2, b)) + 1.0*u(1, rho_bar, index_2d(i+3, b)));

            index = c*n_dof_per_var_ + index_2d(i, b);
            residual(index) += u(c, index_2d(i, b)) - u(1, c, index_2d(i+1, b));
        }

        // Loop over the interior cells
        for(unsigned i = 1; i <= n_interior_cell_1d_; ++i)
        {
            for(unsigned j = 1; j <= n_interior_cell_1d_; ++j)
            {
                MJRFD_TRACE("");
                MJRFD_TRACE("calculate_residual(i = {}, j = {}) [interior]", i, j);

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

                residual(index) += - (local_f_p - local_f_m)/dx_;
                residual(index) += - (local_g_p - local_g_m)/dy_;

                // c
                index = c*n_dof_per_var_ + index_2d(i, j);

                // Time derivative
                residual(index) +=
                    1.0/dt_*(u(c, index_2d(i, j)) - u(1, c, index_2d(i, j)));

                // Using the cell-centred Laplacian stencil, which is the same
                // as regular 5-point stencil in the interior. At
                // boundaries/corners we use the modified stencil from the
                // notes by Long that automatically imposes Neumann BCs (double
                // check the definitions in Long)

                // Laplacian in x-direction
                if(i != 1 && i != n_interior_cell_1d_)
                {
                    for(auto [k, w] : stencil::central_2::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i+k, j))/(dx_*dx_);
                    }
                }
                else if(i == 1)
                {
                    for(auto [k, w] : stencil::first_order::forward_1::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i+k, j))/(dx_*dx_);
                    }
                }
                else if(i == n_interior_cell_1d_)
                {
                    for(auto [k, w] : stencil::first_order::backward_1::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i+k, j))/(dx_*dx_);
                    }
                }

                // Laplacian in y-direction
                if(j != 1 && j != n_interior_cell_1d_)
                {
                    for(auto [k, w] : stencil::central_2::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i, j+k))/(dy_*dy_);
                    }
                }
                else if(j == 1)
                {
                    for(auto [k, w] : stencil::first_order::forward_1::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i, j+k))/(dy_*dy_);
                    }
                }
                else if(j == n_interior_cell_1d_)
                {
                    for(auto [k, w] : stencil::first_order::backward_1::weights)
                    {
                        residual(index) += w*u(1, c, index_2d(i, j+k))/(dy_*dy_);
                    }
                }

                // Reaction terms
                residual(index) +=
                    - p.gamma_c*u(1, c, index_2d(i, j))
                    + p.gamma_rho*u(1, rho_bar, index_2d(i, j));
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

            // c
            index = c*n_dof_per_var_ + index_2d(b, j);
            triplet_list.emplace_back(index, index, 1.0);

            // Right boundary
            unsigned i = n_interior_cell_1d_+1;

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            triplet_list.emplace_back(index, index, 1.0);

            index = c*n_dof_per_var_ + index_2d(i, b);
            triplet_list.emplace_back(index, index, 1.0);

            // Bottom boundary
            j = 0;

            index = rho_bar*n_dof_per_var_ + index_2d(b, j);
            triplet_list.emplace_back(index, index, 1.0);

            index = c*n_dof_per_var_ + index_2d(b, j);
            triplet_list.emplace_back(index, index, 1.0);

            // Left boundary
            i = 0;

            index = rho_bar*n_dof_per_var_ + index_2d(i, b);
            triplet_list.emplace_back(index, index, 1.0);

            index = c*n_dof_per_var_ + index_2d(i, b);
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

                // c
                index = c*n_dof_per_var_ + index_2d(i, j);
                triplet_list.emplace_back(index, index, 1.0/dt_);
            }
        }
    }
};

int main(int argc, char **argv)
{
    Config cf;
    cf.parse_command_line(argc, argv);

    log::set_level("info");

    unsigned n_interior_cell_1d = cf.get_or("n", 3u);
    double dt = cf.get_or("dt", 0.1);
    double t_max = cf.get_or("t_max", 1.0);

    KellerSegelProblem2D problem(n_interior_cell_1d);

    problem.p.chi = cf.get_or("chi", 1.0);
    problem.p.gamma_rho = cf.get_or("gamma_rho", 1.0);
    problem.p.gamma_c = cf.get_or("gamma_c", 1.0);
    
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
        problem.ssp_rk_3_timestep(dt);

        std::sprintf(filename, "output_%05u.csv", i);
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        ++i;
    }

    return 0;
}
