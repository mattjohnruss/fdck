#include "problem.h"
#include "stencil.h"
#include "log.h"
#include "config.h"

#include <tuple>
#include <fstream>

using namespace fdck;

enum Variable
{
    c = 0,
};

class DiffusionProblem2D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    DiffusionProblem2D(const unsigned n_node_1d) :
        Problem(1, n_node_1d*n_node_1d),
        dx_(1.0/(n_node_1d-1)),
        n_node_1d_(n_node_1d)
    {
        Max_residual = 1.0e-11;
    }

    ~DiffusionProblem2D()
    {
    }

    void output(std::ostream &out) const override
    {
        for(unsigned i = 0; i < n_node_1d_; ++i)
        {
            out << "t x y c\n";
            double x = this->x(i);

            for(unsigned j = 0; j < n_node_1d_; ++j)
            {
                double y = this->y(j);

                out << time() << ' '
                    << x << ' '
                    << y << ' '
                    << u(c, i + j*n_node_1d_) << '\n';
            }

            if(i != n_node_1d_-1)
            {
                out << '\n';
            }
        }
    }

    double x(unsigned i) const
    {
        return static_cast<double>(i)/static_cast<double>(n_node_1d_-1);
    }

    double y(unsigned j) const
    {
        return static_cast<double>(j)/static_cast<double>(n_node_1d_-1);
    }

    void set_initial_conditions()
    {
        clear_solution();

        for(unsigned i = 0; i < n_node_1d_; ++i)
        {
            for(unsigned j = 0; j < n_node_1d_; ++j)
            {
                //u(c, index_2d(0, j)) = 1.0;
                double x_c = this->x(i) - 0.5;
                double y_c = this->y(j) - 0.5;

                u(c, index_2d(i, j)) = 1000.0*std::exp(-100.0*(x_c*x_c + y_c*y_c));
            }
        }
    }

private:
    double dx_;
    unsigned n_node_1d_;

    unsigned index_2d(unsigned i, unsigned j) const
    {
        return i + j*n_node_1d_;
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        // Left and right BCs
        for(unsigned j = 0; j < n_node_1d_; ++j)
        {
            //residual(index_2d(0, j)) += (u(c, index_2d(0, j)) - 1.0)*dx_*dx_;
            residual(index_2d(0, j)) += (u(c, index_2d(0, j)))*dx_*dx_;
            residual(index_2d(n_node_1d_-1, j)) += u(c, index_2d(n_node_1d_-1, j))*dx_*dx_;
        }

        // Top and bottom BCs
        for(unsigned i = 0; i < n_node_1d_; ++i)
        {
            residual(index_2d(i, 0)) += u(c, index_2d(i, 0))*dx_*dx_;
            residual(index_2d(i, n_node_1d_-1)) += u(c, index_2d(i, n_node_1d_-1))*dx_*dx_;
        }

        // Bulk
        for(unsigned i = 1; i < n_node_1d_ - 1; ++i)
        {
            for(unsigned j = 1; j < n_node_1d_ - 1; ++j)
            {
                // Time derivative
                residual(index_2d(i, j)) += - (u(c, index_2d(i, j)) - u(1, c, index_2d(i, j)))*dx_*dx_/dt_;

                for(auto [k, w] : stencil::central_2::weights)
                {
                    // Laplacian in x-direction
                    residual(index_2d(i, j)) += w*u(c, index_2d(i+k, j));

                    // Laplacian in y-direction
                    residual(index_2d(i, j)) += w*u(c, index_2d(i, j+k));
                }
            }
        }
    }

    void calculate_jacobian(std::vector<Triplet> &triplet_list) const override
    {
        triplet_list.reserve(4*n_node_1d_ + (1 + 3*2)*(n_node_1d_ - 2)*(n_node_1d_ - 2));

        // Left and right BCs
        for(unsigned j = 0; j < n_node_1d_; ++j)
        {
            triplet_list.emplace_back(index_2d(0, j), index_2d(0, j), dx_*dx_);
            triplet_list.emplace_back(index_2d(n_node_1d_-1, j), index_2d(n_node_1d_-1, j), dx_*dx_);
        }

        // Top and bottom BCs
        for(unsigned i = 0; i < n_node_1d_; ++i)
        {
            triplet_list.emplace_back(index_2d(i, 0), index_2d(i, 0), dx_*dx_);
            triplet_list.emplace_back(index_2d(i, n_node_1d_-1), index_2d(i, n_node_1d_-1), dx_*dx_);
        }

        for(unsigned i = 1; i < n_node_1d_ - 1; ++i)
        {
            for(unsigned j = 1; j < n_node_1d_ - 1; ++j)
            {
                // Time derivative
                triplet_list.emplace_back(index_2d(i, j), index_2d(i, j), -dx_*dx_/dt_);

                for(auto [k, w] : stencil::central_2::weights)
                {
                    // Laplacian in x-direction
                    triplet_list.emplace_back(index_2d(i, j), index_2d(i+k, j), w);

                    // Laplacian in y-direction
                    triplet_list.emplace_back(index_2d(i, j), index_2d(i, j+k), w);
                }
            }
        }
    }
};

int main(int argc, char **argv)
{
    Config cf;
    cf.parse_command_line(argc, argv);

    unsigned n_node_1d = cf.get_or("n_node_1d", 11u);
    double dt = cf.get_or("dt", 0.1);
    double t_max = cf.get_or("t_max", 1.0);

    DiffusionProblem2D problem(n_node_1d);
    
    problem.set_initial_conditions();
    //problem.disable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    std::sprintf(filename, "output_%05u.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() < t_max)
    {
        problem.unsteady_solve(dt);

        std::sprintf(filename, "output_%05u.csv", i);
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        ++i;
    }

    return 0;
}
