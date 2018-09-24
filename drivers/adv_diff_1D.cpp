#include <problem.h>

#include <iostream>
#include <fstream>

using namespace mjrfd;

enum Variable
{
    c = 0,
};

class DiffusionProblem1D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // treat all nodes as dofs in this example
    // we could eliminate the first and last nodes since they are Dirichlet
    // BCs, but this will require a sophisticated dof numbering scheme etc
    DiffusionProblem1D(const unsigned n_node) :
        Problem(1, n_node),
        n_node_(n_node),
        dx_(1.0/(n_node-1))
    {
        std::cout << "n_node = " << n_node_ << '\n';
        std::cout << "dx = " << dx_ << '\n';

        // set zero initial conditions
        for(unsigned i = 0; i < n_node_; ++i)
        {
            u(c,i) = 0.0;
        }

        Max_residual = 1e-14;
    }

    ~DiffusionProblem1D()
    {
    }

    void output(std::ostream &out) const override
    {
        out << "t x c\n";
        for(unsigned i = 0; i < n_node_; ++i)
        {
            out << time() << ' '
                << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' '
                << u(c,i) << '\n';
        }
    }

    void output_exact(std::ostream &out) const override
    {
        out << "t x c\n";
        for(unsigned i = 0; i < n_node_; ++i)
        {
            double x = static_cast<double>(i)/static_cast<double>(n_node_-1);

            out << time() << " " << x << " " << 1.0 - x << '\n';
        }
    }

private:
    unsigned n_node_;
    double dx_;

    void calculate_residual() override
    {
        // set all entries to zero
        // not neccessary here since all entries are set explicitly
        //residual_.setZero();

        residual_(0) = u(c,0)*dx_*dx_;

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // time derivative
            residual_(i) = -(u(c,i) - u(1,c,i))*dx_*dx_/dt_;

            // diffusion (u_ and u_old_ because of Crank-Nicolson)
            residual_(i) += 0.5*(u(c,i-1) - 2*u(c,i) + u(c,i+1));
            residual_(i) += 0.5*(u(1,c,i-1) - 2*u(1,c,i) + u(1,c,i+1));

            // advection (u_ and u_old_ because of Crank-Nicolson)
            residual_(i) -= 0.5*(6.0*(-0.5*u(c,i-1) + 0.5*u(c,i+1)))*dx_;
            residual_(i) -= 0.5*(6.0*(-0.5*u(1,c,i-1) + 0.5*u(1,c,i+1)))*dx_;

            // forcing
            residual_(i) += 8.0*dx_*dx_;
        }

        residual_(n_node_-1) = u(c,n_node_-1)*dx_*dx_;
    }

    void calculate_jacobian() override
    {
        // set all entries to zero
        //jacobian_.setZero();

        std::vector<T> triplet_list;
        triplet_list.reserve(3*(n_node_-2) + 2);

        triplet_list.push_back( T(0, 0, 1.0*dx_*dx_) );

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(i, i-1,  0.5 + 0.5*6.0*0.5*dx_) );
            triplet_list.push_back( T(i, i,   -1.0 - dx_*dx_/dt_) );
            triplet_list.push_back( T(i, i+1,  0.5 - 0.5*6.0*0.5*dx_) );
        }

        triplet_list.push_back( T(n_node_-1, n_node_-1, 1.0*dx_*dx_) );
        
        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }
};

int main()
{
    const unsigned n_node = 1001;
    const double dt = 0.001;

    DiffusionProblem1D problem(n_node);

    char filename[200];
    std::ofstream outfile;

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() < 1.0)
    {
        // solve for current timestep
        problem.unsteady_solve(dt);

        if(i % 10 == 0)
        {
            // output current solution
            std::sprintf(filename, "output_%05i.csv", i/10);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }
}
