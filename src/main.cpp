#include <problem.h>

#include <iostream>
#include <fstream>

using namespace mjrfd;

class DiffusionProblem1D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // treat all nodes as dofs in this example
    // we could eliminate the first and last nodes since they are Dirichlet
    // BCs, but this will require a sophisticated dof numbering scheme etc
    DiffusionProblem1D(unsigned n_node) :
        Problem(n_node), n_node_(n_node), dx_(1.0/(n_node-1))
    {
        std::cout << "n_node = " << n_node_ << '\n';
        std::cout << "dx = " << dx_ << '\n';
    }

    ~DiffusionProblem1D()
    {
    }

    void output(std::ostream &out) const
    {
        for(unsigned i = 0; i < n_node_; ++i)
        {
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' '
                << u_(i) << '\n';
        }
    }

    void output_exact(std::ostream &out) const
    {
        for(unsigned i = 0; i < n_node_; ++i)
        {
            double x = static_cast<double>(i)/static_cast<double>(n_node_-1);

            out << x << " " << 1.0 - x << '\n';
        }
    }

private:
    unsigned n_node_;
    double dx_;

    void calculate_residual()
    {
        // set all entries to zero
        // not neccessary here since all entries are set explicitly
        //residual_.setZero();

        residual_(0) = (u_(0) - 1.0)*dx_*dx_;

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            residual_(i) = (u_(i-1) - 2*u_(i) + u_(i+1));
        }

        residual_(n_node_-1) = u_(n_node_-1)*dx_*dx_;
    }

    void calculate_jacobian()
    {
        // set all entries to zero
        //jacobian_.setZero();

        std::vector<T> triplet_list;
        triplet_list.reserve(3*(n_node_-2) + 2);

        triplet_list.push_back( T(0, 0, 1.0*dx_*dx_) );

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(i, i-1,  1.0) );
            triplet_list.push_back( T(i, i,   -2.0) );
            triplet_list.push_back( T(i, i+1,  1.0) );
        }

        triplet_list.push_back( T(n_node_-1, n_node_-1, 1.0*dx_*dx_) );
        
        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }
};

int main()
{
    const unsigned n_node = 1001;

    DiffusionProblem1D problem(n_node);
    DiffusionProblem1D::Max_residual = 1e-14;

    problem.solve();

    std::ofstream outfile("output.dat");
    outfile.precision(std::numeric_limits<double>::max_digits10);

    problem.output(outfile);
    outfile.close();

    outfile.open("output_exact.dat");
    problem.output_exact(outfile);
    outfile.close();

    //// set up mesh with 101 nodes between 0 and 1
    //Mesh1D mesh(0.0, 1.0, n_node);

    //for(i = 1; i < n_node-1; ++i)
    //{
        //mesh.add_equation()
    //}

    //Problem problem;
    //problem.add_mesh(mesh);
    //problem.solve();
    
    //return test();
}
