#include <problem.h>

#include <iostream>
#include <fstream>

using namespace mjrfd;

struct ChemokinesParams
{
    double p_u;
    double alpha;
    double beta;
    double gamma_u;
    double gamma_b;
    double D_su;
    double D_ju;
    double nu;
    double lambda;
};

class ChemokinesProblem1D : public Problem
{
    typedef Eigen::Triplet<double> T;

public:
    // dofs:
    // C_u - 0--n_node-1
    // C_b - n_node--2*n_node-1
    // C_s - 2*n_node--3*n_node-1
    // phi - 3*n_node--4*n_node-1
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        Problem(4*n_node, dt),
        n_node_(n_node),
        dx_(1.0/(n_node-1))
    {
        std::cout << "n_node = " << n_node_ << '\n';
        std::cout << "n_dof  = " << n_dof_ << '\n';
        std::cout << "dx     = " << dx_ << '\n';
        std::cout << "dt     = " << dt_ << '\n';

        // set zero initial conditions
        // TODO implement actual ICs
        for(unsigned i = 0; i < 4*n_node_; ++i)
        {
            u(i) = 0.0;
        }

        u(0) = 1.0;
        u(4*n_node_-1) = 1.0;
    }

    ~ChemokinesProblem1D()
    {
    }

    ChemokinesParams p;

    double c_u(unsigned t, unsigned i)
    {
        return u(t, i);
    }

    double c_b(unsigned t, unsigned i)
    {
        return u(t, i + n_node_);
    }

    double c_s(unsigned t, unsigned i)
    {
        return u(t, i + 2*n_node_);
    }

    double phi(unsigned t, unsigned i)
    {
        return u(t, i + 3*n_node_);
    }

    void output(std::ostream &out) const
    {
        for(unsigned i = 0; i < n_node_; ++i)
        {
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' '
                << u(i) << ' '
                << u(n_node_   + i) << ' '
                << u(2*n_node_ + i) << ' '
                << u(3*n_node_ + i) << '\n';
        }
    }

private:
    unsigned n_node_;
    double dx_;

    void calculate_residual()
    {
        // set all entries to zero
        // not neccessary here since all entries are set explicitly
        residual_.setZero();

        // C_u
        residual_(0) += (c_u(0,0) - 1.0)*dx_*dx_;

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // time derivative
            residual_(i) += -(c_u(0,i) - c_u(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(i) += 0.5*(c_u(0,i-1) - 2.0*c_u(0,i) + c_u(0,i+1));
            residual_(i) += 0.5*(c_u(1,i-1) - 2.0*c_u(1,i) + c_u(1,i+1));

            // advection (Crank-Nicolson)
            residual_(i) += -0.5*p.p_u*(-0.5*c_u(0,i-1) + 0.5*c_u(0,i+1))*dx_;
            residual_(i) += -0.5*p.p_u*(-0.5*c_u(1,i-1) + 0.5*c_u(1,i+1))*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(i) += 0.5*(-p.alpha*c_u(0,i) + p.beta*c_b(0,i) - p.gamma_u*c_u(0,i)*phi(0,i))*dx_*dx_;
            residual_(i) += 0.5*(-p.alpha*c_u(1,i) + p.beta*c_b(1,i) - p.gamma_u*c_u(1,i)*phi(1,i))*dx_*dx_;
        }

        residual_(n_node_-1) += (c_u(0,n_node_-1) - 0.0)*dx_*dx_;

        // C_b (only ODEs so no boundary conditions)
        for(unsigned i = 0; i < n_node_; ++i)
        {
            // time derivative
            residual_(n_node_ + i) += -(c_b(0,i) - c_b(1,i))*dx_*dx_/dt_;

            // binding/unbinding (Crank-Nicolson)
            residual_(n_node_ + i) += 0.5*(p.alpha*c_u(0,i) - p.beta*c_b(0,i) - p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(n_node_ + i) += 0.5*(p.alpha*c_u(1,i) - p.beta*c_b(1,i) - p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;
        }

        // C_s
        residual_(2*n_node_) += (-1.5*c_s(0,0) + 2.0*c_s(0,1) - 0.5*c_s(0,2))*dx_ - p.p_u*c_s(0,0)*dx_*dx_;

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // time derivative
            residual_(2*n_node_ + i) += -(c_s(0,i) - c_s(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(2*n_node_ + i) += 0.5*p.D_su*(c_s(0,i-1) - 2.0*c_s(0,i) + c_s(0,i+1));
            residual_(2*n_node_ + i) += 0.5*p.D_su*(c_s(1,i-1) - 2.0*c_s(1,i) + c_s(1,i+1));

            // advection (Crank-Nicolson)
            residual_(2*n_node_ + i) += -0.5*p.p_u*(-0.5*c_s(0,i-1) + 0.5*c_s(0,i+1))*dx_;
            residual_(2*n_node_ + i) += -0.5*p.p_u*(-0.5*c_s(1,i-1) + 0.5*c_s(1,i+1))*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(2*n_node_ + i) += 0.5*(p.gamma_u*c_u(0,i)*phi(0,i) + p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(2*n_node_ + i) += 0.5*(p.gamma_u*c_u(1,i)*phi(1,i) + p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;
        }

        residual_(3*n_node_ - 1) += (0.5*c_s(0,n_node_-3) - 2.0*c_s(0,n_node_-2) + 1.5*c_s(0,n_node_-1))*dx_ - p.p_u*c_s(0,n_node_-1)*dx_*dx_;

        // phi
        residual_(3*n_node_) += (-1.5*phi(0,0) + 2.0*phi(0,1) - 0.5*phi(0,2) - p.lambda*phi(0,0))*dx_*dx_;

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // time derivative
            residual_(3*n_node_ + i) += -(phi(0,i) - phi(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(3*n_node_ + i) += 0.5*p.D_ju*(phi(0,i-1) - 2.0*phi(0,i) + phi(0,i+1));
            residual_(3*n_node_ + i) += 0.5*p.D_ju*(phi(1,i-1) - 2.0*phi(1,i) + phi(1,i+1));

            // advection (Chemotaxis?) (Crank-Nicolson)
            residual_(3*n_node_ + i) += -0.5*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1))*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1));
            residual_(3*n_node_ + i) += -0.5*p.nu*(-0.5*phi(1,i-1) + 0.5*phi(1,i+1))*(-0.5*c_b(1,i-1) + 0.5*c_b(1,i+1));

            residual_(3*n_node_ + i) += -0.5*p.nu*phi(0,i)*(c_b(0,i-1) - 2.0*c_b(0,i) + c_b(0,i+1));
            residual_(3*n_node_ + i) += -0.5*p.nu*phi(1,i)*(c_b(1,i-1) - 2.0*c_b(1,i) + c_b(1,i+1));
        }

        residual_(4*n_node_ - 1) += (phi(0,n_node_-1) - 1.0)*dx_*dx_;
    }

    void calculate_jacobian()
    {
        // set all entries to zero
        //jacobian_.setZero();

        std::vector<T> triplet_list;
        // TODO reserve correct amount
        triplet_list.reserve(5*n_node_);

        // --------------------------------------------------------------------

        // C_u / C_u
        triplet_list.push_back( T(0, 0, 1.0*dx_*dx_) );

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(i, i-1,  0.5 + 0.5*p.p_u*0.5*dx_) );
            triplet_list.push_back( T(i, i,   -dx_*dx_/dt_ - 0.5*2.0 - 0.5*(p.alpha + p.gamma_u*phi(0,i))*dx_*dx_) );
            triplet_list.push_back( T(i, i+1,  0.5 - 0.5*p.p_u*0.5*dx_) );
        }

        triplet_list.push_back( T(n_node_-1, n_node_-1, 1.0*dx_*dx_) );

        // C_u / C_b
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(i, n_node_ + i, 0.5*p.beta) );
        }

        // C_u / C_s
        // (no entries)
        
        // C_u / phi
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(i, 3*n_node_ + i, -0.5*p.gamma_u*c_u(0,i)*dx_*dx_) );
        }

        // --------------------------------------------------------------------

        // C_b / C_u
        for(unsigned i = 0; i < n_node_; ++i)
        {
            triplet_list.push_back( T(n_node_ + i, i, 0.5*p.alpha*dx_*dx_) );
        }

        // C_b / C_b
        for(unsigned i = 0; i < n_node_; ++i)
        {
            triplet_list.push_back( T(n_node_ + i, n_node_ + i, -dx_*dx_/dt_ - 0.5*(p.beta + p.gamma_b*phi(0,i))*dx_*dx_) );
        }

        // C_b / C_s
        // (no entries)

        // C_b / phi
        for(unsigned i = 0; i < n_node_; ++i)
        {
            triplet_list.push_back( T(n_node_ + i, 3*n_node_ + i, -0.5*(p.gamma_b*c_b(0,i))*dx_*dx_) );
        }

        // --------------------------------------------------------------------

        // C_s / C_u
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(2*n_node_ + i, i, 0.5*(p.gamma_u*phi(0,i))*dx_*dx_) );
        }

        // C_s / C_b
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(2*n_node_ + i, n_node_ + i, 0.5*(p.gamma_b*phi(0,i)*dx_*dx_)) );
        }

        // C_s / C_s
        triplet_list.push_back( T(2*n_node_, 2*n_node_,  -1.5*dx_ - p.p_u*dx_*dx_) );
        triplet_list.push_back( T(2*n_node_, 2*n_node_+1, 2.0*dx_) );
        triplet_list.push_back( T(2*n_node_, 2*n_node_+2, -0.5*dx_) );
        
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(2*n_node_ + i, 2*n_node_ + i-1, p.D_su*0.5 + 0.5*p.p_u*0.5*dx_) );
            triplet_list.push_back( T(2*n_node_ + i, 2*n_node_ + i,   -dx_*dx_/dt_ - 0.5*2.0*p.D_su) );
            triplet_list.push_back( T(2*n_node_ + i, 2*n_node_ + i+1, p.D_su*0.5 - 0.5*p.p_u*0.5*dx_) );
        }

        triplet_list.push_back( T(3*n_node_-1, 3*n_node_-3,  0.5*dx_) );
        triplet_list.push_back( T(3*n_node_-1, 3*n_node_-2, -2.0*dx_) );
        triplet_list.push_back( T(3*n_node_-1, 3*n_node_-1,  1.5*dx_ - p.p_u*dx_*dx_) );

        // C_s / phi
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(2*n_node_ + i, 3*n_node_ + i, 0.5*(p.gamma_u*c_u(0,i) + p.gamma_b*c_b(0,i))*dx_*dx_) );
        }

        // --------------------------------------------------------------------

        // phi / C_u
        // (no entries)

        // phi / C_b
        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(3*n_node_ + i, n_node_ + i-1,  0.5*0.5*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1)) - 0.5*p.nu*phi(0,i)) );
            triplet_list.push_back( T(3*n_node_ + i, n_node_ + i,    0.5*(2.0*p.nu*phi(0,i))) );
            triplet_list.push_back( T(3*n_node_ + i, n_node_ + i+1, -0.5*0.5*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1)) - 0.5*p.nu*phi(0,i)) );
        }

        // phi / C_s
        // (no entries)
        
        // phi / phi
        triplet_list.push_back( T(3*n_node_, 3*n_node_,   (-1.5 - p.lambda)*dx_*dx_) );
        triplet_list.push_back( T(3*n_node_, 3*n_node_+1, 2.0*dx_*dx_) );
        triplet_list.push_back( T(3*n_node_, 3*n_node_+2, -0.5*dx_*dx_) );

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            triplet_list.push_back( T(3*n_node_ + i, 3*n_node_ + i-1, 0.5*p.D_ju + 0.5*p.nu*0.5*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1))) );
            triplet_list.push_back( T(3*n_node_ + i, 3*n_node_ + i,   -dx_*dx_/dt_ - 2.0*0.5*p.D_ju - 0.5*p.nu*(c_b(0,i-1) - 2.0*c_b(0,i) + c_b(0,i+1))) );
            triplet_list.push_back( T(3*n_node_ + i, 3*n_node_ + i+1, 0.5*p.D_ju - 0.5*p.nu*0.5*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1))) );
        }

        triplet_list.push_back( T(4*n_node_-1, 4*n_node_-1, dx_*dx_) );

        // construct matrix
        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }
};

int main(int argc, char **argv)
{
    if(argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[1]);
    const double dt = std::atof(argv[2]);
    const double t_max = std::atof(argv[3]);

    unsigned output_interval = 1;

    if(argc == 5)
    {
        output_interval = std::atoi(argv[4]);
    }

    ChemokinesProblem1D problem(n_node, dt);
    ChemokinesProblem1D::Max_residual = 1e-8;

    // Holly Jackson's params
    //problem.p.p_u     = 1.0;
    //problem.p.alpha   = 2.0;
    //problem.p.beta    = 1.0;
    //problem.p.gamma_u = 1.0;
    //problem.p.gamma_b = 1.0;
    //problem.p.D_su    = 0.01;
    //problem.p.D_ju    = 1.0;
    //problem.p.nu      = 1.0;
    //problem.p.lambda  = 1.0;

    // checking
    problem.p.p_u     = 1.0;
    problem.p.alpha   = 2.0;
    problem.p.beta    = 1.0;
    problem.p.gamma_u = 1.0;
    problem.p.gamma_b = 1.0;
    problem.p.D_su    = 0.01;
    problem.p.D_ju    = 1.0;
    problem.p.nu      = 1.0;
    problem.p.lambda  = 1.0;

    char filename[200];
    std::ofstream outfile;

    // output initial conditions
    std::sprintf(filename, "output_%i.dat", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    while(problem.time() <= t_max)
    {
        // solve for current timestep
        bool dump = false;
        //bool dump = true;
        problem.unsteady_solve(dump);

        if(i % output_interval == 0)
        {
            std::cout << "Outputting solution at time = " << problem.time() << '\n';
            // output current solution
            std::sprintf(filename, "output_%i.dat", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }
}
