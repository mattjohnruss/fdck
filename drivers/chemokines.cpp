#include <problem.h>
#include <config_file.h>

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
        dx_(1.0/(n_node-1)),
        cn_theta_(1.0),
        time_factor_(1.0),
        c_u_offset_(0),
        c_b_offset_(n_node_),
        c_s_offset_(2*n_node_),
        phi_offset_(3*n_node_)
    {
        std::cout << "n_node = " << n_node_ << '\n';
        std::cout << "n_dof  = " << n_dof_ << '\n';
        std::cout << "dx     = " << dx_ << '\n';
        std::cout << "dt     = " << dt_ << "\n\n";
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // set zero initial conditions
        for(unsigned i = 0; i < n_node_; ++i)
        {
            c_u(0,i) = 0.0;
            c_b(0,i) = 0.0;
            c_s(0,i) = 0.0;
            phi(0,i) = 0.0;
        }

        c_u(0,0) = 1.0;
        phi(0,n_node_-1) = 1.0;
    }

    ChemokinesParams p;

    const double c_u(unsigned t, unsigned i) const
    {
        return u(t, i + c_u_offset_);
    }

    const double c_b(unsigned t, unsigned i) const
    {
        return u(t, i + c_b_offset_);
    }

    const double c_s(unsigned t, unsigned i) const
    {
        return u(t, i + c_s_offset_);
    }

    const double phi(unsigned t, unsigned i) const
    {
        return u(t, i + phi_offset_);
    }

    double& c_u(unsigned t, unsigned i)
    {
        return u(t, i + c_u_offset_);
    }

    double& c_b(unsigned t, unsigned i)
    {
        return u(t, i + c_b_offset_);
    }

    double& c_s(unsigned t, unsigned i)
    {
        return u(t, i + c_s_offset_);
    }

    double& phi(unsigned t, unsigned i)
    {
        return u(t, i + phi_offset_);
    }

    const double d_c_u_dx(unsigned t, unsigned i) const
    {
        return d_u_dx_helper(t, i, c_u_offset_);
    }

    const double d_c_b_dx(unsigned t, unsigned i) const
    {
        return d_u_dx_helper(t, i, c_b_offset_);
    }

    const double d_c_s_dx(unsigned t, unsigned i) const
    {
        return d_u_dx_helper(t, i, c_s_offset_);
    }

    const double d_phi_dx(unsigned t, unsigned i) const
    {
        return d_u_dx_helper(t, i, phi_offset_);
    }

    void make_steady() override
    {
        Problem::make_steady();
        time_factor_ = 0.0;
        old_cn_theta_ = cn_theta_;
        cn_theta_ = 1.0;
    }

    void make_unsteady() override
    {
        Problem::make_unsteady();
        time_factor_ = 1.0;
        cn_theta_ = old_cn_theta_;
    }

    void output(std::ostream &out) const override
    {
        for(unsigned i = 0; i < n_node_; ++i)
        {
            // Output values
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' '
                << c_u(0,i) << ' '
                << c_b(0,i) << ' '
                << c_s(0,i) << ' '
                << phi(0,i) << ' ';

            // Output derivatives
            out << d_c_u_dx(0,i) << ' '
                << d_c_b_dx(0,i) << ' '
                << d_c_s_dx(0,i) << ' '
                << d_phi_dx(0,i) << '\n';
        }
    }

private:
    const unsigned n_node_;
    const double dx_;
    double cn_theta_;
    double time_factor_;
    double old_cn_theta_;

    const unsigned c_u_offset_;
    const unsigned c_b_offset_;
    const unsigned c_s_offset_;
    const unsigned phi_offset_;

    const double d_u_dx_helper(unsigned t, unsigned i, unsigned offset) const
    {
        if(i == 0)
        {
            // Left boundary
            return (-1.5*u(t,offset+0) + 2.0*u(t,offset+1) - 0.5*u(t,offset+2))/dx_;
        }
        else if(i == (n_node_-1))
        {
            // Right boundary
            return (0.5*u(t,offset+n_node_-3) - 2.0*u(t,offset+n_node_-2) + 1.5*u(t,offset+n_node_-1))/dx_;
        }
        else
        {
            // Sanity check (debug only)
            assert(i > 0 && i < (n_node_-1));

            // Bulk
            return (-0.5*u(t,offset+i-1) + 0.5*u(t,offset+i+1))/dx_;
        }
    }

    void calculate_residual()
    {
        // set all entries to zero
        // not neccessary here since all entries are set explicitly
        residual_.setZero();

        // LHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u
        residual_(c_u_offset_+0) += (c_u(0,0) - 1.0)*dx_*dx_;

        // C_b (full equation, not just BC)
        // time derivative
        residual_(c_b_offset_+0) += -time_factor_*(c_b(0,0) - c_b(1,0))*dx_*dx_/dt_;

        // binding/unbinding (Crank-Nicolson)
        residual_(c_b_offset_+0) += cn_theta_*(p.alpha*c_u(0,0) - p.beta*c_b(0,0) - p.gamma_b*c_b(0,0)*phi(0,0))*dx_*dx_;
        residual_(c_b_offset_+0) += (1.0-cn_theta_)*(p.alpha*c_u(1,0) - p.beta*c_b(1,0) - p.gamma_b*c_b(1,0)*phi(1,0))*dx_*dx_;

        // C_s
        residual_(c_s_offset_+0) += (-1.5*c_s(0,0) + 2.0*c_s(0,1) - 0.5*c_s(0,2))*dx_ - p.p_u*c_s(0,0)*dx_*dx_;

        // phi
        residual_(phi_offset_+0) += (-1.5*phi(0,0) + 2.0*phi(0,1) - 0.5*phi(0,2))*dx_ - p.lambda*phi(0,0)*dx_*dx_;

        // Bulk equations
        // --------------------------------------------------------------------

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // C_u
            // time derivative
            residual_(c_u_offset_+i) += -time_factor_*(c_u(0,i) - c_u(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(c_u_offset_+i) += cn_theta_*(c_u(0,i-1) - 2.0*c_u(0,i) + c_u(0,i+1));
            residual_(c_u_offset_+i) += (1.0-cn_theta_)*(c_u(1,i-1) - 2.0*c_u(1,i) + c_u(1,i+1));

            // advection (Crank-Nicolson)
            residual_(c_u_offset_+i) += -cn_theta_*p.p_u*(-0.5*c_u(0,i-1) + 0.5*c_u(0,i+1))*dx_;
            residual_(c_u_offset_+i) += -(1.0-cn_theta_)*p.p_u*(-0.5*c_u(1,i-1) + 0.5*c_u(1,i+1))*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_u_offset_+i) += cn_theta_*(-p.alpha*c_u(0,i) + p.beta*c_b(0,i) - p.gamma_u*c_u(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_u_offset_+i) += (1.0-cn_theta_)*(-p.alpha*c_u(1,i) + p.beta*c_b(1,i) - p.gamma_u*c_u(1,i)*phi(1,i))*dx_*dx_;

            // C_b
            // time derivative
            residual_(c_b_offset_+i) += -time_factor_*(c_b(0,i) - c_b(1,i))*dx_*dx_/dt_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_b_offset_+i) += cn_theta_*(p.alpha*c_u(0,i) - p.beta*c_b(0,i) - p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_b_offset_+i) += (1.0-cn_theta_)*(p.alpha*c_u(1,i) - p.beta*c_b(1,i) - p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;

            // C_s
            // time derivative
            residual_(c_s_offset_+i) += -time_factor_*(c_s(0,i) - c_s(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(c_s_offset_+i) += cn_theta_*p.D_su*(c_s(0,i-1) - 2.0*c_s(0,i) + c_s(0,i+1));
            residual_(c_s_offset_+i) += (1.0-cn_theta_)*p.D_su*(c_s(1,i-1) - 2.0*c_s(1,i) + c_s(1,i+1));

            // advection (Crank-Nicolson)
            residual_(c_s_offset_+i) += -cn_theta_*p.p_u*(-0.5*c_s(0,i-1) + 0.5*c_s(0,i+1))*dx_;
            residual_(c_s_offset_+i) += -(1.0-cn_theta_)*p.p_u*(-0.5*c_s(1,i-1) + 0.5*c_s(1,i+1))*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_s_offset_+i) += cn_theta_*(p.gamma_u*c_u(0,i)*phi(0,i) + p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_s_offset_+i) += (1.0-cn_theta_)*(p.gamma_u*c_u(1,i)*phi(1,i) + p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;

            // phi
            // time derivative
            residual_(phi_offset_+i) += -time_factor_*(phi(0,i) - phi(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(phi_offset_+i) += cn_theta_*p.D_ju*(phi(0,i-1) - 2.0*phi(0,i) + phi(0,i+1));
            residual_(phi_offset_+i) += (1.0-cn_theta_)*p.D_ju*(phi(1,i-1) - 2.0*phi(1,i) + phi(1,i+1));

            // advection (Chemotaxis?) (Crank-Nicolson)
            residual_(phi_offset_+i) += -cn_theta_*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1))*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1));
            residual_(phi_offset_+i) += -(1.0-cn_theta_)*p.nu*(-0.5*phi(1,i-1) + 0.5*phi(1,i+1))*(-0.5*c_b(1,i-1) + 0.5*c_b(1,i+1));

            residual_(phi_offset_+i) += -cn_theta_*p.nu*phi(0,i)*(c_b(0,i-1) - 2.0*c_b(0,i) + c_b(0,i+1));
            residual_(phi_offset_+i) += -(1.0-cn_theta_)*p.nu*phi(1,i)*(c_b(1,i-1) - 2.0*c_b(1,i) + c_b(1,i+1));
        }

        // RHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u
        residual_(c_u_offset_+n_node_-1) += (c_u(0,n_node_-1) - 0.0)*dx_*dx_;

        // C_b (full equation, not just BC)
        // time derivative
        residual_(c_b_offset_+n_node_-1) += -time_factor_*(c_b(0,n_node_-1) - c_b(1,n_node_-1))*dx_*dx_/dt_;

        // binding/unbinding (Crank-Nicolson)
        residual_(c_b_offset_+n_node_-1) += cn_theta_*(p.alpha*c_u(0,n_node_-1) - p.beta*c_b(0,n_node_-1) - p.gamma_b*c_b(0,n_node_-1)*phi(0,n_node_-1))*dx_*dx_;
        residual_(c_b_offset_+n_node_-1) += (1.0-cn_theta_)*(p.alpha*c_u(1,n_node_-1) - p.beta*c_b(1,n_node_-1) - p.gamma_b*c_b(1,n_node_-1)*phi(1,n_node_-1))*dx_*dx_;

        // C_s
        residual_(c_s_offset_+n_node_-1) += (0.5*c_s(0,n_node_-3) - 2.0*c_s(0,n_node_-2) + 1.5*c_s(0,n_node_-1))*dx_ - p.p_u*c_s(0,n_node_-1)*dx_*dx_;

        // phi
        residual_(phi_offset_+n_node_-1) += (phi(0,n_node_-1) - 1.0)*dx_*dx_;
    }

    void calculate_jacobian()
    {
        // set all entries to zero
        //jacobian_.setZero();

        std::vector<T> triplet_list;
        triplet_list.reserve(10 + 20*(n_node_-2) + 8);

        // LHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u / C_u
        triplet_list.push_back( T(c_u_offset_+0, c_u_offset_+0, 1.0*dx_*dx_) );

        // C_b / C_u
        triplet_list.push_back( T(c_b_offset_+0, c_u_offset_+0, cn_theta_*p.alpha*dx_*dx_) );

        // C_b / C_b
        triplet_list.push_back( T(c_b_offset_+0, c_b_offset_+0, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,0))*dx_*dx_) );

        // C_b / C_s
        // (no entries)

        // C_b / phi
        triplet_list.push_back( T(c_b_offset_+0, phi_offset_+0, -cn_theta_*(p.gamma_b*c_b(0,0))*dx_*dx_) );

        // C_s / C_s
        triplet_list.push_back( T(c_s_offset_+0, c_s_offset_+0, -1.5*dx_ - p.p_u*dx_*dx_) );
        triplet_list.push_back( T(c_s_offset_+0, c_s_offset_+1, 2.0*dx_) );
        triplet_list.push_back( T(c_s_offset_+0, c_s_offset_+2, -0.5*dx_) );

        // phi / phi
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+0, -1.5*dx_ - p.lambda*dx_*dx_) );
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+1, 2.0*dx_) );
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+2, -0.5*dx_) );

        // Bulk equations
        // --------------------------------------------------------------------

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // C_u / C_u
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-1,  cn_theta_*(1.0 + p.p_u*0.5*dx_)) );
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i,   -time_factor_*dx_*dx_/dt_ - cn_theta_*(2.0 - (p.alpha + p.gamma_u*phi(0,i))*dx_*dx_)) );
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+1,  cn_theta_*(1.0 - p.p_u*0.5*dx_)) );

            // C_u / C_b
            triplet_list.push_back( T(c_u_offset_+i, c_b_offset_+i, cn_theta_*p.beta*dx_*dx_) );

            // C_u / C_s
            // (no entries)

            // C_u / phi
            triplet_list.push_back( T(c_u_offset_+i, phi_offset_+i, -cn_theta_*p.gamma_u*c_u(0,i)*dx_*dx_) );

            // C_b / C_u
            triplet_list.push_back( T(c_b_offset_+i, c_u_offset_+i, cn_theta_*p.alpha*dx_*dx_) );

            // C_b / C_b
            triplet_list.push_back( T(c_b_offset_+i, c_b_offset_+i, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,i))*dx_*dx_) );

            // C_b / C_s
            // (no entries)

            // C_b / phi
            triplet_list.push_back( T(c_b_offset_+i, phi_offset_+i, -cn_theta_*(p.gamma_b*c_b(0,i))*dx_*dx_) );

            // C_s / C_u
            triplet_list.push_back( T(c_s_offset_+i, c_u_offset_+i, cn_theta_*(p.gamma_u*phi(0,i))*dx_*dx_) );

            // C_s / C_b
            triplet_list.push_back( T(c_s_offset_+i, c_b_offset_+i, cn_theta_*(p.gamma_b*phi(0,i)*dx_*dx_)) );

            // C_s / C_s
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-1, cn_theta_*(p.D_su + p.p_u*0.5*dx_)) );
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i,   -time_factor_*dx_*dx_/dt_ - cn_theta_*2.0*p.D_su) );
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+1, cn_theta_*(p.D_su - p.p_u*0.5*dx_)) );

            // C_s / phi
            triplet_list.push_back( T(c_s_offset_+i, phi_offset_+i, cn_theta_*(p.gamma_u*c_u(0,i) + p.gamma_b*c_b(0,i))*dx_*dx_) );

            // phi / C_u
            // (no entries)

            // phi / C_b
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i-1,  cn_theta_*( 0.5*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1)) - p.nu*phi(0,i) )) );
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i,    cn_theta_*2.0*p.nu*phi(0,i)) );
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i+1, -cn_theta_*( 0.5*p.nu*(-0.5*phi(0,i-1) + 0.5*phi(0,i+1)) + p.nu*phi(0,i) )) );

            // phi / C_s
            // (no entries)

            // phi / phi
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-1, cn_theta_*(p.D_ju + p.nu*0.5*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1)))) );
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i,   -time_factor_*dx_*dx_/dt_ + cn_theta_*(-2.0*p.D_ju - p.nu*(c_b(0,i-1) - 2.0*c_b(0,i) + c_b(0,i+1)))) );
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+1, cn_theta_*(p.D_ju - p.nu*0.5*(-0.5*c_b(0,i-1) + 0.5*c_b(0,i+1)))) );
        }

        // RHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u / C_u
        triplet_list.push_back( T(c_u_offset_+n_node_-1, c_u_offset_+n_node_-1, 1.0*dx_*dx_) );

        // C_b / C_u
        triplet_list.push_back( T(c_b_offset_+n_node_-1, c_u_offset_+n_node_-1, cn_theta_*p.alpha*dx_*dx_) );

        // C_b / C_b
        triplet_list.push_back( T(c_b_offset_+n_node_-1, c_b_offset_+n_node_-1, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,n_node_-1))*dx_*dx_) );

        // C_b / C_s
        // (no entries)

        // C_b / phi
        triplet_list.push_back( T(c_b_offset_+n_node_-1, phi_offset_+n_node_-1, -cn_theta_*(p.gamma_b*c_b(0,n_node_-1))*dx_*dx_) );

        // C_s / C_s
        triplet_list.push_back( T(c_s_offset_+n_node_-1, c_s_offset_+n_node_-3,  0.5*dx_) );
        triplet_list.push_back( T(c_s_offset_+n_node_-1, c_s_offset_+n_node_-2, -2.0*dx_) );
        triplet_list.push_back( T(c_s_offset_+n_node_-1, c_s_offset_+n_node_-1,  1.5*dx_ - p.p_u*dx_*dx_) );

        // phi / phi
        triplet_list.push_back( T(phi_offset_+n_node_-1, phi_offset_+n_node_-1, dx_*dx_) );

        // construct matrix
        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }
};

int main(int argc, char **argv)
{
    if(argc != 5 && argc != 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " config_file n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[2]);
    const double dt = std::atof(argv[3]);
    const double t_max = std::atof(argv[4]);

    unsigned output_interval = 1;

    if(argc == 6)
    {
        output_interval = std::atoi(argv[5]);
    }

    ChemokinesProblem1D problem(n_node, dt);
    ChemokinesProblem1D::Max_residual = 1e-14;

    std::ifstream config_file(argv[1]);
    ConfigFile cf(config_file);

    problem.p.p_u     = cf.get<double>("p_u");
    problem.p.alpha   = cf.get<double>("alpha");
    problem.p.beta    = cf.get<double>("beta");
    problem.p.gamma_u = cf.get<double>("gamma_u");
    problem.p.gamma_b = cf.get<double>("gamma_b");
    problem.p.D_su    = cf.get<double>("D_su");
    problem.p.D_ju    = cf.get<double>("D_ju");
    problem.p.nu      = cf.get<double>("nu");
    problem.p.lambda  = cf.get<double>("lambda");

    std::cout << "p_u     = " << problem.p.p_u     << '\n';
    std::cout << "alpha   = " << problem.p.alpha   << '\n';
    std::cout << "beta    = " << problem.p.beta    << '\n';
    std::cout << "gamma_u = " << problem.p.gamma_u << '\n';
    std::cout << "gamma_b = " << problem.p.gamma_b << '\n';
    std::cout << "D_su    = " << problem.p.D_su    << '\n';
    std::cout << "D_ju    = " << problem.p.D_ju    << '\n';
    std::cout << "nu      = " << problem.p.nu      << '\n';
    std::cout << "lambda  = " << problem.p.lambda  << '\n';

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    bool dump = false;
    //bool dump = true;

    problem.Max_newton_iterations = 100;

    // perform a steady solve and output it
    problem.steady_solve(dump);
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    // timestepping loop
    while(problem.time() <= t_max)
    {
        // solve for current timestep
        problem.unsteady_solve(dump);

        if(i % output_interval == 0)
        {
            // output current solution
            //std::cout << "Outputting solution at time = " << problem.time() << '\n';
            std::cout << ";\tOutputting";
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }
}
