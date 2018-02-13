#include <advection_diffusion_reaction_problem.h>
#include <config_file.h>

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

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(4, n_node, dt)
    {
        enable_bc(Boundary::Left, { 0, 2, 3 });
        enable_bc(Boundary::Right, { 0, 2, 3 });
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Then set c_u and phi at the appropriate boundaries
        u(0, 0, 0) = 1.0;
        u(0, 3, n_node_-1) = 1.0;
    }

    ChemokinesParams p;

private:
    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[0] = 1.0;       a2[0] = 0.0; a3[0] = 1.0;
            a1[2] = 1.0;       a2[2] = 0.0; a3[2] = 0.0;
            a1[3] = -p.lambda; a2[3] = 1.0; a3[3] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[0] = 1.0; a2[0] = 0.0; a3[0] = 0.0;
            a1[2] = 1.0; a2[2] = 0.0; a3[2] = 0.0;
            a1[3] = 1.0; a2[3] = 0.0; a3[3] = 1.0;
        }
    }

    void get_d(std::vector<double> &d) const override
    {
        d[0] = 1.0;
        d[1] = 0.0;
        d[2] = p.D_su;
        d[3] = p.D_ju;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        if(var == 0)
        {
            v = p.p_u;
        }
        else if(var == 1)
        {
            v = 0.0;
        }
        else if(var == 2)
        {
            v = p.p_u;
        }
        else if(var == 3)
        {
            std::unordered_map<int, double> weights;

            if(i == 0)
            {
                weights = stencil::forward_1::weights;
            }
            else if(i == n_node_-1)
            {
                weights = stencil::backward_1::weights;
            }
            else
            {
                weights = stencil::central_1::weights;
            }

            v = 0;

            for(const auto& [j, w] : weights)
            {
                v += w*p.nu*u(0, 1, i+j);
            }
        }
    }

    void get_dv_du(const unsigned i,
                   const std::vector<double> &u,
                   std::vector<std::vector<double>> &dv_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dv_du[var][var2] = 0.0;
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[0] = -p.alpha*u[0] + p.beta*u[1] - p.gamma_u*u[3]*u[0];
        r[1] = p.alpha*u[0] - p.beta*u[1] - p.gamma_b*u[3]*u[1];
        r[2] = p.gamma_u*u[3]*u[0] + p.gamma_b*u[3]*u[1];
        r[4] = 0.0;
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dr_du[var][var2] = 0.0;
            }
        }
    }
};

int main(int argc, char **argv)
{
    unsigned n_node = 11;

    if(argc == 3)
    {
        n_node = std::atoi(argv[2]);
    }

    ChemokinesProblem1D problem(n_node, 0.01);

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

    problem.enable_fd_jacobian();
    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve and output it
    problem.steady_solve();
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    // set initial conditions (zero)
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    double t_max = 10.0;
    unsigned output_interval = 1;

    // timestepping loop
    while(problem.time() <= t_max)
    {
        // solve for current timestep
        problem.unsteady_solve();

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
