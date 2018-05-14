#include <advection_diffusion_reaction_problem.h>
#include <config_file.h>

#include <fstream>

using namespace mjrfd;

struct NonlinearDiffusionParams
{
    double D;
    double M;
};

enum Variable
{
    c = 0
};

class NonlinearDiffusionProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    NonlinearDiffusionProblem1D(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(1, n_node, dt)
    {
        enable_bc(Boundary::Left,  { c });
        enable_bc(Boundary::Right, { c });

        set_variable_names({ "c" });

        Max_residual = 1.0e-14;
        Max_newton_iterations = 100;
    }

    ~NonlinearDiffusionProblem1D()
    {
    }

    const double front_location(const double time) const
    {
        return std::cbrt(9*p.M*p.D*time);
    }

    void exact_solution(const double time,
                        const double x,
                        std::vector<double> &sol) const
    {
        const double X = front_location(time);

        if(x <= X)
        {
            sol[c] = (X*X - x*x)/(6.0*p.D*time);
        }
        else
        {
            sol[c] = 0.0;
        }
    }

    void set_initial_conditions()
    {
        // set to the exact solution at the nodes, at the current time
        std::vector<double> sol(1);

        for(unsigned i = 0; i < n_node_; ++i)
        {
            const double x = this->x(i);
            exact_solution(time(), x, sol);
            u(0, c, i) = sol[0];
        }
    }

    NonlinearDiffusionParams p;

private:
    const std::unordered_map<int, double>& stencil_1_helper(unsigned i) const
    {
        if(i == 0)
            return stencil::forward_1::weights;
        else if(i == n_node_-1)
            return stencil::backward_1::weights;
        else
            return stencil::central_1::weights;
    }

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c] = 0.0; a2[c] = 1.0; a3[c] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[c] = 1.0; a2[c] = 0.0; a3[c] = 0.0;
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[c] = p.D*u(t, c, i);
    }

    void get_dd_du(const unsigned t,
                   const unsigned i,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        dd_du[c][c] = p.D;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;
    }

    void get_dv_du(const unsigned i,
                   const unsigned var,
                   const unsigned i2,
                   const unsigned var2,
                   double &dv_du) const override
    {
        dv_du = 0.0;
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c] = 0.0;
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c][c] = 0.0;
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

    NonlinearDiffusionProblem1D problem(n_node, dt);

    std::ifstream config_file(argv[1]);
    ConfigFile cf(config_file);
    cf.print_all();

    problem.p.D = cf.get<double>("D");
    problem.p.M = cf.get<double>("M");

    problem.enable_terse_logging();

    bool fd_jacobian = cf.get<bool>("fd_jacobian");

    if(fd_jacobian)
    {
        std::cout << "using fd jacobian\n";
        problem.enable_fd_jacobian();
    }
    else
    {
        std::cout << "using exact jacobian\n";
    }

    char filename[200];
    std::ofstream outfile;

    // set the initial time to a nonzero value so that we can impose the exact
    // solution and then perform time-evolution from there
    problem.time() = 0.001;

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

#ifdef _WIN32
    std::cout.sync_with_stdio(false);
#endif

    // timestepping loop
    while(problem.time() <= t_max)
    {
        // solve for current timestep
        problem.unsteady_solve();

        if(i % output_interval == 0)
        {
            // output current solution
            std::cout << ";\tOutputting";
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            outfile << std::setprecision(16);
            problem.output(outfile);
            outfile.close();

            std::sprintf(filename, "output_exact_%05i.csv", i/output_interval);
            outfile.open(filename);
            outfile << std::setprecision(16);
            problem.output_exact(outfile);
            outfile.close();
        }

        ++i;
    }

    std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
}
