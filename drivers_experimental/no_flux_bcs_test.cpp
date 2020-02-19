#include "advection_diffusion_reaction_problem.h"
#include "config.h"

#include <fstream>

using namespace fdck;

struct Params
{
    double v_0;
    double r_0;
};

enum Variable
{
    c_0 = 0,
};

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(1, n_node)
    {
        enable_bc(Boundary::Left,  { c_0 });
        enable_bc(Boundary::Right, { c_0 });

        enable_spatial_terms({ c_0 });

        set_variable_names({ "c_0" });

        Max_residual = 1.0e-14;
        Max_newton_iterations = 50;
    }

    ~TestProblem()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Set c_0 to 1.0 at the left-hand boundary
        u(0, c_0, 0) = 1.0;

        for(unsigned i = 0; i < n_node_; ++i)
        {
            u(0, c_0, i) = 1.0 - x(i);
        }
    }

    Params p;

private:
    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        unsigned i;

        if(b == Boundary::Left)
        {
            i = 0;
        }
        else if(b == Boundary::Right)
        {
            i = n_node_-1;
        }

        // no flux conditions

        // get the diffusivity (unfortunately all of them...)
        std::vector<double> d(n_var_);
        get_d(0, i, d);

        // get the velocity
        double v_0 = 0;
        get_v(0, c_0, v_0);

        // set the flux to zero (overall sign is irrelevant since flux is zero)
        a1[c_0] = v_0; a2[c_0] = -d[c_0]; a3[c_0] = 0.0;
    }

    void get_dbc_du(Boundary b,
                    const unsigned i2,
                    std::vector<std::vector<double>> &da1_du,
                    std::vector<std::vector<double>> &da2_du,
                    std::vector<std::vector<double>> &da3_du) const override
    {
        // first set all the derivatives to zero since this will the be case
        // most for most entries

        for(unsigned var = 0; var < n_var_; ++var)
        {
            std::fill(da1_du[var].begin(), da1_du[var].end(), 0.0);
            std::fill(da2_du[var].begin(), da2_du[var].end(), 0.0);
            std::fill(da3_du[var].begin(), da3_du[var].end(), 0.0);
        }

        // a1[phi_m] terms

        //unsigned i = 0;

        //if(b == Boundary::Left)
            //i = 0;
        //else if(b == Boundary::Right)
            //i = n_node_-1;

        //double dv_du_var2 = 0.0;
        //get_dv_du(i, phi_m, i2, var2, dv_du_var2);
        //da1_du[phi_m][var2] = dv_du_var2;
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[c_0] = 1.0;
    }

    void get_dd_du(const unsigned t,
                   const unsigned i,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        for(auto var : { c_0 })
        {
            std::fill(dd_du[var].begin(), dd_du[var].end(), 0.0);
        }
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        if(var == c_0)
            v = p.v_0;
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
        r[c_0] = p.r_0;
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_0][c_0] = 0.0;
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

    TestProblem problem(n_node);

    std::ifstream config_file(argv[1]);
    Config cf(config_file);
    cf.print_all();

    problem.p.v_0 = cf.get<double>("v_0");
    problem.p.r_0 = cf.get<double>("r_0");

    const bool do_steady_solve = cf.get<bool>("steady");
    const bool do_time_evolution = cf.get<bool>("time_evo");

    const bool fd_jacobian = cf.get<bool>("fd_jacobian");

    if(fd_jacobian == true)
    {
        problem.enable_fd_jacobian();
    }

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    if(do_time_evolution)
    {
        std::cout << "\nTime evolution:";

        problem.enable_exit_on_solve_fail();

        // set initial condition
        problem.set_initial_conditions();

        // output initial conditions
        std::sprintf(filename, "output_%05u.csv", 0);
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        unsigned i = 1;

        // timestepping loop
        while(problem.time() <= t_max)
        {
            // solve for current timestep
            problem.unsteady_solve(dt);

            if(i % output_interval == 0)
            {
                // output current solution
                //std::cout << "Outputting solution at time = " << problem.time() << '\n';
                std::cout << ";\tOutputting";
                std::sprintf(filename, "output_%05u.csv", i/output_interval);
                outfile.open(filename);
                problem.output(outfile);
                outfile.close();
            }

            ++i;
        }

        std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
    }

    if(do_steady_solve)
    {
        std::cout << "\nSteady solve:\n";

        problem.disable_exit_on_solve_fail();

        // set initial conditions again - required since phi_i is sensitive to
        // ICs even at steady state
        problem.set_initial_conditions();

        // perform a steady solve and output it
        problem.steady_solve();
        std::sprintf(filename, "output_steady.csv");
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        std::cout << '\n';
    }

    return 0;
}

