#include "advection_diffusion_reaction_problem.h"
#include "config.h"
#include "functions.h"
#include "log.h"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <fstream>

#include <fenv.h>

using namespace fdck;

struct KellerSegelParams
{
    double rho_init;
    double c_init;
    double D;
    double r;
    std::unique_ptr<DifferentiableFunction> chi;
};

enum Variable
{
    rho = 0,
    c = 1,
};

class KellerSegelProblem : public AdvectionDiffusionReactionProblem
{
public:
    KellerSegelProblem(const unsigned n_node, const double L) :
        AdvectionDiffusionReactionProblem(2, n_node, 0, 0.0, L)
    {
        enable_bc(Boundary::Left,  { rho, c });
        enable_bc(Boundary::Right, { rho, c });

        enable_spatial_terms({ rho, c });

        set_variable_names({ "rho", "c" });

        Max_residual = 1.0e-14;
        Max_newton_iterations = 50;

        trace_file_.open("trace.dat");
        trace_header_ = "t rho_total c_total";
        trace_file_ << trace_header_ << '\n';
    }

    ~KellerSegelProblem()
    {
    }

    void set_initial_conditions()
    {
        clear_solution();

        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<double> dist{0.0, 0.001};

        for(unsigned i = 0; i < n_node_; ++i)
        {
            //double x = this->x(i);

            //u(rho, i) = 100.0*std::exp(-100.0*x*x);
            //u(c, i)   = 50.0*std::exp(-50.0*x*x);

            u(rho, i) = p.rho_init + dist(gen);
            u(c, i)   = p.c_init + dist(gen);
        }
    }

    void trace()
    {
        trace_file_ << time() << " "
                    << integrate_solution(rho) << " "
                    << integrate_solution(c) << '\n';
    }

    KellerSegelParams p;

private:
    std::string trace_header_;
    std::ofstream trace_file_;

    void actions_after_timestep() override
    {
        trace();
    }

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        // No flux conditions on both boundaries

        // set the node number i according to the boundary
        unsigned i = 0;

        if(b == Boundary::Left)
        {
            i = 0;
        }
        else if(b == Boundary::Right)
        {
            i = n_node_-1;
        }

        // get the rho velocity at the boundary
        double v_rho = 0.0;
        get_v(i, rho, v_rho);

        // get the diffusivities at the boundary
        std::vector<double> d(n_var_);
        get_d(0, i, d);

        a1[rho] = v_rho; a2[rho] = -d[rho]; a3[rho] = 0.0;
        a1[c]   = 0.0;   a2[c]   = -d[c];   a3[c]   = 0.0;
    }

    // NOTE:
    // Use default FD impl of get_dbc_du

    void get_d(const unsigned,
               const unsigned,
               std::vector<double> &d) const override
    {
        d[rho] = p.D;
        d[c]   = 1.0;
    }

    void get_dd_du(const unsigned,
                   const unsigned,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        dd_du[rho][rho] = 0.0;
        dd_du[rho][c] = 0.0;
        dd_du[c][rho] = 0.0;
        dd_du[c][c] = 0.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == rho)
        {
            double dc_dx = 0.0;

            for(const auto& [j, w] : central_1_stencil_weights(i))
            {
                dc_dx += w*u(c, i+j)/dx_;
            }

            v = p.chi->value(u(rho, i))*dc_dx;
        }
    }

    // NOTE:
    // Use default FD impl of get_dv_du

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[rho] = p.r*u[rho]*(1.0 - u[rho]);
        r[c]   = - u[c] + u[rho];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[rho][rho] = p.r*(1.0 - 2.0*u[rho]);
        dr_du[rho][c] = 0.0;
        dr_du[c][rho] = 1.0;
        dr_du[c][c] = -1.0;
    }
};

int main(int argc, char **argv)
{
    //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0]
                  << " config_file [ --param1 value1 --param2 value 2 ... ]\n";
        std::exit(1);
    }

    std::ifstream config_file(argv[1]);
    Config cf;
    cf.parse_config_file(config_file);
    cf.parse_command_line(argc, argv);
    cf.print_all();

    const double L = cf.get_or<double>("L", 1.0);
    const unsigned n_node = cf.get_or<unsigned>("n_node", 101);
    const double dt = cf.get_or<double>("dt", 0.01);
    const double t_max = cf.get_or<double>("t_max", 10);
    const unsigned output_interval = cf.get_or<unsigned>("output_interval", 1);

    log::set_level("trace");

    KellerSegelProblem problem(n_node, L);

    problem.enable_terse_logging();

    problem.p.rho_init = cf.get<double>("rho_init");
    problem.p.c_init = cf.get<double>("c_init");

    problem.p.D = cf.get<double>("D");
    problem.p.r = cf.get<double>("r");

    if(std::string chi_type = cf.get<std::string>("chi"); chi_type == "constant" || chi_type == "const")
    {
        problem.p.chi =
            std::make_unique<ConstantFunction>(cf.get<double>("chi_const_val"));
    }

    char filename[200];
    std::ofstream outfile;

    const bool do_time_evolution = cf.get_or<bool>("time_evo", true);
    const bool do_steady_solve = cf.get_or<bool>("steady", true);

    if(do_time_evolution)
    {
        FDCK_INFO("Time evolution:");

        problem.enable_exit_on_solve_fail();

        problem.set_initial_conditions();

        problem.trace();

        std::sprintf(filename, "output_%05u.csv", 0);
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        unsigned i = 1;

        while(problem.time() < t_max)
        {
            problem.unsteady_solve(dt);

            if(i % output_interval == 0)
            {
                FDCK_TRACE("Outputting");
                std::sprintf(filename, "output_%05u.csv", i/output_interval);
                outfile.open(filename);
                problem.output(outfile);
                outfile.close();
            }

            ++i;
        }

        FDCK_INFO("Reached t > t_max ({}) after performing {} timesteps", t_max, i-1);
    }

    if(do_steady_solve)
    {
        FDCK_INFO("Steady solve:");

        if(do_time_evolution == false)
        {
            problem.set_initial_conditions();
        }

        problem.steady_solve();
        std::sprintf(filename, "output_steady.csv");
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();
    }

    return 0;
}
