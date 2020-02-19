#include "advection_diffusion_reaction_problem.h"
#include "config.h"
#include "utilities.h"
#include "log.h"

#include <fstream>

using namespace fdck;

struct ChemokinesParams
{
    double D;
    double pe_u;
    double timescale_ratio;
};

enum Variable
{
    c_u = 0,
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(1, n_node),
        inlet_poly_coeffs{ 0.076787020335399,
                           0.254758928187113,
                           1.195706261658923,
                          -3.972963699044380,
                           4.612795695646865,
                          -1.882844223154888 },

        outlet_poly_coeffs{ 0.023250104348371,
                            0.089837404766442,
                           -0.513577991056639,
                            1.274529558134060,
                           -0.250287771102272,
                           -0.462170292297035 }
    {
        enable_bc(Boundary::Left,  { c_u });
        enable_bc(Boundary::Right, { c_u });

        enable_spatial_terms({ c_u });

        disable_output_time_column();

        set_variable_names({ "c_u" });

        Max_residual = 1.0e-14;
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        double c_start = utilities::evaluate_polynomial(0.0, inlet_poly_coeffs);
        double c_end   = utilities::evaluate_polynomial(0.0, outlet_poly_coeffs);

        for(unsigned i = 0; i < n_node_; ++i)
        {
            double m = (c_end - c_start)/(b_ - a_);
            double x = this->x(i);

            u(0, c_u, i) = m*x + c_start;

            FDCK_DEBUG("x = {}, u(x, 0) = {}", x, u(0, c_u, i));
        }
    }

    ChemokinesParams p;

    // Coeffs for the polynomial fits of the inlet/outlet BC functions
    const std::vector<double> inlet_poly_coeffs;
    const std::vector<double> outlet_poly_coeffs;

private:
    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c_u] = 1.0;
            a2[c_u] = 0.0;
            a3[c_u] = utilities::evaluate_polynomial(time(), inlet_poly_coeffs);
        }
        if(b == Boundary::Right)
        {
            a1[c_u] = 1.0;
            a2[c_u] = 0.0;
            a3[c_u] = utilities::evaluate_polynomial(time(), outlet_poly_coeffs);
        }
    }

    void get_d(const unsigned,
               const unsigned,
               std::vector<double> &d) const override
    {
        d[c_u] = 1.0;
    }

    void get_v(const unsigned,
               const unsigned,
               double &v) const override
    {
        v = p.pe_u;
    }

    void get_r(const std::vector<double> &,
               std::vector<double> &r) const override
    {
        r[c_u] = 0.0;
    }

    void get_timescale_ratio(double &timescale_ratio) const override
    {
        timescale_ratio = p.timescale_ratio;
    }
};

int main(int argc, char **argv)
{
    if(argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " config_file n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[2]);
    const double dt = std::atof(argv[3]);
    const double t_max = std::atof(argv[4]);

    unsigned output_interval = 1;

    if(argc >= 6)
    {
        output_interval = std::atoi(argv[5]);
    }

    std::ifstream config_file(argv[1]);
    Config cf;
    cf.parse_config_file(config_file);
    cf.parse_command_line(argc, argv);
    cf.print_all();

    log::set_level(cf.get_or<std::string>("log_level", "info"));

    ChemokinesProblem1D problem(n_node);

    problem.p.D               = cf.get_or<double>("D", 1.0);
    problem.p.pe_u            = cf.get<double>("pe_u");
    problem.p.timescale_ratio = cf.get_or<double>("timescale_ratio", 1.0);

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve and output it
    //problem.steady_solve();
    //std::sprintf(filename, "output_steady.csv");
    //outfile.open(filename);
    //problem.output(outfile);
    //outfile.close();

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05u.csv", 0);
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
        problem.unsteady_solve(dt);

        if(i % output_interval == 0)
        {
            // output current solution
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
