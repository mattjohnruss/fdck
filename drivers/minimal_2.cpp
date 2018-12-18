#include <advection_diffusion_reaction_problem.h>
#include <config.h>
#include <log.h>

#include <fstream>

using namespace mjrfd;

struct ChemokinesParams
{
    double pe_u;
    double alpha;
    double beta;
    double bc_time_factor;
};

enum Variable
{
    c_u = 0,
    c_b = 1,
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(2, n_node),
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

        set_variable_names({ "c_u", "c_b" });

        Max_residual = 1.0e-14;
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Set c_u to 1.0 at the left-hand boundary
        //u(0, c_u, 0) = 1.0;

        // Set c_u to the BC functions evalutated at t = 0
        //u(c_u, 0)         = evaluate_polynomial(0.0, inlet_poly_coeffs);
        //u(c_u, n_node_-1) = evaluate_polynomial(0.0, outlet_poly_coeffs);
        u(c_u, 0)         = 1.0;
        u(c_u, n_node_-1) = evaluate_polynomial(0.0, outlet_poly_coeffs)/evaluate_polynomial(0.0, inlet_poly_coeffs);
    }

    ChemokinesParams p;

    // Coeffs for the polynomial fits of the inlet/outlet BC functions
    const std::vector<double> inlet_poly_coeffs;
    const std::vector<double> outlet_poly_coeffs;

    // Evaluate a polynomial with coeffs in ascending power order at x using
    // Horner's method
    static double evaluate_polynomial(const double x,
                                      const std::vector<double> &coeffs)
    {
        double result = 0.0;
        unsigned n = coeffs.size();

        for(int i = n-1; i >= 0; --i)
        {
            assert(i <= (n-1) && i >= 0);
            result = result*x + coeffs[i];
        }

        return result;
    }

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
            a3[c_u] = evaluate_polynomial(p.bc_time_factor*time(), inlet_poly_coeffs)/evaluate_polynomial(0.0, inlet_poly_coeffs);
        }
        if(b == Boundary::Right)
        {
            a1[c_u] = 1.0;
            a2[c_u] = 0.0;
            a3[c_u] = evaluate_polynomial(p.bc_time_factor*time(), outlet_poly_coeffs)/evaluate_polynomial(0.0, inlet_poly_coeffs);
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[c_u] = 1.0;
        d[c_b] = 0.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == c_u)
            v = p.pe_u;
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
        r[c_u] = -p.alpha*u[c_u] + p.beta*u[c_b];
        r[c_b] =  p.alpha*u[c_u] - p.beta*u[c_b];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_u][c_u] = -p.alpha;
        dr_du[c_u][c_b] = p.beta;

        dr_du[c_b][c_u] = p.alpha;
        dr_du[c_b][c_b] = -p.beta;
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

    ChemokinesProblem1D problem(n_node);

    std::ifstream config_file(argv[1]);
    Config cf(config_file);
    cf.print_all();

    problem.p.pe_u  = cf.get<double>("pe_u");
    problem.p.alpha = cf.get<double>("alpha");
    problem.p.beta  = cf.get<double>("beta");
    problem.p.bc_time_factor = cf.get<double>("bc_time_factor");

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve and output it
    problem.steady_solve();
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
            MJRFD_TRACE("Outputting");
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }

    MJRFD_INFO("Reached t > t_max ({}) after performing {} timesteps", t_max, i-1);
}
