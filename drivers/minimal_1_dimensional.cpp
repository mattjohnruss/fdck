#include <advection_diffusion_reaction_problem.h>
#include <config.h>
#include <utilities.h>

#include <fstream>

using namespace mjrfd;

struct ChemokinesParams
{
    double D;
    double u;
};

enum Variable
{
    c_u = 0,
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(1, n_node, 0, 0.0, 0.9217),
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
                           -0.462170292297035 },
        start_time_(30.0),
        end_time_(2730.0)
    {
        enable_bc(Boundary::Left,  { c_u });
        enable_bc(Boundary::Right, { c_u });

        enable_spatial_terms({ c_u });

        disable_output_time_column();

        set_variable_names({ "c_u" });

        Max_residual = 1.0e-14;

        //trace_file_.open("trace.dat");
        //trace_header_ =
            //"t a3_{left} a3_{right} c_u(0) c_u(1)";
        //trace_file_ << trace_header_ << '\n';
    }

    ~ChemokinesProblem1D()
    {
        //trace_file_.close();
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Initial time is 30s
        time() = start_time_;

        // Set c_u to 1.0 at the left-hand boundary
        //u(0, c_u, 0) = 1.0;

        // Set c_u to the BC functions evalutated at t = 30
        u(c_u, 0) =
            utilities::evaluate_polynomial(map_time(time()), inlet_poly_coeffs);
        u(c_u, n_node_-1) =
            utilities::evaluate_polynomial(map_time(time()), outlet_poly_coeffs);
    }

    ChemokinesParams p;

    // Coeffs for the polynomial fits of the inlet/outlet BC functions
    const std::vector<double> inlet_poly_coeffs;
    const std::vector<double> outlet_poly_coeffs;

private:
    //std::string trace_header_;
    //std::ofstream trace_file_;

    double start_time_;
    double end_time_;

    double map_time(double time) const
    {
        return (time - start_time_)/(end_time_ - start_time_);
    }

    //void actions_after_timestep() override
    //{
        //std::vector<double> a1_left(n_var_);
        //std::vector<double> a2_left(n_var_);
        //std::vector<double> a3_left(n_var_);

        //std::vector<double> a1_right(n_var_);
        //std::vector<double> a2_right(n_var_);
        //std::vector<double> a3_right(n_var_);

        //get_bc(Boundary::Left, a1_left, a2_left, a3_left);
        //get_bc(Boundary::Right, a1_right, a2_right, a3_right);

        ////trace_file_ << time() << " "
                    ////<< a3_left[c_u] << " "
                    ////<< a3_right[c_u] << " "
                    ////<< u(c_u, 0) << " "
                    ////<< u(c_u, n_node_-1) << '\n'
                    ////<< std::flush;
    //}

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c_u] = 1.0;
            a2[c_u] = 0.0;
            a3[c_u] = utilities::evaluate_polynomial(map_time(time()), inlet_poly_coeffs);
        }
        if(b == Boundary::Right)
        {
            a1[c_u] = 1.0;
            a2[c_u] = 0.0;
            a3[c_u] = utilities::evaluate_polynomial(map_time(time()), outlet_poly_coeffs);
        }
    }

    void get_d(const unsigned,
               const unsigned,
               std::vector<double> &d) const override
    {
        d[c_u] = p.D;
    }

    void get_v(const unsigned,
               const unsigned,
               double &v) const override
    {
        v = p.u;
    }

    //void get_dv_du(const unsigned i,
                   //const unsigned var,
                   //const unsigned i2,
                   //const unsigned var2,
                   //double &dv_du) const override
    //{
        //dv_du = 0.0;
    //}

    void get_r(const std::vector<double> &,
               std::vector<double> &r) const override
    {
        r[c_u] = 0.0;
    }

    //void get_dr_du(const std::vector<double> &u,
                   //std::vector<std::vector<double>> &dr_du) const override
    //{
        //dr_du[c_u][c_u] = 0.0;
    //}
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

    problem.p.u = cf.get<double>("u");
    problem.p.D = cf.get<double>("D");

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
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

#ifdef _WIN32
    std::cout.sync_with_stdio(false);
#endif

    // timestepping loop
    while(problem.time() < t_max)
    {
        // solve for current timestep
        problem.unsteady_solve(dt);

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

    std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
}
