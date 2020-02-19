#include "advection_diffusion_reaction_problem.h"
#include "config.h"

#include <fstream>
#include <sstream>

using namespace fdck;

struct ChemokinesParams
{
    double D_iu;
    double phi_i_init;
    double R;
    double M;
    double J_i_right;
};

enum Variable
{
    phi_i = 0,
};

class ImmatureDCsProblem : public AdvectionDiffusionReactionProblem
{
public:
    ImmatureDCsProblem(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(1, n_node)
    {
        enable_bc(Boundary::Left,  { phi_i });
        enable_bc(Boundary::Right, { phi_i });

        enable_spatial_terms({ phi_i });

        set_variable_names({ "phi_i" });

        Max_residual = 1.0e-14;
    }

    ~ImmatureDCsProblem()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Set phi_i to its uniform initial density
        for(unsigned i = 0; i < n_node_; ++i)
        {
            u(0, phi_i, i) = p.phi_i_init;
        }
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
            a1[phi_i] = 0.0; a2[phi_i] = 1.0; a3[phi_i] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[phi_i] = 0.0; a2[phi_i] = 1.0; a3[phi_i] = -p.J_i_right;
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[phi_i] = p.D_iu;
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
        r[phi_i] = p.R*u[phi_i]*(1.0 - u[phi_i]) - p.M*u[phi_i];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[phi_i][phi_i] = p.R*(1.0 - 2.0*u[phi_i]) - p.M;
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

    ImmatureDCsProblem problem(n_node);

    std::ifstream config_file(argv[1]);
    Config cf(config_file);
    cf.print_all();

    problem.p.D_iu       = cf.get<double>("D_iu");
    problem.p.phi_i_init = cf.get<double>("phi_i_init");
    problem.p.R          = cf.get<double>("R");
    problem.p.M          = cf.get<double>("M");
    problem.p.J_i_right  = cf.get<double>("J_i_right");

    problem.enable_terse_logging();

    //problem.enable_fd_jacobian();

    char filename[200];
    std::ofstream outfile;
    std::stringstream outstream;

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    //std::sprintf(filename, "output_%05u.csv", 0);
    //outfile.open(filename);
    //problem.output(outfile);
    //outfile.close();
    problem.output(outstream);

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
            //std::sprintf(filename, "output_%05u.csv", i/output_interval);
            //outfile.open(filename);
            //problem.output(outfile);
            //outfile.close();
            outstream << '\n';
            problem.output(outstream);
        }

        ++i;
    }

    outfile.open("output_all.csv");
    outfile << outstream.rdbuf();
    outfile.close();

    std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";

    // set initial conditions
    //problem.set_initial_conditions();

    // perform a steady solve and output it
    problem.steady_solve();
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    std::cout << '\n';

}
