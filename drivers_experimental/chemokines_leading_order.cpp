#include "advection_diffusion_reaction_problem.h"
#include "config.h"

#include <fstream>

using namespace fdck;

struct ChemokinesParams
{
    double pe_u;
    double alpha;
    double beta;
    double gamma;
    double nu_u;
    double nu_b;
    double R;
    double M;
    double phi_m0_left;
    double phi_m0_right;
};

enum Variable
{
    c_u0   = 0,
    phi_m0 = 1
    //c_b   = 1,
    //c_s   = 2,
    //phi_i = 3,
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(2, n_node, dt)
    {
        enable_bc(Boundary::Left,  { c_u0, phi_m0 });
        enable_bc(Boundary::Right, { c_u0, phi_m0 });

        set_variable_names({ "c_u0", "phi_m0" });

        Max_residual = 1.0e-14;
        //Max_newton_iterations = 50;
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        //// Set c_u to 1.0 at the left-hand boundary
        //u(0, c_u, 0) = 1.0;

        //// Set phi_i to its uniform initial density
        //for(unsigned i = 0; i < n_node_; ++i)
        //{
            //u(0, phi_i, i) = p.phi_i_init;
        //}
    }

    ChemokinesParams p;

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
            // TODO proper flux bcs for phi_m0
            a1[c_u0]   = 1.0;       a2[c_u0]   = 0.0; a3[c_u0]   = 1.0;
            a1[phi_m0] = -p.lambda; a2[phi_m0] = 1.0; a3[phi_m0] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[c_u0]   = 1.0; a2[c_u0]   = 0.0; a3[c_u0]   = 0.0;
            a1[phi_m0] = 0.0; a2[phi_m0] = 1.0; a3[phi_m0] = 0.0;
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[c_u0]   = 1.0;
        d[phi_m0] = 0.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == c_u0)
            v = p.pe_u;
        else if(var == phi_m0)
        {
            for(const auto& [j, w] : stencil_1_helper(i))
            {
                v += w*p.nu_u*u(0, c_u0, i+j)/dx_;
                v += w*p.nu_b*u(0, c_b, i+j)/dx_;
                v += w*p.nu_s*u(0, c_s, i+j)/dx_;
            }
        }
    }

    void get_dv_du(const unsigned i,
                   const unsigned var,
                   const unsigned i2,
                   const unsigned var2,
                   double &dv_du) const override
    {
        dv_du = 0.0;

        if(var == phi_m0 && (var2 == c_u0 || var2 == c_b || var2 == c_s))
        {
            double nu;

            if(var2 == c_u0)
                nu = p.nu_u;
            else if(var2 == c_b)
                nu = p.nu_b;
            else if(var2 == c_s)
                nu = p.nu_s;
            else
                std::cerr << "var2 is not any of the expected values!\n";

            for(const auto& [j, w] : stencil_1_helper(i))
            {
                if(i+j == i2)
                {
                    // this condition will be true at most once per loop so we
                    // can break
                    dv_du = w*nu/dx_;
                    break;
                }
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c_u0]   = -p.alpha*u[c_u0] + p.beta*u[c_b] - p.gamma_ui*u[phi_i]*u[c_u0] - p.gamma_um*u[phi_m0]*u[c_u0] - p.q_u*u[phi_i]*u[c_u0];
        r[c_b]   = p.alpha*u[c_u0] - p.beta*u[c_b] - p.gamma_bi*u[phi_i]*u[c_b] - p.gamma_bm*u[phi_m0]*u[c_b] - p.q_b*u[phi_i]*u[c_b];
        r[c_s]   = p.gamma_ui*u[phi_i]*u[c_u0] + p.gamma_um*u[phi_m0]*u[c_u0] + p.gamma_bi*u[phi_i]*u[c_b] + p.gamma_bm*u[phi_m0]*u[c_b] - p.q_s*u[phi_i]*u[c_s];
        r[phi_i] = p.R*u[phi_i]*(1.0 - u[phi_i]) - p.M*u[phi_i];
        r[phi_m0] = p.M*u[phi_i];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_u0][c_u0]   = -p.alpha - p.gamma_ui*u[phi_i] - p.gamma_um*u[phi_m0] - p.q_u*u[phi_i];
        dr_du[c_u0][c_b]   = p.beta;
        dr_du[c_u0][c_s]   = 0.0;
        dr_du[c_u0][phi_i] = -p.gamma_ui*u[c_u0] - p.q_u*u[c_u0];
        dr_du[c_u0][phi_m0] = -p.gamma_um*u[c_u0];

        dr_du[c_b][c_u0]   = p.alpha;
        dr_du[c_b][c_b]   = -p.beta - p.gamma_bi*u[phi_i] - p.gamma_bm*u[phi_m0] - p.q_b*u[phi_i];
        dr_du[c_b][c_s]   = 0.0;
        dr_du[c_b][phi_i] = -p.gamma_bi*u[c_b] - p.q_b*u[c_b];
        dr_du[c_b][phi_m0] = -p.gamma_bm*u[c_b];

        dr_du[c_s][c_u0]   = p.gamma_ui*u[phi_i] + p.gamma_um*u[phi_m0];
        dr_du[c_s][c_b]   = p.gamma_bi*u[phi_i] + p.gamma_bm*u[phi_m0];
        dr_du[c_s][c_s]   = -p.q_s*u[phi_i];
        dr_du[c_s][phi_i] = p.gamma_ui*u[c_u0] + p.gamma_bi*u[c_b] - p.q_s*u[c_s];
        dr_du[c_s][phi_m0] = p.gamma_um*u[c_u0] + p.gamma_bm*u[c_b];

        dr_du[phi_i][c_u0]   = 0.0;
        dr_du[phi_i][c_b]   = 0.0;
        dr_du[phi_i][c_s]   = 0.0;
        dr_du[phi_i][phi_i] = p.R*(1.0 - 2.0*u[phi_i]) - p.M;
        dr_du[phi_i][phi_m0] = 0.0;

        dr_du[phi_m0][c_u0]   = 0.0;
        dr_du[phi_m0][c_b]   = 0.0;
        dr_du[phi_m0][c_s]   = 0.0;
        dr_du[phi_m0][phi_i] = p.M;
        dr_du[phi_m0][phi_m0] = 0.0;
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

    std::ifstream config_file(argv[1]);
    Config cf(config_file);
    cf.print_all();

    problem.p.pe_u       = cf.get<double>("pe_u");
    problem.p.alpha      = cf.get<double>("alpha");
    problem.p.beta       = cf.get<double>("beta");
    problem.p.gamma_ui   = cf.get<double>("gamma_ui");
    problem.p.gamma_bi   = cf.get<double>("gamma_bi");
    problem.p.gamma_um   = cf.get<double>("gamma_um");
    problem.p.gamma_bm   = cf.get<double>("gamma_bm");
    problem.p.q_u        = cf.get<double>("q_u");
    problem.p.q_b        = cf.get<double>("q_b");
    problem.p.q_s        = cf.get<double>("q_s");
    problem.p.D_su       = cf.get<double>("D_su");
    problem.p.D_iu       = cf.get<double>("D_iu");
    problem.p.D_mu       = cf.get<double>("D_mu");
    problem.p.nu_u       = cf.get<double>("nu_u");
    problem.p.nu_b       = cf.get<double>("nu_b");
    problem.p.nu_s       = cf.get<double>("nu_s");
    problem.p.lambda     = cf.get<double>("lambda");
    problem.p.phi_i_init = cf.get<double>("phi_i_init");
    problem.p.R          = cf.get<double>("R");
    problem.p.M          = cf.get<double>("M");

    const bool do_steady_solve = cf.get<bool>("steady");
    const bool do_time_evolution = cf.get<bool>("time_evo");

    problem.enable_terse_logging();
    //problem.disable_exit_on_solve_fail();

    //problem.enable_fd_jacobian();

    char filename[200];
    std::ofstream outfile;

    if(do_steady_solve)
    {
        std::cout << "\nSteady solve:\n";

        // set initial conditions - required since phi_i is sensitive to ICs even
        // at steady state
        problem.set_initial_conditions();

        // perform a steady solve and output it
        problem.steady_solve();
        std::sprintf(filename, "output_steady.csv");
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        std::cout << '\n';
    }

    if(do_time_evolution)
    {
        std::cout << "\nTime evolution:";

        problem.enable_exit_on_solve_fail();

        // set initial conditions again
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
            problem.unsteady_solve();

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

    return 0;
}

