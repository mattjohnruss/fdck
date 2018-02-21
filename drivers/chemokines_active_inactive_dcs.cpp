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
    double phi_i_max;
};

enum Variable
{
    c_u   = 0,
    c_b   = 1,
    c_s   = 2,
    phi_a = 3,
    phi_i = 4
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(5, n_node, dt)
    {
        enable_bc(Boundary::Left,  { c_u, c_s, phi_a, phi_i });
        enable_bc(Boundary::Right, { c_u, c_s, phi_a, phi_i });

        set_variable_names({ "c_u", "c_b", "c_s", "phi_a", "phi_i" });

        Max_residual = 1.0e-14;
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Then set c_u and phi_a at the appropriate boundaries
        u(0, c_u, 0) = 1.0;
        u(0, phi_a, n_node_-1) = 1.0;
        u(0, phi_i, n_node_-1) = 1.0;
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
            a1[c_u]   = 1.0;       a2[c_u]   = 0.0; a3[c_u]   = 1.0;
            a1[c_s]   = 1.0;       a2[c_s]   = 0.0; a3[c_s]   = 0.0;
            a1[phi_a] = -p.lambda; a2[phi_a] = 1.0; a3[phi_a] = 0.0;
            a1[phi_i] = -p.lambda; a2[phi_i] = 1.0; a3[phi_i] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[c_u]   = 1.0; a2[c_u]   = 0.0; a3[c_u]   = 0.0;
            a1[c_s]   = 1.0; a2[c_s]   = 0.0; a3[c_s]   = 0.0;
            a1[phi_a] = 1.0; a2[phi_a] = 0.0; a3[phi_a] = 1.0;
            a1[phi_i] = 1.0; a2[phi_i] = 0.0; a3[phi_i] = 1.0;
        }
    }

    void get_d(std::vector<double> &d) const override
    {
        d[c_u]   = 1.0;
        d[c_b]   = 0.0;
        d[c_s]   = p.D_su;
        d[phi_a] = p.D_ju;
        d[phi_i] = p.D_ju;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == c_u)
            v = p.p_u;
        else if(var == c_s)
            v = p.p_u;
        else if(var == phi_a)
        {
            for(const auto& [j, w] : stencil_1_helper(i))
            {
                v += w*p.nu*u(0, c_b, i+j)/dx_;
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

        if(var == phi_a && var2 == c_b)
        {
            for(const auto& [j, w] : stencil_1_helper(i))
            {
                if(i+j == i2)
                {
                    // this condition will be true at most once per loop so we
                    // can break
                    dv_du = w*p.nu/dx_;
                    break;
                }
            }
        }
    }

    static double k(const double c_u)       { return 10.0*c_u; }
    static double dk_dc_u(const double c_u) { return 10.0; }

    static double k_a(const double c_u)       { return 5.0*c_u; }
    static double dk_a_dc_u(const double c_u) { return 5.0; }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c_u]   = -p.alpha*u[c_u] + p.beta*u[c_b] - p.gamma_u*u[phi_a]*u[c_u];
        r[c_b]   = p.alpha*u[c_u] - p.beta*u[c_b] - p.gamma_b*u[phi_a]*u[c_b];
        r[c_s]   = p.gamma_u*u[phi_a]*u[c_u] + p.gamma_b*u[phi_a]*u[c_b];
        r[phi_a] = k_a(u[c_u])*u[phi_i] - 1.0*u[phi_a];
        r[phi_i] = k(u[c_u])*u[phi_i]*(p.phi_i_max - u[phi_i]) - k_a(u[c_u])*u[phi_i];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_u][c_u]   = -p.alpha - p.gamma_u*u[phi_a];
        dr_du[c_u][c_b]   = p.beta;
        dr_du[c_u][c_s]   = 0.0;
        dr_du[c_u][phi_a] = -p.gamma_u*u[c_u];
        dr_du[c_u][phi_i] = 0.0;

        dr_du[c_b][c_u]   = p.alpha;
        dr_du[c_b][c_b]   = -p.beta - p.gamma_b*u[phi_a];
        dr_du[c_b][c_s]   = 0.0;
        dr_du[c_b][phi_a] = -p.gamma_b*u[c_b];
        dr_du[c_b][phi_i] = 0.0;

        dr_du[c_s][c_u]   = p.gamma_u*u[phi_a];
        dr_du[c_s][c_b]   = p.gamma_b*u[phi_a];
        dr_du[c_s][c_s]   = 0.0;
        dr_du[c_s][phi_a] = p.gamma_u*u[c_u] + p.gamma_b*u[c_b];
        dr_du[c_s][phi_i] = 0.0;

        dr_du[phi_a][c_u]   = dk_a_dc_u(u[c_u])*u[phi_i];
        dr_du[phi_a][c_b]   = 0.0;
        dr_du[phi_a][c_s]   = 0.0;
        dr_du[phi_a][phi_a] = -1.0;
        dr_du[phi_a][phi_i] = k_a(u[c_u]);

        dr_du[phi_i][c_u]   = dk_dc_u(u[c_u])*u[phi_i]*(p.phi_i_max - u[phi_i]) - dk_a_dc_u(u[c_u])*u[phi_i];
        dr_du[phi_i][c_b]   = 0.0;
        dr_du[phi_i][c_s]   = 0.0;
        dr_du[phi_i][phi_a] = 0.0;
        dr_du[phi_i][phi_i] = k(u[c_u])*(p.phi_i_max - 2.0*u[phi_i]) - k_a(u[c_u]);
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
    ConfigFile cf(config_file);
    cf.print_all();

    problem.p.p_u        = cf.get<double>("p_u");
    problem.p.alpha      = cf.get<double>("alpha");
    problem.p.beta       = cf.get<double>("beta");
    problem.p.gamma_u    = cf.get<double>("gamma_u");
    problem.p.gamma_b    = cf.get<double>("gamma_b");
    problem.p.D_su       = cf.get<double>("D_su");
    problem.p.D_ju       = cf.get<double>("D_ju");
    problem.p.nu         = cf.get<double>("nu");
    problem.p.lambda     = cf.get<double>("lambda");
    problem.p.phi_i_max  = cf.get<double>("phi_i_max");

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

    std::cout << '\n';
}
