#include "advection_diffusion_reaction_problem.h"
#include "config.h"

#include <fstream>
#include <iomanip>

using namespace fdck;

struct Params
{
    double v0;
    double v1;
};

enum Variable
{
    c_0 = 0,
    c_1 = 1,
};

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(2, n_node, dt)
    {
        enable_bc(Boundary::Left,  { c_0, c_1 });
        enable_bc(Boundary::Right, { c_0, c_1 });

        enable_spatial_terms({ c_0, c_1 });

        set_variable_names({ "c_0 num", "c_1 num" });

        Max_residual = 1.0e-14;
    }

    ~TestProblem()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Set c_1 to 1.0 at the left-hand boundary
        u(0, c_1, 0) = 1.0;
    }

    Params p;

private:
    const std::unordered_map<int, double>& stencil_1_helper(const unsigned i) const
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
            double c_1_deriv = 0.0;

            for(const auto& [j, w] : stencil::forward_1::weights)
            {
                c_1_deriv += w*u(0, c_1, 0+j)/dx_;
            }

            //std::cout << c_1_deriv << '\n';

            //a1[c_0] = p.v0;             a2[c_0] = -1.0; a3[c_0] = 1.0 + c_1_deriv;
            //a1[c_1] = 2.0*u(0, c_0, 0); a2[c_1] = 0.0;  a3[c_1] = 1.0;

            a1[c_0] = c_1_deriv; a2[c_0] = 0.0; a3[c_0] = 1.0;
            a1[c_1] = 1.0;       a2[c_1] = 0.0; a3[c_1] = 1.0;
        }
        if(b == Boundary::Right)
        {
            a1[c_0] = 1.0; a2[c_0] = 0.0; a3[c_0] = 0.0;
            a1[c_1] = 1.0; a2[c_1] = 0.0; a3[c_1] = 0.0;
        }
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

        if(b == Boundary::Left)
        {
            for(const auto& [j, w] : stencil::forward_1::weights)
            {
                if(i2 == 0+j)
                {
                    da1_du[c_0][c_1] = w/dx_;
                }
            }
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        const double x = this->x(i);
        static constexpr double pi = std::acos(-1.0);

        d[c_0] = 1.0;
        d[c_1] = 1.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        if(var == c_0)
            v = p.v0;
        if(var == c_1)
            v = p.v1;
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c_0] = 0.0;
        r[c_1] = 0.0;
    }

    void exact_solution(const double time,
                        const double x,
                        std::vector<double> &sol) const override
    {
        using std::exp;

        //sol[0] = (exp(p.v1) - 1.0)*(exp(p.v0) - exp(p.v0*x))/(exp(p.v0)*((exp(p.v1) - 1.0)*p.v0 + 2.0*p.v1) - 2.0*p.v1);
        //sol[1] = 2.0*(exp(p.v0) - 1.0)*(exp(p.v1) - exp(p.v1*x))/(exp(p.v0)*((exp(p.v1) - 1.0)*p.v0 + 2.0*p.v1) - 2.0*p.v1);

        //sol[c_0] = 2.0*(exp(p.v0) - exp(p.v0*x))/(exp(p.v0) - 1.0);
        //sol[c_1] = (exp(p.v1) - exp(p.v1*x))/(exp(p.v1) - 1.0);

        sol[c_0] = (exp(p.v1) - 1.0)*(exp(p.v0*x) - exp(p.v0))/(p.v1*(exp(p.v0) - 1.0));
        sol[c_1] = (exp(p.v1) - exp(p.v1*x))/(exp(p.v1) - 1.0);
    }
};

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " config_file n_node\n";
        std::exit(1);
    }

    unsigned n_node = std::atoi(argv[2]);

    TestProblem problem(n_node, 0.01);

    std::ifstream config_file(argv[1]);
    Config cf(config_file);
    cf.print_all();

    problem.p.v0 = cf.get<double>("v0");
    problem.p.v1 = cf.get<double>("v1");

    problem.enable_terse_logging();
    problem.disable_exit_on_solve_fail();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve using exact jacobian and output it
    std::cout << "Solve using exact jacobian:\n";
    problem.disable_fd_jacobian();
    //problem.enable_dump_jacobian("ex_");

    problem.set_initial_conditions();

    problem.steady_solve();
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    outfile << std::setprecision(16);
    problem.output(outfile);
    outfile.close();

    std::cout << std::endl;

    //// perform a steady solve using fd jacobian and output it
    //std::cout << "Solve using fd jacobian:\n";
    //problem.enable_fd_jacobian();
    ////problem.enable_dump_jacobian("fd_");

    //problem.set_initial_conditions();

    //problem.steady_solve();
    //std::sprintf(filename, "output_steady_fd.csv");
    //outfile.open(filename);
    //outfile << std::setprecision(16);
    //problem.output(outfile);
    //outfile.close();

    // output the exact solution
    std::sprintf(filename, "output_exact.csv");
    outfile.open(filename);
    outfile << std::setprecision(16);
    problem.output_exact(outfile);
    outfile.close();

    std::cout << '\n';

    //// set initial conditions (zero)
    //problem.clear_solution();

    //// output initial conditions
    //std::sprintf(filename, "output_%05u.csv", 0);
    //outfile.open(filename);
    //problem.output(outfile);
    //outfile.close();

    //unsigned i = 1;

    //double t_max = 10.0;
    //unsigned output_interval = 1;

    //// timestepping loop
    //while(problem.time() <= t_max)
    //{
        //// solve for current timestep
        //problem.unsteady_solve();

        //if(i % output_interval == 0)
        //{
            //// output current solution
            ////std::cout << "Outputting solution at time = " << problem.time() << '\n';
            //std::cout << ";\tOutputting";
            //std::sprintf(filename, "output_%05u.csv", i/output_interval);
            //outfile.open(filename);
            //problem.output(outfile);
            //outfile.close();
        //}

        //++i;
    //}

    //std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
}
