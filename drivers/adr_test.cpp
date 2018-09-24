#include <advection_diffusion_reaction_problem.h>

#include <fstream>
#include <iomanip>

#include <fenv.h>

using namespace mjrfd;

enum Variable
{
    c_0 = 0,
    c_1 = 1,
};

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(2, n_node)
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

private:
    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c_0] = 1.0; a2[c_0] = 0.0; a3[c_0] = 1.0;
            a1[c_1] = 1.0; a2[c_1] = 0.0; a3[c_1] = 1.0;
        }
        if(b == Boundary::Right)
        {
            a1[c_0] = 1.0; a2[c_0] = 0.0; a3[c_0] = 0.0;
            a1[c_1] = 1.0; a2[c_1] = 0.0; a3[c_1] = 0.0;
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        const double x = this->x(i);

        d[c_0] = 1.0;
        d[c_1] = 2.0 - 5.0*x*u(t, c_0, i);
    }

    void get_dd_du(const unsigned t,
                   const unsigned i,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        const double x = this->x(i);

        dd_du[c_0][c_0] = 0.0;
        dd_du[c_0][c_1] = 0.0;

        dd_du[c_1][c_0] = -5.0*x;
        dd_du[c_1][c_1] = 0.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == c_1)
        {
            v = 1.0;
        }
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
        r[c_0] = 0.0;
        r[c_1] = 0.0;
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_0][c_0] = 0.0;
        dr_du[c_0][c_1] = 0.0;

        dr_du[c_1][c_1] = 0.0;
        dr_du[c_1][c_1] = 0.0;
    }

    void exact_solution(const double time,
                        const double x,
                        std::vector<double> &sol) const override
    {
        using std::atan;
        using std::exp;
        using std::sqrt;

        sol[0] = 1.0 - x;
        sol[1] = (exp((4.0*atan(sqrt(5.0/3.0)))/sqrt(15.0)) - exp((2.0*(atan(sqrt(5.0/3.0)) + atan(sqrt(5.0/3.0)*(-1.0 + 2.0*x))))/sqrt(15.0)))/(-1.0 + exp((4.0*atan(sqrt(5.0/3.0)))/sqrt(15.0)));
    }
};

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    unsigned n_node = 11;

    if(argc == 2)
    {
        n_node = std::atoi(argv[1]);
    }

    //double dt = 0.01;

    TestProblem problem(n_node);

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve using exact jacobian and output it
    std::cout << "Solve using exact jacobian:\n";
    problem.disable_fd_jacobian();
    //problem.enable_dump_jacobian("ex_");

    problem.steady_solve();
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    outfile << std::setprecision(16);
    problem.output(outfile);
    outfile.close();

    std::cout << std::endl;

    std::cout << "integral of c_0 = "
              << problem.integrate_solution(c_0)
              << '\n';

    std::cout << "integral of c_1 = "
              << problem.integrate_solution(c_1)
              << '\n';

    // perform a steady solve using fd jacobian and output it
    std::cout << "Solve using fd jacobian:\n";
    problem.enable_fd_jacobian();
    //problem.enable_dump_jacobian("fd_");

    problem.clear_solution();

    problem.steady_solve();
    std::sprintf(filename, "output_steady_fd.csv");
    outfile.open(filename);
    outfile << std::setprecision(16);
    problem.output(outfile);
    outfile.close();

    std::cout << std::endl;

    std::cout << "integral of c_0 = "
              << problem.integrate_solution(c_0)
              << '\n';

    std::cout << "integral of c_1 = "
              << problem.integrate_solution(c_1)
              << '\n';

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
    //std::sprintf(filename, "output_%05i.csv", 0);
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
        //problem.unsteady_solve(dt);

        //if(i % output_interval == 0)
        //{
            //// output current solution
            ////std::cout << "Outputting solution at time = " << problem.time() << '\n';
            //std::cout << ";\tOutputting";
            //std::sprintf(filename, "output_%05i.csv", i/output_interval);
            //outfile.open(filename);
            //problem.output(outfile);
            //outfile.close();
        //}

        //++i;
    //}

    //std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
}
