#include <advection_diffusion_reaction_problem.h>

#include <fstream>

using namespace mjrfd;

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(2, n_node, dt)
    {
    }

    ~TestProblem()
    {
    }

private:
    void get_d(std::vector<double> &d) const override
    {
        d[0] = 1.0;
        d[1] = 0.1;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        //static const double pi = std::acos(-1.0);
        //v = 10*std::sin(10.0*pi*x(i));
        if(var == 0)
        {
            v = 10*x(i);
        }
        else if(var == 1)
        {
            v = 1.5;
        }
    }

    void get_dv_du(const unsigned i,
                   const std::vector<double> &u,
                   std::vector<std::vector<double>> &dv_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dv_du[var2][var2] = 0.0;
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[0] = 3.0*u[1];
        r[1] = 1.0*u[0];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dr_du[var][var2] = 0.0;
            }
        }
    }
};

int main(int argc, char **argv)
{
    unsigned n_node = 11;

    if(argc == 2)
    {
        n_node = std::atoi(argv[1]);
    }

    TestProblem problem(n_node, 0.005);

    problem.enable_fd_jacobian();
    //problem.enable_dump_jacobian();
    //problem.steady_solve();

    //unsigned n_timestep = 100;

    //problem.clear_solution();

    //for(unsigned t = 0; t < n_timestep; ++t)
    //{
        //problem.unsteady_solve();
    //}

    //std::ofstream outfile("output.dat");
    //problem.output(outfile);
    //outfile.close();

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve and output it
    problem.steady_solve();
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    // set initial conditions (zero)
    problem.clear_solution();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    double t_max = 1.0;
    unsigned output_interval = 1;

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
}
