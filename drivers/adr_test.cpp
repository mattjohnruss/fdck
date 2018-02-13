#include <advection_diffusion_reaction_problem.h>

#include <fstream>

using namespace mjrfd;

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(1, n_node, dt)
    {
        //enable_bc(Boundary::Left, { 0, 1 });
        //enable_bc(Boundary::Right, { 0, 1 });
        enable_bc(Boundary::Left, { 0 });
        enable_bc(Boundary::Right, { 0 });
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
            a1[0] = 1.0; a2[0] = 0.0; a3[0] = 1.0;
            //a1[1] = 0.0; a2[1] = 1.0; a3[1] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[0] = 1.0; a2[0] = 0.0; a3[0] = 0.0;
            //a1[1] = 1.0; a2[1] = 1.0; a3[1] = 1.0;
        }
    }

    void get_d(std::vector<double> &d) const override
    {
        d[0] = 1.0;
        //d[1] = 0.1;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        if(var == 0)
        {
            //v = 10*x(i);
            v = 6.5;
        }
        //else if(var == 1)
        //{
            //v = 1.5;
        //}
    }

    void get_dv_du(const unsigned i,
                   const std::vector<double> &u,
                   std::vector<std::vector<double>> &dv_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dv_du[var][var2] = 0.0;
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[0] = -1.0*u[0];
        //r[1] = 1.0*u[0];
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

    TestProblem problem(n_node, 0.01);

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

    double t_max = 10.0;
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
