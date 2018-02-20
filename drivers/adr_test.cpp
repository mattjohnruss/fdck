#include <advection_diffusion_reaction_problem.h>

#include <fstream>

using namespace mjrfd;

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(2, n_node, dt)
    {
        enable_bc(Boundary::Left, { 0, 1 });
        enable_bc(Boundary::Right, { 0, 1 });
        //enable_bc(Boundary::Left, { 0 });
        //enable_bc(Boundary::Right, { 0 });
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
            a1[1] = 1.0; a2[1] = 0.0; a3[1] = 0.0;
        }
        if(b == Boundary::Right)
        {
            a1[0] = 1.0; a2[0] = 0.0; a3[0] = 0.0;
            a1[1] = 1.0; a2[1] = 0.0; a3[1] = 1.0;
        }
    }

    void get_d(std::vector<double> &d) const override
    {
        d[0] = 1.0;
        d[1] = 0.1;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        if(var == 0)
        {
            //v = 10*x(i);
            //v = 6.5;

            v = 0.0;

            // v_0 now depends on u_1!
            for(const auto& [j, w] : stencil::central_1::weights)
            {
                v += 10.0*w*u(0, 1, i+j)/dx_;
            }
        }
        else if(var == 1)
        {
            v = 1.5;
        }
    }

    void get_dv_du(const unsigned i,
                   const unsigned var,
                   const unsigned i2,
                   const unsigned var2,
                   double &dv_du) const override
    {
        dv_du = 0.0;

        // only contribution is from the case d(u_0)/d(u_1)
        if(var == 0 && var2 == 1)
        {
            // loop over the velocity stencil points
            for(const auto& [j, w] : stencil::central_1::weights)
            {
                if(i+j == i2)
                {
                    dv_du += 10.0*w/dx_;
                }
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[0] = -1.0*u[1];
        r[1] = 1.0*u[0];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[0][0] = 0.0;
        dr_du[0][1] = -1.0;

        dr_du[1][0] = 1.0;
        dr_du[1][1] = 0.0;
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

    //problem.enable_fd_jacobian();

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
