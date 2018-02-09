#include <advection_diffusion_reaction_problem.h>

#include <fstream>

using namespace mjrfd;

class TestProblem : public AdvectionDiffusionReactionProblem
{
public:
    TestProblem(const unsigned n_node,
                const unsigned dt) :
        AdvectionDiffusionReactionProblem(1, n_node, dt)
    {
    }

    ~TestProblem()
    {
    }

private:
    void get_d(std::vector<double> &d) const override
    {
        for(unsigned i = 0; i < n_var_; ++i)
        {
            d[i] = 1.0;
        }
    }

    void get_v(const unsigned node,
               const std::vector<double> &u,
               std::vector<double> &v) const override
    {
        for(unsigned i = 0; i < n_var_; ++i)
        {
            v[i] = 0.0;
        }
    }

    void get_dv_du(const unsigned node,
                   const std::vector<double> &u,
                   std::vector<std::vector<double>> &dv_du) const override
    {
        for(unsigned i = 0; i < n_var_; ++i)
        {
            for(unsigned j = 0; j < n_var_; ++j)
            {
                dv_du[i][j] = 0.0;
            }
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        for(unsigned i = 0; i < n_var_; ++i)
        {
            r[i] = 0.0;
        }
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        for(unsigned i = 0; i < n_var_; ++i)
        {
            for(unsigned j = 0; j < n_var_; ++j)
            {
                dr_du[i][j] = 0.0;
            }
        }
    }
};

int main(int argc, char **argv)
{
    TestProblem problem(11, 0.01);

    //bool dump = true;
    bool dump = false;

    std::cout << "solving with the exact jacobian:\n";
    problem.steady_solve(dump);

    std::cout << "solving with the finite differenced jacobian:\n";
    problem.enable_fd_jacobian();
    problem.clear_solution();
    problem.steady_solve(dump);

    std::ofstream outfile("output.dat");
    problem.output(outfile);
    outfile.close();
}
