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

    bool dump = true;

    problem.steady_solve(dump);

    std::ofstream outfile("output.dat");
    problem.output(outfile);
    outfile.close();

    //Eigen::VectorXd x(3);
    //x(0) = 3.4;
    //x(1) = 5.1;
    //x(2) = 9.3;

    //Eigen::VectorXd x_plus(3), x_minus(3);

    //double dx = 0.1;
    //x_plus = x.array() - dx;

    //std::cout << x << '\n';
    //std::cout << x_plus << '\n';
}
