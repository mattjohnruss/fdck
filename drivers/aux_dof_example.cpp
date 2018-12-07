#include <advection_diffusion_reaction_problem.h>

#include <iostream>
#include <fstream>

enum Variable
{
    c_0 = 0,
    c_1 = 1,
};

enum AuxVariable
{
    y = 0,
};

using namespace mjrfd;

class TransportWithAuxDofProblem : public AdvectionDiffusionReactionProblem
{
public:
    TransportWithAuxDofProblem(const unsigned n_node) :
        AdvectionDiffusionReactionProblem(2, n_node, 1),
        trace_header_{"t c_0(0) c_0(1) c_1(0) c_1(1) y(t)"},
        trace_file_{"trace.dat"}
        //, trace_count_(0)
    {
        enable_bc(Boundary::Left,  { c_0, c_1 });
        enable_bc(Boundary::Right, { c_0, c_1 });

        enable_spatial_terms({ c_0, c_1 });

        set_variable_names({ "c_0", "c_1" });

        Max_residual = 1.0e-14;
    }

    ~TransportWithAuxDofProblem()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        //for(unsigned i = 0; i < n_node_; ++i)
        //{
            //u(0, c_0, i) = 0;
        //}

        u_aux(y) = 1.0;
    }

    void trace()
    {
        trace_file_ << time() << " "
                    << u(c_0, 0) << " "
                    << u(c_0, n_node_-1) << " "
                    << u(c_1, 0) << " "
                    << u(c_1, n_node_-1) << " "
                    << u_aux(y) << '\n';
        trace_file_.flush();
    }

    //void trace(Eigen::VectorXd &solution, Eigen::VectorXd &exact_solution)
    //{
        //solution(trace_count_) = u_aux(y);
        //exact_solution(trace_count_) = this->exact_solution(time());

        //trace_count_++;
    //}

private:
    //unsigned trace_count_;
    std::string trace_header_;
    std::ofstream trace_file_;

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            a1[c_0] = 1.0; a2[c_0] = 0.0; a3[c_0] = u_aux(y);
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
        d[c_0] = 1.0;
        d[c_1] = 1.0;
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 1.0;
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c_0] = 0.0;
        r[c_1] = 0.0;
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        // Add the residuals for the transport pdes
        AdvectionDiffusionReactionProblem::calculate_residual(residual);

        // Add the residual for the ode
        residual(aux_dof_index(y)) += (u_aux(0, y) - u_aux(1, y)) - dt_*u_aux(y)*u(c_1, n_node_/2);
    }

    void calculate_jacobian(std::vector<Triplet> &triplet_list) const override
    {
        // Add the jacobian entries for the transport pdes
        AdvectionDiffusionReactionProblem::calculate_jacobian(triplet_list);

        // Add the jacobian entry for the ode
        triplet_list.emplace_back(aux_dof_index(y),
                                  aux_dof_index(y),
                                  1.0 - dt_*u(c_1, n_node_/2));

        // Add the jacobian entries for the coupling
        // We have to include the dx^2 factor and -1 sign since this is how
        // AdvectionDiffusionReactionProblem constructs the BC residuals - this
        // is low level and a bit messy but works fine. We could maybe create
        // some ODE helper functions to hide these details later...
        triplet_list.emplace_back(nodal_dof_index(c_0, 0),
                                  aux_dof_index(y),
                                  -1.0*dx_*dx_);

        triplet_list.emplace_back(aux_dof_index(y),
                                  nodal_dof_index(c_1, n_node_/2),
                                  -dt_*u_aux(y));
    }
};

int main(int argc, char **argv)
{
    if(argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[1]);
    const double dt = std::atof(argv[2]);
    const double t_max = std::atof(argv[3]);

    unsigned output_interval = 1;

    if(argc == 5)
    {
        output_interval = std::atoi(argv[4]);
    }

    TransportWithAuxDofProblem problem(n_node);

    problem.enable_terse_logging();
    //problem.enable_fd_jacobian();

    char filename[200];
    std::ofstream outfile;

    // perform a steady solve and output it
    //problem.steady_solve();
    //std::sprintf(filename, "output_steady.csv");
    //outfile.open(filename);
    //problem.output(outfile);
    //outfile.close();

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    problem.trace();

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
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();

            problem.trace();
        }

        ++i;
    }

    std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
}
