#include <advection_diffusion_reaction_problem.h>
#include <config.h>
#include <functions.h>

#include <fstream>
#include <memory>

using namespace mjrfd;

struct ChemokinesParams
{
    double pe_u;
    double alpha;
    double beta;
    double gamma_ui;
    double gamma_bi;
    double gamma_um;
    double gamma_bm;
    double q_u;
    double q_b;
    double q_s;
    double D_su;
    double D_iu;
    std::unique_ptr<DifferentiableFunction> D_mu;
    double nu_u;
    double nu_b;
    double nu_s;
    double phi_i_init;
    double R;
    double M;
    double M_a;
    double M_b;
    double t1;
    double t2;
    double t3;
    double t4;
    double J_m_left_prop;
    double J_m_right_prop;
    double J_m_left_abs;
    double J_m_right_abs;
    double J_i_right;
    std::unique_ptr<DifferentiableFunction> chi;
};

enum Variable
{
    c_u   = 0,
    c_b   = 1,
    c_s   = 2,
    phi_i = 3,
    phi_m = 4
};

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        AdvectionDiffusionReactionProblem(5, n_node, dt)
    {
        enable_bc(Boundary::Left,  { c_u, c_s, phi_i, phi_m });
        enable_bc(Boundary::Right, { c_u, c_s, phi_i, phi_m });

        enable_spatial_terms({ c_u, c_s, phi_i, phi_m });

        set_variable_names({ "c_u", "c_b", "c_s", "phi_i", "phi_m" });

        Max_residual = 1.0e-14;
        Max_newton_iterations = 50;

        trace_file_.open("trace.dat");
        trace_header_ =
            "t phi_m|_0 dphi_{m}\\\\_dx|_0 M total\\\\_c_u total\\\\_c_b total\\\\_c_s total\\\\_phi_i total\\\\_phi_m";
        trace_file_ << trace_header_ << '\n';
    }

    ~ChemokinesProblem1D()
    {
        trace_file_.close();
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        // Set c_u to 1.0 at the left-hand boundary
        //u(0, c_u, 0) = 1.0;

        // Set phi_i to its uniform initial density
        for(unsigned i = 0; i < n_node_; ++i)
        {
            u(0, phi_i, i) = p.phi_i_init;
        }

        for(unsigned i = 0; i < n_node_; ++i)
        {
            //u(0, c_u, i) = 0.5 + 0.5*std::cos(M_PI*x(i));
            u(0, c_u, i) = 1.0 - x(i);
        }
    }

    ChemokinesParams p;

    static double cos_ramp_down(const double time)
    {
        assert(time >= 0.0 && time <= 1.0);

        static const double pi = std::acos(-1.0);
        return 0.5*(1.0 + std::cos(pi*time));
    }

    static double cos_ramp_up(const double time)
    {
        assert(time >= 0.0 && time <= 1.0);

        static const double pi = std::acos(-1.0);
        return 1.0 - 0.5*(1.0 + std::cos(pi*time));
    }

private:
    std::string trace_header_;
    std::ofstream trace_file_;

    void actions_before_timestep() override
    {
        // update the maturation parameter before solving for the new time
        p.M = maturation_piecewise_ramp(p.M_a, p.M_b, p.t1, p.t2, p.t3, p.t4);
    }

    void actions_after_timestep() override
    {
        std::vector<double> a1(n_var_);
        std::vector<double> a2(n_var_);
        std::vector<double> a3(n_var_);

        get_bc(Boundary::Left, a1, a2, a3);

        double dphi_m_dx = 0.0;
        for(auto& [j, w] : stencil::forward_1::weights)
        {
            dphi_m_dx += w*u(0, phi_m, 0+j)/dx_;
        }

        if(is_steady() == true)
        {
            // if we're doing a steady solve after timestepping, output the
            // trace data to a new data block to keep gnuplot from plotting it
            // on the unsteady data
            trace_file_ << "\n\n" << trace_header_ << '\n';
        }

        trace_file_ << time() << " "
                    << u(0, phi_m, 0) << " "
                    << dphi_m_dx << " "
                    << p.M << " "
                    << integrate_solution(c_u) << " "
                    << integrate_solution(c_b) << " "
                    << integrate_solution(c_s) << " "
                    << integrate_solution(phi_i) << " "
                    << integrate_solution(phi_m) << '\n'
                    << std::flush;
    }

    double maturation_piecewise_ramp(const double a,  const double b,
                                     const double t1, const double t2,
                                     const double t3, const double t4) const
    {
        assert(t1 <= t2);
        assert(t2 <= t3);
        assert(t3 <= t4);

        double m = 0.0;
        double time = this->time();

        if(time < t1 || time > t4)
            m = a;
        else if(time >= t1 && time <= t2)
            m = a + (b - a)*cos_ramp_up((time - t1)/(t2 - t1));
        else if(time > t2 && time < t3)
            m = b;
        else if(time >= t3 && time <= t4)
            m = a + (b - a)*cos_ramp_down((time - t3)/(t4 - t3));

        return m;
    }

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        if(b == Boundary::Left)
        {
            // get the phi_m advection velocity at the left boundary
            double v_phi_m_left = 0;
            get_v(0, phi_m, v_phi_m_left);

            // get the phi_m diffusivity at the left boundary
            // this gets all components of the diffusivity but oh well
            std::vector<double> d_left(n_var_);
            get_d(0, 0, d_left);

            a1[c_u]   = 1.0;                            a2[c_u]   = 0.0;            a3[c_u]   = 1.0;
            a1[c_s]   = 1.0;                            a2[c_s]   = 0.0;            a3[c_s]   = 0.0;
            a1[phi_i] = 0.0;                            a2[phi_i] = 1.0;            a3[phi_i] = 0.0;
            a1[phi_m] = v_phi_m_left + p.J_m_left_prop; a2[phi_m] = -d_left[phi_m]; a3[phi_m] = -p.J_m_left_abs;
        }
        if(b == Boundary::Right)
        {
            // get the phi_m advection velocity at the left boundary
            double v_phi_m_right = 0;
            get_v(n_node_-1, phi_m, v_phi_m_right);

            // get the phi_m diffusivity at the right boundary
            // this gets all components of the diffusivity but oh well
            std::vector<double> d_right(n_var_);
            get_d(0, n_node_-1, d_right);

            a1[c_u]   = 1.0;                              a2[c_u]   = 0.0;             a3[c_u]   = 0.0;
            a1[c_s]   = 1.0;                              a2[c_s]   = 0.0;             a3[c_s]   = 0.0;
            a1[phi_i] = 0.0;                              a2[phi_i] = 1.0;             a3[phi_i] = -p.J_i_right*p.M;
            a1[phi_m] = v_phi_m_right - p.J_m_right_prop; a2[phi_m] = -d_right[phi_m]; a3[phi_m] = p.J_m_right_abs;
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

        unsigned i = 0;

        if(b == Boundary::Left)
            i = 0;
        else if(b == Boundary::Right)
            i = n_node_-1;

        std::vector<std::vector<double>> dd_du(n_var_, std::vector<double>(n_var_));
        get_dd_du(0, i, dd_du);

        // loop over the variables we know that a1[phi_m], a2[phi_m] depend on
        for(auto var2 : { c_u, c_s, c_b })
        {
            // terms from v_phi_m
            double dv_du_var2 = 0.0;
            get_dv_du(i, phi_m, i2, var2, dv_du_var2);
            da1_du[phi_m][var2] += dv_du_var2;

            // terms from D_mu
            da2_du[phi_m][var2] += -dd_du[phi_m][var2];
        }
    }

    void get_d(const unsigned t,
               const unsigned i,
               std::vector<double> &d) const override
    {
        d[c_u]   = 1.0;
        d[c_b]   = 0.0;
        d[c_s]   = p.D_su;
        d[phi_i] = p.D_iu;
        d[phi_m] = p.D_mu->value(u(t, c_b, i));
    }

    void get_dd_du(const unsigned t,
                   const unsigned i,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        for(auto var : { c_u, c_b, c_s, phi_i, phi_m })
        {
            std::fill(dd_du[var].begin(), dd_du[var].end(), 0.0);
        }

        dd_du[phi_m][phi_m] = p.D_mu->deriv(u(t, c_b, i));
    }

    void get_v(const unsigned i,
               const unsigned var,
               double &v) const override
    {
        v = 0.0;

        if(var == c_u)
            v = p.pe_u;
        else if(var == c_s)
            v = p.pe_u;
        else if(var == phi_m)
        {
            //double dc_u_dx = 0.0;
            double dc_b_dx = 0.0;
            //double dc_s_dx = 0.0;

            for(const auto& [j, w] : central_1_stencil_weights(i))
            {
                //dc_u_dx += w*u(0, c_u, i+j)/dx_;
                dc_b_dx += w*u(0, c_b, i+j)/dx_;
                //dc_s_dx += w*u(0, c_s, i+j)/dx_;
            }

            //v += p.nu_u*p.chi->value(u(0, c_u, i))*dc_u_dx;
            v += p.nu_b*p.chi->value(u(0, c_b, i))*dc_b_dx;
            //v += p.nu_s*p.chi->value(u(0, c_s, i))*dc_s_dx;
        }
    }

    void get_dv_du(const unsigned i,
                   const unsigned var,
                   const unsigned i2,
                   const unsigned var2,
                   double &dv_du) const override
    {
        dv_du = 0.0;

        //if(var == phi_m && (var2 == c_u || var2 == c_b || var2 == c_s))
        if(var == phi_m && var2 == c_b)
        {
            double nu;

            if(var2 == c_u)
                nu = p.nu_u;
            else if(var2 == c_b)
                nu = p.nu_b;
            else if(var2 == c_s)
                nu = p.nu_s;
            else
                std::cerr << "var2 is not any of the expected values!\n";

            double dvar2_dx = 0.0;

            for(const auto& [j, w] : central_1_stencil_weights(i))
            {
                if(i+j == i2)
                {
                    // this condition will be true at most once per loop so we
                    // could break here, but we need to continue the loop to
                    // calculate the derivative below
                    dv_du += w*nu*p.chi->value(u(0, var2, i))/dx_;
                }

                dvar2_dx += w*u(0, var2, i+j)/dx_;
            }

            dv_du += nu*p.chi->deriv(u(0, var2, i))*dvar2_dx;
        }
    }

    void get_r(const std::vector<double> &u,
               std::vector<double> &r) const override
    {
        r[c_u]   = -p.alpha*u[c_u] + p.beta*u[c_b] - p.gamma_ui*u[phi_i]*u[c_u] - p.gamma_um*u[phi_m]*u[c_u] - p.q_u*u[phi_i]*u[c_u];
        r[c_b]   = p.alpha*u[c_u] - p.beta*u[c_b] - p.gamma_bi*u[phi_i]*u[c_b] - p.gamma_bm*u[phi_m]*u[c_b] - p.q_b*u[phi_i]*u[c_b];
        r[c_s]   = p.gamma_ui*u[phi_i]*u[c_u] + p.gamma_um*u[phi_m]*u[c_u] + p.gamma_bi*u[phi_i]*u[c_b] + p.gamma_bm*u[phi_m]*u[c_b] - p.q_s*u[phi_i]*u[c_s];
        r[phi_i] = p.R*u[phi_i]*(1.0 - u[phi_i]) - p.M*u[phi_i];
        r[phi_m] = p.M*u[phi_i];
    }

    void get_dr_du(const std::vector<double> &u,
                   std::vector<std::vector<double>> &dr_du) const override
    {
        dr_du[c_u][c_u]   = -p.alpha - p.gamma_ui*u[phi_i] - p.gamma_um*u[phi_m] - p.q_u*u[phi_i];
        dr_du[c_u][c_b]   = p.beta;
        dr_du[c_u][c_s]   = 0.0;
        dr_du[c_u][phi_i] = -p.gamma_ui*u[c_u] - p.q_u*u[c_u];
        dr_du[c_u][phi_m] = -p.gamma_um*u[c_u];

        dr_du[c_b][c_u]   = p.alpha;
        dr_du[c_b][c_b]   = -p.beta - p.gamma_bi*u[phi_i] - p.gamma_bm*u[phi_m] - p.q_b*u[phi_i];
        dr_du[c_b][c_s]   = 0.0;
        dr_du[c_b][phi_i] = -p.gamma_bi*u[c_b] - p.q_b*u[c_b];
        dr_du[c_b][phi_m] = -p.gamma_bm*u[c_b];

        dr_du[c_s][c_u]   = p.gamma_ui*u[phi_i] + p.gamma_um*u[phi_m];
        dr_du[c_s][c_b]   = p.gamma_bi*u[phi_i] + p.gamma_bm*u[phi_m];
        dr_du[c_s][c_s]   = -p.q_s*u[phi_i];
        dr_du[c_s][phi_i] = p.gamma_ui*u[c_u] + p.gamma_bi*u[c_b] - p.q_s*u[c_s];
        dr_du[c_s][phi_m] = p.gamma_um*u[c_u] + p.gamma_bm*u[c_b];

        dr_du[phi_i][c_u]   = 0.0;
        dr_du[phi_i][c_b]   = 0.0;
        dr_du[phi_i][c_s]   = 0.0;
        dr_du[phi_i][phi_i] = p.R*(1.0 - 2.0*u[phi_i]) - p.M;
        dr_du[phi_i][phi_m] = 0.0;

        dr_du[phi_m][c_u]   = 0.0;
        dr_du[phi_m][c_b]   = 0.0;
        dr_du[phi_m][c_s]   = 0.0;
        dr_du[phi_m][phi_i] = p.M;
        dr_du[phi_m][phi_m] = 0.0;
    }
};

int main(int argc, char **argv)
{
    if(argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " config_file n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[2]);
    const double dt = std::atof(argv[3]);
    const double t_max = std::atof(argv[4]);

    unsigned output_interval = 1;

    if(argc >= 6)
    {
        output_interval = std::atoi(argv[5]);
    }

    ChemokinesProblem1D problem(n_node, dt);

    std::ifstream config_file(argv[1]);
    Config cf;
    cf.parse_config_file(config_file);
    cf.parse_command_line(argc, argv);
    cf.print_all();

    problem.p.pe_u           = cf.get<double>("pe_u");
    problem.p.alpha          = cf.get<double>("alpha");
    problem.p.beta           = cf.get<double>("beta");
    problem.p.gamma_ui       = cf.get<double>("gamma_ui");
    problem.p.gamma_bi       = cf.get<double>("gamma_bi");
    problem.p.gamma_um       = cf.get<double>("gamma_um");
    problem.p.gamma_bm       = cf.get<double>("gamma_bm");
    problem.p.q_u            = cf.get<double>("q_u");
    problem.p.q_b            = cf.get<double>("q_b");
    problem.p.q_s            = cf.get<double>("q_s");
    problem.p.D_su           = cf.get<double>("D_su");
    problem.p.D_iu           = cf.get<double>("D_iu");
    problem.p.nu_u           = cf.get<double>("nu_u");
    problem.p.nu_b           = cf.get<double>("nu_b");
    problem.p.nu_s           = cf.get<double>("nu_s");
    problem.p.phi_i_init     = cf.get<double>("phi_i_init");
    problem.p.R              = cf.get<double>("R");
    problem.p.M_a            = cf.get<double>("M_a");
    problem.p.M_b            = cf.get<double>("M_b");
    problem.p.M              = problem.p.M_a;
    problem.p.t1             = cf.get<double>("t1");
    problem.p.t2             = cf.get<double>("t2");
    problem.p.t3             = cf.get<double>("t3");
    problem.p.t4             = cf.get<double>("t4");
    problem.p.J_m_left_prop  = cf.get<double>("J_m_left_prop");
    problem.p.J_m_right_prop = cf.get<double>("J_m_right_prop");
    problem.p.J_m_left_abs   = cf.get<double>("J_m_left_abs");
    problem.p.J_m_right_abs  = cf.get<double>("J_m_right_abs");
    problem.p.J_i_right      = cf.get<double>("J_i_right");

    // Construct the type of DifferentiableFunction in the config for chi
    // the std::unique_ptr will destroy the object for us when it goes out of scope
    std::string chi_type = cf.get<std::string>("chi");

    if(chi_type == "constant")
    {
        // TODO remove
        std::cout << "chi is constant with value " << cf.get<double>("chi_const_val") << '\n';

        problem.p.chi =
            std::make_unique<ConstantFunction>(cf.get<double>("chi_const_val"));
    }
    else if(chi_type == "hill")
    {
        // TODO remove
        std::cout << "chi is hill with "
                  << "a = " << cf.get<double>("chi_hill_a") << ", "
                  << "n = " << cf.get<double>("chi_hill_n") << ", "
                  << "min = " << cf.get<double>("chi_hill_min") << ", "
                  << "max = " << cf.get<double>("chi_hill_max") << '\n';;

        problem.p.chi =
            std::make_unique<HillFunction>(cf.get<double>("chi_hill_a"),
                                           cf.get<double>("chi_hill_n"),
                                           cf.get<double>("chi_hill_min"),
                                           cf.get<double>("chi_hill_max"));
    }
    else
    {
        std::cerr << "unrecognised value of \"chi\" parameter:" << chi_type << '\n';
        std::exit(1);
    }

    // Do a similar thing for D_mu
    std::string D_mu_type = cf.get<std::string>("D_mu");

    if(D_mu_type == "constant")
    {
        // TODO remove
        std::cout << "D_mu is constant with value " << cf.get<double>("D_mu_const_val") << '\n';

        problem.p.D_mu = std::make_unique<ConstantFunction>(cf.get<double>("D_mu_const_val"));
    }
    else if(D_mu_type == "hill")
    {
        // TODO remove
        std::cout << "D_mu is hill with "
                  << "a = "   << cf.get<double>("D_mu_hill_a") << ", "
                  << "n = "   << cf.get<double>("D_mu_hill_n") << ", "
                  << "min = " << cf.get<double>("D_mu_hill_min") << ", "
                  << "max = " << cf.get<double>("D_mu_hill_max") << '\n';;

        problem.p.D_mu =
            std::make_unique<HillFunction>(cf.get<double>("D_mu_hill_a"),
                                           cf.get<double>("D_mu_hill_n"),
                                           cf.get<double>("D_mu_hill_min"),
                                           cf.get<double>("D_mu_hill_max"));
    }
    else
    {
        std::cerr << "unrecognised value of \"D_mu\" parameter:" << D_mu_type << '\n';
        std::exit(1);
    }

    const bool do_steady_solve = cf.get<bool>("steady");
    const bool do_time_evolution = cf.get<bool>("time_evo");

    const bool fd_jacobian = cf.get<bool>("fd_jacobian");

    problem.enable_terse_logging();

    if(fd_jacobian == true)
    {
        problem.enable_fd_jacobian();
    }

    char filename[200];
    std::ofstream outfile;

    if(do_time_evolution)
    {
        std::cout << "\nTime evolution:";

        problem.enable_exit_on_solve_fail();

        // set initial condition
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

        std::cout << "\n\nReached t > t_max (" << t_max << ") after performing " << i-1 << " timesteps\n";
    }

    if(do_steady_solve)
    {
        std::cout << "\nSteady solve:\n";

        problem.disable_exit_on_solve_fail();

        // set initial conditions again - required since phi_i is sensitive to
        // ICs even at steady state
        // Only do this if we haven't timestepped - if we have been
        // timestepping, we want to use the most recent solution as the initial
        // guess for the steady solve, or it is much less likely to converge
        if(do_time_evolution == false)
        {
            problem.set_initial_conditions();
        }

        // perform a steady solve and output it
        problem.steady_solve();
        std::sprintf(filename, "output_steady.csv");
        outfile.open(filename);
        problem.output(outfile);
        outfile.close();

        std::cout << '\n';
    }

    return 0;
}
