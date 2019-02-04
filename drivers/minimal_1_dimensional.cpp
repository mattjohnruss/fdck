#include "advection_diffusion_reaction_problem.h"
#include "config.h"
#include "utilities.h"
#include "log.h"

#include <fstream>

using namespace mjrfd;

struct ChemokinesParams
{
    double D;
    double u;
    std::string ics; // Allowed values: "zero", "polynomial", and any other is treated as a data file
    std::string bcs; // Allowed values: "polynomial", and any other is treated as a data file
};

enum Variable
{
    c_u = 0,
};

struct ChemokinesBuilder
{
    double a;
    double b;
    std::vector<double> inlet_poly_coeffs;
    std::vector<double> outlet_poly_coeffs;
    std::vector<double> ic_poly_coeffs;
    double start_time;
    double end_time;
};

namespace data
{
    // version ignoring the time points before t=120
    //ChemokinesBuilder proposal_builder {
        //// a
        //0.0,
        //// b
        //0.9100772827,
        //// inlet_poly_coeffs
        //{ 0.0920171184030373,
          //0.281788395976262,
          //0.899875950593092,
         //-3.35207139698531,
          //4.07826456814076,
         //-1.71607161439073 },
        //// outlet_poly_coeffs
        //{ 0.0252628593733173,
          //0.0702608135830629,
         //-0.434142866323765,
          //1.33417142646105,
         //-0.565008183472908,
         //-0.268762291710336 },
        //// ic_poly_coeffs
        //{ 0.0903660062452174,
         //-0.232998153044166,
          //0.382087522035561,
         //-0.301365627938283,
          //0.0901094767962003,
         //-0.00177964573302294 },
        //// start_time
        //120.0,
        //// end_time
        //2730.0
    //};

    // version with all time points, for doing simulations starting at t=0 with
    // made up zero ICs and a linear ramp up to starting concentration
    ChemokinesBuilder proposal_builder {
        // a
        0.0,
        // b
        0.9100772827,
        // inlet_poly_coeffs
        { 0.079725755309440,
          0.276550445893222,
          1.044702075477787,
         -3.649471430632135,
          4.314373224651199,
         -1.781640087070480 },
        // outlet_poly_coeffs
        { 0.024176997335253,
          0.078142098914668,
         -0.461436891550468,
          1.221789288442903,
         -0.263762150324452,
         -0.437328330024978 },
        // ic_poly_coeffs
        { 0.0903660062452174,
         -0.232998153044166,
          0.382087522035561,
         -0.301365627938283,
          0.0901094767962003,
         -0.00177964573302294 },
        // start_time
        30.0,
        // end_time
        2730.0
    };

    ChemokinesBuilder april_2018_builder {
        // a
        0.0,
        // b
        0.9135985775,
        // inlet_poly_coeffs
        { 1.717652246214328e-01,
         -3.640172690898086e+00,
          5.997093506520408e+01,
         -2.935621551243856e+02,
          7.080243268663148e+02,
         -9.015067569946547e+02,
          5.808596311440184e+02,
         -1.493448744259853e+02 },
        // outlet_poly_coeffs
        { 8.868545844710741e-04,
          2.745136282303021e-01,
         -5.270461098637592e+00,
          4.978467452548591e+01,
         -1.537692402545972e+02,
          2.250891943459496e+02,
         -1.597531135191434e+02,
          4.429141978683641e+01 },
        // ic_poly_coeffs
        { 0.165334038478026,
         -0.517953087933811,
          0.120090488110154,
          1.394193954411455,
         -1.916598633553111,
          0.756766820009559 },
        // start_time
        30.0,
        // end_time
        2340.0
    };

    ChemokinesBuilder october_2018_builder {
        // a
        0.0,
        // b
        0.883700440500000,
        // inlet_poly_coeffs
        { 2.089469978370483e+00,
         -6.288794624170881e+00,
          3.895363132976494e+01,
         -3.561922955612051e+02,
          1.294768779070999e+03,
         -2.080383892606891e+03,
          1.544539202175533e+03,
         -4.353784103844123e+02 },
        // outlet_poly_coeffs
        { 2.410682232453739e-01,
         -6.776078238804784e-02,
         -3.018376713796205e-01,
         -2.185644182382399e-01,
          6.552097236925701e+00,
         -1.502540557005286e+01,
          1.338255076496678e+01,
         -4.207229174888799e+00 },
        // ic_poly_coeffs
        { 2.085761209893685e+00,
         -9.180585776152208e-01,
         -3.158942817647431e+00,
         -1.555222595486961e+01,
          4.035122488579531e+01,
         -1.623152322157534e+01,
         -1.819122005046435e+01,
          1.186140204264494e+01 },
        // start time
        210.0,
        // end_time
        600.0
    };
}

class ChemokinesProblem1D : public AdvectionDiffusionReactionProblem
{
public:
    ChemokinesProblem1D(const unsigned n_node,
                        const ChemokinesBuilder &&builder) :
        AdvectionDiffusionReactionProblem(1, n_node, 0, builder.a, builder.b),
        inlet_poly_coeffs_{builder.inlet_poly_coeffs},
        outlet_poly_coeffs_{builder.outlet_poly_coeffs},
        ic_poly_coeffs_{builder.ic_poly_coeffs},
        start_time{builder.start_time},
        end_time{builder.end_time}
    {
        enable_bc(Boundary::Left,  { c_u });
        enable_bc(Boundary::Right, { c_u });

        enable_spatial_terms({ c_u });

        disable_output_time_column();

        set_variable_names({ "c_u" });

        Max_residual = 1.0e-14;

        trace_file_.open("trace.dat");
        trace_header_ =
            "t c_u(0) c_u(1)";
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

        if(p.ics == "zero")
        {
            // Initial time is zero
            time() = 0.0;

            // solution has been cleared above and everything is zero at this point

            //// set the initial conditions to zero
            //for(unsigned i = 0; i < n_node_; ++i)
            //{
                //u(0, c_u, i) = 0.0;
            //}
        }
        else if(p.ics == "polynomial")
        {
            // Initial time is the "start time" of the experiment
            time() = start_time;

            // set the initial conditions from the polynomial fit
            for(unsigned i = 0; i < n_node_; ++i)
            {
                double scaled_x = map_space(x(i));
                u(0, c_u, i) = utilities::evaluate_polynomial(scaled_x, ic_poly_coeffs_);
            }
        }
        // Assume p.ics is the name of the ics CSV file
        else
        {
            std::ifstream ics_file(p.ics);

            if(ics_file)
            {
                auto [m_vec, n_rows, n_cols] =
                    utilities::read_csv_to_flat_vector(ics_file, ' ');

                typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> S;

                Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, S>
                    m(m_vec.data(), n_rows, n_cols, S(1, n_cols));

                Eigen::VectorXd x_interp = Eigen::VectorXd::LinSpaced(n_node_, a_, b_);
                Eigen::VectorXd c_u_interp(n_node_);

                utilities::lerp_mesh(m.col(0), m.col(1), x_interp, c_u_interp);

                for(unsigned i = 0; i < n_node_; ++i)
                {
                    u(c_u, i) = c_u_interp(i);
                }
            }
            else
            {
                MJRFD_FATAL("Error opening ICs file \"{}\"! Exiting", p.ics);
                std::exit(1);
            }
        }
    }

    ChemokinesParams p;

    void trace()
    {
        trace_file_ << time() << " "
                    << u(c_u, 0) << " "
                    << u(c_u, n_node_-1) << '\n'
                    << std::flush;
    }

private:
    std::string trace_header_;
    std::ofstream trace_file_;

    // Coeffs for the polynomial fits of the inlet/outlet BC functions and ICs
    // The polynomials are fitted to the data after scaling the "x" axis to [0,1]
    const std::vector<double> inlet_poly_coeffs_;
    const std::vector<double> outlet_poly_coeffs_;
    const std::vector<double> ic_poly_coeffs_;

public:
    double start_time;
    double end_time;

private:
    // map the spatial variable from its actual range to the interval [0,1]
    double map_space(double x) const
    {
        return (x - a_)/(b_ - a_);
    }

    // map the time variable from its actual range to the interval [0,1]
    double map_time(double time) const
    {
        return (time - start_time)/(end_time - start_time);
    }

    void actions_after_timestep() override
    {
        trace();
    }

    void get_bc(Boundary b,
                std::vector<double> &a1,
                std::vector<double> &a2,
                std::vector<double> &a3) const override
    {
        a1[c_u] = 1.0;
        a2[c_u] = 0.0;

        if(b == Boundary::Left)
        {
            // if the flag is set and we haven't reached start_time, do a linear interpolation between (0,0) and (start_time, c(start_time))
            if(p.ics == "zero" && time() < start_time)
            {
                // get the (interpolated) concentration at t = start_time
                const static double start_c =
                    utilities::evaluate_polynomial(map_time(start_time), inlet_poly_coeffs_);

                // calculate the slope of the line
                const static double m = start_c/start_time;

                a3[c_u] = m*time();
            }
            else
            {
                a3[c_u] = utilities::evaluate_polynomial(map_time(time()), inlet_poly_coeffs_);
            }
        }
        if(b == Boundary::Right)
        {
            // if the flag is set and we haven't reached start_time, do a linear interpolation between (0,0) and (start_time, c(start_time))
            if(p.ics == "zero" && time() < start_time)
            {
                // get the (interpolated) concentration at t = start_time
                const static double start_c =
                    utilities::evaluate_polynomial(map_time(start_time), outlet_poly_coeffs_);

                // calculate the slope of the line
                const static double m = start_c/start_time;

                a3[c_u] = m*time();
            }
            else
            {
                a3[c_u] = utilities::evaluate_polynomial(map_time(time()), outlet_poly_coeffs_);
            }
        }
    }

    void get_dbc_du(Boundary,
                    const unsigned,
                    std::vector<std::vector<double>> &da1_du,
                    std::vector<std::vector<double>> &da2_du,
                    std::vector<std::vector<double>> &da3_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                da1_du[var][var2] = 0.0;
                da2_du[var][var2] = 0.0;
                da3_du[var][var2] = 0.0;
            }
        }
    }

    void get_d(const unsigned,
               const unsigned,
               std::vector<double> &d) const override
    {
        d[c_u] = p.D;
    }

    void get_dd_du(const unsigned,
                   const unsigned,
                   std::vector<std::vector<double>> &dd_du) const override
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dd_du[var][var2] = 0.0;
            }
        }
    }

    void get_v(const unsigned,
               const unsigned,
               double &v) const override
    {
        v = p.u;
    }

    void get_dv_du(const unsigned,
                   const unsigned,
                   const unsigned,
                   const unsigned,
                   double &dv_du) const override
    {
        dv_du = 0.0;
    }

    void get_r(const std::vector<double> &,
               std::vector<double> &r) const override
    {
        r[c_u] = 0.0;
    }

    void get_dr_du(const std::vector<double> &,
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

    log::disable_labels();

    std::ifstream config_file(argv[1]);
    Config cf;
    cf.parse_config_file(config_file);
    cf.parse_command_line(argc, argv);
    cf.print_all();

    log::set_level(cf.get_or<std::string>("log_level", "info"));

    ChemokinesBuilder builder;

    std::string dataset = cf.get_or<std::string>("dataset", "");

    if(dataset == "proposal")
    {
        builder = std::move(data::proposal_builder);
    }
    else if(dataset == "april2018")
    {
        builder = std::move(data::april_2018_builder);
    }
    else if(dataset == "october2018")
    {
        builder = std::move(data::october_2018_builder);
    }
    else
    {
        MJRFD_FATAL("Unrecognised dataset `{}` specified! Exiting", dataset);
        std::exit(1);
    }

    ChemokinesProblem1D problem(n_node, std::move(builder));

    problem.p.u = cf.get<double>("u");
    problem.p.D = cf.get<double>("D");
    problem.p.ics = cf.get_or<std::string>("ics", "polynomial");

    problem.enable_terse_logging();

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

#ifdef _WIN32
    std::cout.sync_with_stdio(false);
#endif

    // timestepping loop
    while(problem.time() <= t_max && problem.time() <= problem.end_time)
    {
        // solve for current timestep
        problem.unsteady_solve(dt);

        if(i % output_interval == 0)
        {
            // output current solution
            MJRFD_TRACE("Outputting");
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }

    MJRFD_INFO("Reached t > min(t_max = {}, end_time = {}) after performing {} timesteps", t_max, problem.end_time, i-1);
}
