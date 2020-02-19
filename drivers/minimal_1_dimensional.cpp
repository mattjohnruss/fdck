#include "advection_diffusion_reaction_problem.h"
#include "config.h"
#include "utilities.h"
#include "log.h"

#include <fstream>

using namespace fdck;

#define FDCK_TRACE_FILE
#define FDCK_DIST_MODE

struct ChemokinesParams
{
    double D;
    double u;
    std::string ics; // Allowed values: "zero", "polynomial", and any other is treated as a data file
    //std::string bcs; // Allowed values: "polynomial", and any other is treated as a data file
    char csv_delimiter;
    unsigned csv_skip_rows;
    double start_time;
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
    double fit_start_time;
    double fit_end_time;
};

namespace data
{
    // version ignoring the time points before t=120
    ChemokinesBuilder proposal_builder {
        // a
        0.0,
        // b
        0.9100772827,
        // inlet_poly_coeffs
        { 0.0920171184030373,
          0.281788395976262,
          0.899875950593092,
         -3.35207139698531,
          4.07826456814076,
         -1.71607161439073 },
        // outlet_poly_coeffs
        { 0.0252628593733173,
          0.0702608135830629,
         -0.434142866323765,
          1.33417142646105,
         -0.565008183472908,
         -0.268762291710336 },
        // ic_poly_coeffs
        { 0.0903660062452174,
         -0.232998153044166,
          0.382087522035561,
         -0.301365627938283,
          0.0901094767962003,
         -0.00177964573302294 },
        // fit_start_time
        120.0,
        // fit_end_time
        2730.0
    };

    // version with all time points, for doing simulations starting at t=0 with
    // made up zero ICs and a linear ramp up to starting concentration
    //ChemokinesBuilder proposal_zero_builder {
        //// a
        //0.0,
        //// b
        //0.9100772827,
        //// inlet_poly_coeffs
        //{ 0.079725755309440,
          //0.276550445893222,
          //1.044702075477787,
         //-3.649471430632135,
          //4.314373224651199,
         //-1.781640087070480 },
        //// outlet_poly_coeffs
        //{ 0.024176997335253,
          //0.078142098914668,
         //-0.461436891550468,
          //1.221789288442903,
         //-0.263762150324452,
         //-0.437328330024978 },
        //// ic_poly_coeffs
        //{ 0.0903660062452174,
         //-0.232998153044166,
          //0.382087522035561,
         //-0.301365627938283,
          //0.0901094767962003,
         //-0.00177964573302294 },
        //// fit_start_time
        //30.0,
        //// fit_end_time
        //2730.0
    //};

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
        // fit_start_time
        30.0,
        // fit_end_time
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
        // fit_end_time
        600.0
    };

    ChemokinesBuilder june_2019_builder {
        // a
        0.0,
        // b
        0.4966,
        // inlet_poly_coeffs
        { 1.274860854700616e-01,
         -2.518604210006898e-01,
          7.710766377203647e+00,
         -4.252175480340080e+01,
          1.418179464387050e+02,
         -2.375102919110578e+02,
          1.885378484147790e+02,
         -5.691021409522835e+01 },
        // outlet_poly_coeffs
        { 2.633576534354847e-06,
         -2.528577356545149e-03,
          3.101007808239337e-01,
         -1.641112220170440e+00,
          4.177919401936037e+00,
         -5.256685249420444e+00,
          3.191692523139487e+00,
         -7.367432289520273e-01 },
        // ic_poly_coeffs
        { 1.254966070500784e-01,
         -5.113492603989770e-01,
          1.557271081373574e+00,
         -2.833601879386733e+00,
          2.497843589400179e+00,
         -8.359940722753865e-01 },
        // start time
        0.0,
        // fit_end_time
        240.0
    };

    ChemokinesBuilder october_2019_builder {
        // a
        0.0,
        // b
        0.8284,
        // inlet_poly_coeffs
        { 8.659898215947195e-02,
          5.816072353204450e-01,
         -1.839989360108578e+00,
          2.122982326848622e+01,
         -7.466692445303730e+01,
          1.235604949819846e+02,
         -9.750757257841730e+01,
          2.953834540361141e+01 },
        // outlet_poly_coeffs
        { 8.077267859772258e-03,
          1.739979517326950e-01,
          2.010604972517853e-01,
         -4.417652941108896e+00,
          1.999608499820827e+01,
         -3.677705257034391e+01,
          3.093370737262772e+01,
         -9.824767532415690e+00 },
        // ic_poly_coeffs
        { 8.985858360360878e-02,
         -1.653276524259068e-01,
          2.441368877205544e-03,
          1.993037678999328e-01,
         -1.622117529248056e-01,
          4.959499333720381e-02 },
        // start time
        0.0,
        // fit_end_time
        3540.0
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
        fit_start_time{builder.fit_start_time},
        fit_end_time{builder.fit_end_time}
    {
        enable_bc(Boundary::Left,  { c_u });
        enable_bc(Boundary::Right, { c_u });

        enable_spatial_terms({ c_u });

#ifdef FDCK_DIST_MODE
        disable_output_time_column();
#endif

        set_variable_names({ "c_u" });

        Max_residual = 1.0e-14;

#ifdef FDCK_TRACE_FILE
        trace_file_.open("trace.dat");
        trace_header_ = "t c_u(0) c_u(1)";
        trace_file_ << trace_header_ << '\n';
#endif
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        if(p.ics == "zero")
        {
            // Initial time is zero. We do a linear interpolation from t = 0
            // and uniformly zero BCs to the first actual time point/BC values
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
            time() = p.start_time;

            // set the initial conditions from the polynomial fit
            for(unsigned i = 0; i < n_node_; ++i)
            {
                const double scaled_x = map_space(x(i));

                // get the interpolated value
                const double c_u_interp =
                    utilities::evaluate_polynomial(scaled_x, ic_poly_coeffs_);

                // set to the interpolated value, avoiding negative
                // concentrations by taking the max with zero
                u(0, c_u, i) = std::max(0.0, c_u_interp);
            }
        }
        // Assume p.ics is the name of the ics CSV file
        else
        {
            // Initial time is the "start time" of the experiment
            time() = p.start_time;

            // Open the ICs file
            std::ifstream ics_file(p.ics);

            // If the file opens successfully, read the CSV data and set the
            // ICs by lerping onto the mesh
            if(ics_file)
            {
                auto [m_vec, n_rows, n_cols] =
                    utilities::read_csv_to_flat_vector(ics_file, p.csv_delimiter, p.csv_skip_rows);

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
                FDCK_FATAL("Error opening ICs file \"{}\"! Exiting", p.ics);
                std::exit(1);
            }
        }
    }

    ChemokinesParams p;

#ifdef FDCK_TRACE_FILE
    void trace()
    {
        trace_file_ << time() << " "
                    << u(c_u, 0) << " "
                    << u(c_u, n_node_-1) << '\n'
                    << std::flush;
    }
#endif

private:
#ifdef FDCK_TRACE_FILE
    std::string trace_header_;
    std::ofstream trace_file_;
#endif

    // Coeffs for the polynomial fits of the inlet/outlet BC functions and ICs
    // The polynomials are fitted to the data after scaling the "x" axis to [0,1]
    const std::vector<double> inlet_poly_coeffs_;
    const std::vector<double> outlet_poly_coeffs_;
    const std::vector<double> ic_poly_coeffs_;

public:
    // The lower bound of the BC polynomial fit
    double fit_start_time;

    // The upper bound of the BC polynomial fit
    double fit_end_time;

private:
    // Map the spatial variable from its actual range to the interval [0,1]
    double map_space(double x) const
    {
        return (x - a_)/(b_ - a_);
    }

    // Map the time variable from its actual range to the interval [0,1].
    // This uses the fit start/end times since these define the actual range of
    // times that are valid.
    double map_time(double time) const
    {
        return (time - fit_start_time)/(fit_end_time - fit_start_time);
    }

    void actions_after_timestep() override
    {
#ifdef FDCK_TRACE_FILE
        trace();
#endif
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
            // if the flag is set and we haven't reached fit_start_time, do a
            // linear interpolation between (0,0) and (fit_start_time, c(fit_start_time))
            if(p.ics == "zero" && time() < fit_start_time)
            {
                // get the (interpolated) concentration at t = fit_start_time
                const static double start_c =
                    utilities::evaluate_polynomial(map_time(fit_start_time), inlet_poly_coeffs_);

                // calculate the slope of the line
                const static double m = start_c/fit_start_time;

                a3[c_u] = m*time();
            }
            else
            {
                // get the interpolated value
                const double c_u_interp =
                    utilities::evaluate_polynomial(map_time(time()), inlet_poly_coeffs_);

                // set to the interpolated value, avoiding negative
                // concentrations by taking the max with zero
                a3[c_u] = std::max(0.0, c_u_interp);
            }
        }
        if(b == Boundary::Right)
        {
            // if the flag is set and we haven't reached fit_start_time, do a
            // linear interpolation between (0,0) and (fit_start_time, c(fit_start_time))
            if(p.ics == "zero" && time() < fit_start_time)
            {
                // get the (interpolated) concentration at t = fit_start_time
                const static double start_c =
                    utilities::evaluate_polynomial(map_time(fit_start_time), outlet_poly_coeffs_);

                // calculate the slope of the line
                const static double m = start_c/fit_start_time;

                a3[c_u] = m*time();
            }
            else
            {
                // get the interpolated value
                const double c_u_interp =
                    utilities::evaluate_polynomial(map_time(time()), outlet_poly_coeffs_);

                // set to the interpolated value, avoiding negative
                // concentrations by taking the max with zero
                a3[c_u] = std::max(0.0, c_u_interp);
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

#ifdef FDCK_DIST_MODE
    log::disable_labels();
#endif

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
    else if(dataset == "june2019")
    {
        builder = std::move(data::june_2019_builder);
    }
    else if(dataset == "october2019")
    {
        builder = std::move(data::october_2019_builder);
    }
    else
    {
        FDCK_FATAL("Unrecognised dataset `{}` specified! Exiting", dataset);
        std::exit(1);
    }

    ChemokinesProblem1D problem(n_node, std::move(builder));

    problem.p.u = cf.get<double>("u");
    problem.p.D = cf.get<double>("D");
    problem.p.ics = cf.get_or<std::string>("ics", "polynomial");
    problem.p.csv_delimiter = cf.get_or<char>("csv_delimiter", ',');
    problem.p.csv_skip_rows = cf.get_or<unsigned>("csv_skip_rows", 0);
    problem.p.start_time = cf.get_or<double>("start_time", builder.fit_start_time);

    if(problem.p.ics == "zero" && std::abs(problem.p.start_time) > 1.0e-15)
    {
        FDCK_ERROR("ics = {} but start_time = {}; start_time must be 0.0 if zero ICs are used",
                    problem.p.ics,
                    problem.p.start_time);
    }

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
    std::sprintf(filename, "output_%05u.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

#ifdef FDCK_TRACE_FILE
    problem.trace();
#endif

    unsigned i = 1;

#ifdef _WIN32
    std::cout.sync_with_stdio(false);
#endif

    // timestepping loop
    while(problem.time() <= t_max && problem.time() <= problem.fit_end_time)
    {
        // solve for current timestep
        problem.unsteady_solve(dt);

        if(i % output_interval == 0)
        {
            // output current solution
            FDCK_TRACE("Outputting");
            std::sprintf(filename, "output_%05u.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }

    FDCK_INFO("Reached t > min(t_max = {}, fit_end_time = {}) after performing {} timesteps", t_max, problem.fit_end_time, i-1);
}
