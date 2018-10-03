#include <advection_diffusion_reaction_problem.h>
#include <utilities.h>

namespace mjrfd
{
    /// Constructor
    AdvectionDiffusionReactionProblem::AdvectionDiffusionReactionProblem(
        const unsigned n_var,
        const unsigned n_node,
        const double a,
        const double b) :
        Problem(n_var, n_node),
        n_node_(n_node),
        a_(a),
        b_(b),
        dx_((b - a)/(n_node-1)),
        cn_theta_(1.0),
        left_bc_(n_var, false),
        right_bc_(n_var, false),
        spatial_terms_(n_var, false),
        output_time_column_(true),
        var_names_(n_var)
    {
        // Set the default variable names
        for(unsigned var = 0; var < n_var_; ++var)
        {
            var_names_[var] = "c" + std::to_string(var);
        }
    }

    /// Destructor
    AdvectionDiffusionReactionProblem::~AdvectionDiffusionReactionProblem()
    {
    }

    void AdvectionDiffusionReactionProblem::make_steady()
    {
        Problem::make_steady();
        cn_theta_backup_ = cn_theta_;
        cn_theta_ = 1.0;
    }

    void AdvectionDiffusionReactionProblem::make_unsteady()
    {
        Problem::make_unsteady();
        cn_theta_ = cn_theta_backup_;
    }

    /// Output the solution
    void AdvectionDiffusionReactionProblem::output(std::ostream &out) const
    {
        if(output_time_column_ == true)
        {
            out << "t ";
        }

        out << "x ";

        for(unsigned var = 0; var < n_var_; ++var)
        {
            out << '\"' << var_names_[var] << '\"';
            out << (var < (n_var_-1) ? ' ' : '\n');
        }

        for(unsigned i = 0; i < n_node_; ++i)
        {
            if(output_time_column_ == true)
            {
                out << time() << ' ';
            }

            out << x(i) << ' ';

            for(unsigned var = 0; var < n_var_; ++var)
            {
                out << u(0, var, i);
                out << (var < (n_var_-1) ? ' ' : '\n');
            }
        }
    }

    void AdvectionDiffusionReactionProblem::output_exact(std::ostream &out) const
    {
        if(output_time_column_ == true)
        {
            out << "t ";
        }

        out << "x ";

        for(unsigned var = 0; var < n_var_; ++var)
        {
            out << '\"' << var_names_[var] << '\"';
            out << (var < (n_var_-1) ? ' ' : '\n');
        }

        double x;
        std::vector<double> sol(2);

        for(unsigned i = 0; i < n_node_; ++i)
        {
            if(output_time_column_ == true)
            {
                out << time() << ' ';
            }

            x = this->x(i);
            out << x << ' ';

            exact_solution(time(), x, sol);

            for(unsigned var = 0; var < n_var_; ++var)
            {
                out << sol[var];
                out << (var < (n_var_-1) ? ' ' : '\n');
            }
        }
    }

    /// Calculate the residual vector
    void AdvectionDiffusionReactionProblem::calculate_residual(Eigen::VectorXd &residual)
    {
        // Set the residuals to zero
        residual.setZero();

        // Storage for the diffusion coefficients
        std::vector<double> d(n_var_);
        std::vector<double> d_plus(n_var_);
        std::vector<double> d_minus(n_var_);

        // Get the left boundary condition coefficients
        std::vector<double> a1_left(n_var_);
        std::vector<double> a2_left(n_var_);
        std::vector<double> a3_left(n_var_);
        get_bc(Boundary::Left, a1_left, a2_left, a3_left);

        // Get the right boundary condition coefficients
        std::vector<double> a1_right(n_var_);
        std::vector<double> a2_right(n_var_);
        std::vector<double> a3_right(n_var_);
        get_bc(Boundary::Right, a1_right, a2_right, a3_right);

        // Storage for the advection velocity at a node
        double v = 0;

        std::vector<double> r(n_var_);
        std::vector<double> r_old(n_var_);

        // Storage for u (all vars) at a node
        std::vector<double> u_at_node(n_var_);
        std::vector<double> u_old_at_node(n_var_);

        // Loop over the nodes
        for(unsigned i = 0; i < n_node_; ++i)
        {
            for(unsigned var = 0; var < n_var_; ++var)
            {
                u_at_node[var]     = u(0, var, i);
                u_old_at_node[var] = u(1, var, i);
            }

            // Get the reaction term at the current node
            get_r(u_at_node, r);
            get_r(u_old_at_node, r_old);

            // Loop over the variables
            for(unsigned var = 0; var < n_var_; ++var)
            {
                // Calculate the index of the current unknown
                const unsigned index = var*n_node_ + i;

                // Boundary conditions
                if(i == 0 && left_bc_[var] == true)
                {
                    residual(index) += (a1_left[var]*u(0, var, i) - a3_left[var])*dx_*dx_;

                    for(const auto& [j, w] : stencil::forward_1::weights)
                    {
                        residual(index) += w*a2_left[var]*u(0, var, i+j)*dx_;
                    }
                }
                else if(i == n_node_-1 && right_bc_[var] == true)
                {
                    residual(index) += (a1_right[var]*u(0, var, i) - a3_right[var])*dx_*dx_;

                    for(const auto& [j, w] : stencil::backward_1::weights)
                    {
                        residual(index) += w*a2_right[var]*u(0, var, i+j)*dx_;
                    }
                }
                else
                {
                    // Time derivatives
                    if(is_steady() == false)
                    {
                        residual(index) += dx_*dx_/dt_*(u(0, var, i) - u(1, var, i));
                    }

                    // Add the spatial terms to the residual for this variable
                    // if they are enabled
                    if(spatial_terms_[var] == true)
                    {
                        // Diffusion

                        // TODO same probably goes for diffusion as noted below for
                        // velocities regarding old and new values and cn_theta

                        // Get the diffusion coefficients at the current node and
                        // the ones surrounding it so we can interpolate
                        // TODO this is a bit inefficient because it gets the
                        // coeffs for all variables every time

                        // Double check that we haven't got to this point with
                        // the first or last node, since the following would
                        // explode
                        assert(i != 0 && i != (n_node_-1));

                        // the i+-1 here is a bit hardcoded, but we know that these
                        // indices will always be appropriate and in range since
                        // we're using central differences in the outer loop below
                        get_d(0, i,   d);
                        get_d(0, i+1, d_plus);
                        get_d(0, i-1, d_minus);

                        // Here we do a nested central difference with step
                        // size dx/2 (which is where the factor of 4 comes
                        // from) Loop over the outer stencil points
                        for(const auto& [j, w] : stencil::central_1::weights)
                        {
                            // linearly interpolate the diffusion coeff to
                            // the mid-node locations
                            double d_lerp_var = 0.0;

                            // interpolate according to which side of the
                            // outer stencil point the inner stencil point
                            // falls on
                            if(j > 0)
                            {
                                d_lerp_var = 0.5*(d[var] + d_plus[var]);
                            }
                            else
                            {
                                d_lerp_var = 0.5*(d_minus[var] + d[var]);
                            }

                            double du_dx_0 = 0.0;
                            double du_dx_1 = 0.0;

                            // Loop over the inner stencil points
                            for(const auto& [k, w2] : stencil::central_1::weights)
                            {
                                // This integer division is correct since
                                // j+k is always +-2 or 0 (the +-1 from the
                                // individual stencils either combine or
                                // cancel)
                                du_dx_0 += w2*u(0, var, i + (j+k)/2);
                                du_dx_1 += w2*u(1, var, i + (j+k)/2);
                            }

                            // Only add diffusion terms if the coefficient is strictly positive
                            if(d_lerp_var > 0.0)
                            {
                                residual(index) += -4.0*d_lerp_var*w*(cn_theta_*du_dx_0 + (1.0-cn_theta_)*du_dx_1);
                            }
                        }

                        // Advection (or chemotaxis etc)

                        // Get the velocity at the location of the target node
                        double v_at_node = 0.0;
                        get_v(i, var, v_at_node);

                        // TODO when cn_theta != 1 the velocity residuals are
                        // incorrect because we use the same v for both the
                        // current and previous timestep terms - must add time
                        // arg to get_v and use it here

                        // Loop over the stencil points, using the upwind
                        // stencil helper, and get the offset j (relative to i)
                        // and the weight w
                        for(const auto& [j, w] : upwind_stencil_weights(i, v_at_node))
                        {
                            // Get the velocity at the current stencil point
                            get_v(i+j, var, v);

                            // Sum the contributions from the stencil points
                            // Only add advection terms if the velocity
                            // magnitude is strictly positive
                            if(std::abs(v) > 0.0)
                            {
                                residual(index) += dx_*v*w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
                            }
                        }
                    }

                    // Reaction
                    residual(index) += -dx_*dx_*(cn_theta_*r[var] + (1.0-cn_theta_)*r_old[var]);
                }
            }
        }
    }

    /// Calculate the jacobian matrix
    void AdvectionDiffusionReactionProblem::calculate_jacobian(
        Eigen::SparseMatrix<double> &jacobian)
    {
        jacobian.setZero();

        // Vector of triplets for constructing the sparse jacobian matrix
        std::vector<T> triplet_list;

        // TODO check this! It appears to be ~6x the actual number of entries
        // Reserve the correct number of entries
        triplet_list.reserve(n_node_*n_var_*(15*n_var_ + 7));

        // Storage for the diffusion coefficients
        std::vector<double> d(n_var_);
        std::vector<double> d_plus(n_var_);
        std::vector<double> d_minus(n_var_);

        // Storage for the derivatives w.r.t. dofs of the diffusion coefficients
        std::vector<std::vector<double>> dd_du(n_var_,       std::vector<double>(n_var_));
        std::vector<std::vector<double>> dd_plus_du(n_var_,  std::vector<double>(n_var_));
        std::vector<std::vector<double>> dd_minus_du(n_var_, std::vector<double>(n_var_));

        // Get the left boundary condition coefficients
        std::vector<double> a1_left(n_var_);
        std::vector<double> a2_left(n_var_);
        std::vector<double> a3_left(n_var_);
        get_bc(Boundary::Left, a1_left, a2_left, a3_left);

        // Get the right boundary condition coefficients
        std::vector<double> a1_right(n_var_);
        std::vector<double> a2_right(n_var_);
        std::vector<double> a3_right(n_var_);
        get_bc(Boundary::Right, a1_right, a2_right, a3_right);

        // Get the left boundary condition coefficient derivatives
        std::vector<std::vector<double>> da1_du_left(n_var_, std::vector<double>(n_var_));
        std::vector<std::vector<double>> da2_du_left(n_var_, std::vector<double>(n_var_));
        std::vector<std::vector<double>> da3_du_left(n_var_, std::vector<double>(n_var_));

        // Get the right boundary condition coefficient derivatives
        std::vector<std::vector<double>> da1_du_right(n_var_, std::vector<double>(n_var_));
        std::vector<std::vector<double>> da2_du_right(n_var_, std::vector<double>(n_var_));
        std::vector<std::vector<double>> da3_du_right(n_var_, std::vector<double>(n_var_));

        // Storage for the advection velocity at a node
        double v = 0.0;
        double dv_du = 0.0;

        std::vector<double> r(n_var_);
        std::vector<std::vector<double>> dr_du(n_var_, std::vector<double>(n_var_));

        // Storage for u (all vars) at a node
        std::vector<double> u_at_node(n_var_);

        // Loop over the nodes
        for(unsigned i = 0; i < n_node_; ++i)
        {
            for(unsigned var = 0; var < n_var_; ++var)
            {
                u_at_node[var] = u(0, var, i);
            }

            // Get the reaction term and derivatives at the current node
            get_r(u_at_node, r);
            get_dr_du(u_at_node, dr_du);

            // Loop over the variables
            for(unsigned var = 0; var < n_var_; ++var)
            {
                // Calculate the index of the current unknown
                const unsigned index = var*n_node_ + i;

                if(i == 0 && left_bc_[var] == true)
                {
                    // First term
                    if(std::abs(a1_left[var]) > 0.0)
                    {
                        triplet_list.emplace_back(index, index, a1_left[var]*dx_*dx_);
                    }

                    if(std::abs(a2_left[var]) > 0.0)
                    {
                        for(const auto& [j, w] : stencil::forward_1::weights)
                        {
                            triplet_list.emplace_back(index, index+j, w*a2_left[var]*dx_);
                        }
                    }

                    // Second term

                    // Loop over the 3 nodes that the BC coeffs might depend
                    // on. This covers const coeffs and coeffs that involve
                    // first derivatives (forward differences) of other
                    // variables. Can increase range if needed.

                    // NOTE: i == 0 always here - use i rather than 0 to
                    // increase similarity between inlet/outlet BC and bulk
                    // equations
                    for(unsigned k = 0; k <= 2; ++k)
                    {
                        const unsigned i2 = i+k;

                        // Get the BC coeff derivatives wrt all vars at node i2
                        get_dbc_du(Boundary::Left, i2,
                                   da1_du_left, da2_du_left, da3_du_left);

                        // Loop over the variables again
                        for(unsigned var2 = 0; var2 < n_var_; ++var2)
                        {
                            const unsigned index2 = var2*n_node_ + i2;

                            if(std::abs(da1_du_left[var][var2]) > 0.0)
                            {
                                triplet_list.emplace_back(index, index2, da1_du_left[var][var2]*u(0, var, i)*dx_*dx_);
                            }

                            if(std::abs(da2_du_left[var][var2]) > 0.0)
                            {
                                double deriv = 0.0;

                                for(const auto& [j, w] : stencil::forward_1::weights)
                                {
                                    deriv += w*u(0, var, i+j)*dx_;
                                }

                                triplet_list.emplace_back(index, index2, da2_du_left[var][var2]*deriv);
                            }

                            if(std::abs(da3_du_left[var][var2]) > 0.0)
                            {
                                triplet_list.emplace_back(index, index2, -da3_du_left[var][var2]*dx_*dx_);
                            }
                        }
                    }
                }
                else if(i == n_node_-1 && right_bc_[var] == true)
                {
                    // First term
                    if(std::abs(a1_right[var]) > 0.0)
                    {
                        triplet_list.emplace_back(index, index, a1_right[var]*dx_*dx_);
                    }

                    if(std::abs(a2_right[var]) > 0.0)
                    {
                        for(const auto& [j, w] : stencil::backward_1::weights)
                        {
                            triplet_list.emplace_back(index, index+j, w*a2_right[var]*dx_);
                        }
                    }

                    // Second term
                    for(int k = -2; k <= 0; ++k)
                    {
                        const unsigned i2 = i+k;

                        // Get the BC coeff derivatives wrt all vars at node i2
                        get_dbc_du(Boundary::Right, i2,
                                   da1_du_right, da2_du_right, da3_du_right);

                        // Loop over the variables again
                        for(unsigned var2 = 0; var2 < n_var_; ++var2)
                        {
                            const unsigned index2 = var2*n_node_ + i2;

                            if(std::abs(da1_du_right[var][var2]) > 0.0)
                            {
                                triplet_list.emplace_back(index, index2, da1_du_right[var][var2]*u(0, var, i)*dx_*dx_);
                            }

                            if(std::abs(da2_du_right[var][var2]) > 0.0)
                            {
                                double deriv = 0.0;

                                for(const auto& [j, w] : stencil::forward_1::weights)
                                {
                                    deriv += w*u(0, var, i+j)*dx_;
                                }

                                triplet_list.emplace_back(index, index2, da2_du_right[var][var2]*deriv);
                            }

                            if(std::abs(da3_du_right[var][var2]) > 0.0)
                            {
                                triplet_list.emplace_back(index, index2, -da3_du_right[var][var2]*dx_*dx_);
                            }
                        }
                    }
                }
                else
                {
                    // Time derivatives
                    if(is_steady() == false)
                    {
                        triplet_list.emplace_back(index, index, dx_*dx_/dt_);
                    }

                    // Add the spatial terms to the residual for this variable
                    // if they are enabled
                    if(spatial_terms_[var] == true)
                    {
                        // Diffusion

                        // Get the diffusion coeffs and their derivatives at the
                        // required locations
                        get_d(0, i,   d);
                        get_d(0, i+1, d_plus);
                        get_d(0, i-1, d_minus);

                        get_dd_du(0, i,   dd_du);
                        get_dd_du(0, i+1, dd_minus_du);
                        get_dd_du(0, i-1, dd_plus_du);

                        // Loop over the outer stencil points
                        for(const auto& [j, w] : stencil::central_1::weights)
                        {
                            // linearly interpolate the diffusion coeff to the
                            // mid-node locations
                            double d_lerp_var = 0.0;

                            // interpolate according to which side of the outer
                            // stencil point the inner stencil point falls on
                            if(j > 0)
                            {
                                d_lerp_var = 0.5*(d[var] + d_plus[var]);
                            }
                            else
                            {
                                d_lerp_var = 0.5*(d_minus[var] + d[var]);
                            }

                            // Loop over the inner stencil points
                            for(const auto& [k, w2] : stencil::central_1::weights)
                            {
                                // This integer division is correct since j+k is
                                // always +-2 or 0 (the +-1 from the individual
                                // stencils either combine or cancel)

                                // First term - add a jacobian entry for each of
                                // the unknowns that appear in the nested stencils
                                if(d_lerp_var > 0.0)
                                {
                                    triplet_list.emplace_back(index, index+(j+k)/2, -4.0*d_lerp_var*w*w2*cn_theta_);
                                }

                                // Second term - terms with derivatives of the
                                // diffusion coefficients
                                for(unsigned var2 = 0; var2 < n_var_; ++var2)
                                {
                                    if(std::abs(dd_du[var][var2]) > 0.0)
                                    {
                                        const unsigned index2 = var2*n_node_ + i;

                                        triplet_list.emplace_back(index, index2, -4.0*0.5*dd_du[var][var2]*w*w2*cn_theta_*u(0, var, i+(j+k)/2));
                                    }

                                    const unsigned index2 = var2*n_node_ + i + j;

                                    if(j > 0)
                                    {
                                        if(std::abs(dd_plus_du[var][var2]) > 0.0)
                                        {
                                            triplet_list.emplace_back(index, index2, -4.0*0.5*dd_plus_du[var][var2]*w*w2*cn_theta_*u(0, var, i+(j+k)/2));
                                        }
                                    }
                                    else // j < 0
                                    {
                                        if(std::abs(dd_minus_du[var][var2]) > 0.0)
                                        {
                                            triplet_list.emplace_back(index, index2, -4.0*0.5*dd_minus_du[var][var2]*w*w2*cn_theta_*u(0, var, i+(j+k)/2));
                                        }
                                    }
                                } // end of loop over var2
                            } // end of loop over inner stencil
                        } // end of loop over outer stencil

                        // Advection (or chemotaxis etc)
                        double v_at_node = 0.0;
                        get_v(i, var, v_at_node);

                        for(const auto& [j, w] : upwind_stencil_weights(i, v_at_node))
                        {
                            get_v(i+j, var, v);

                            // If the velocity is nonzero, add the term to the
                            // jacobian. Note that dv_du may be != 0 even when
                            // v == 0 (I think), so we always check it
                            // separately. Even if this is wrong, it doesn't
                            // hurt
                            if(std::abs(v) > 0.0)
                            {
                                triplet_list.emplace_back(index, index+j, dx_*w*cn_theta_*v);
                            }

                            // Loop over the variables again
                            for(unsigned var2 = 0; var2 < n_var_; ++var2)
                            {
                                // Loop over the 5 stencil points that might be
                                // used in any reasonable definition of get_v
                                // This should cover forward/backward/central
                                // first derivative stencils
                                for(int k = -2; k <= 2; ++k)
                                {
                                    int i2 = i+j+k;

                                    if(i2 < 0 || i2 >= static_cast<int>(n_node_))
                                    {
                                        // if the index we're trying is out of
                                        // the possible range then skip it
                                        continue;
                                    }

                                    const unsigned index2 = var2*n_node_ + i2;

                                    get_dv_du(i+j, var, i2, var2, dv_du);

                                    if(std::abs(dv_du) > 0.0)
                                    {
                                        triplet_list.emplace_back(index, index2, dx_*w*cn_theta_*u(0, var, i+j)*dv_du);
                                    }
                                }
                            }
                        }
                    }

                    // Reaction

                    // Loop over the variables again
                    for(unsigned var2 = 0; var2 < n_var_; ++var2)
                    {
                        const unsigned index2 = var2*n_node_ + i;

                        if(std::abs(dr_du[var][var2]) > 0.0)
                        {
                            triplet_list.emplace_back(index, index2, -dx_*dx_*cn_theta_*dr_du[var][var2]);
                        }
                    }
                }
            }
        }

        jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian.makeCompressed();
    }

    void AdvectionDiffusionReactionProblem::enable_bc(Boundary b, const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            if(b == Boundary::Left)
            {
                left_bc_[var] = true;
            }
            else if(b == Boundary::Right)
            {
                right_bc_[var] = true;
            }
        }
    }

    void AdvectionDiffusionReactionProblem::disable_bc(Boundary b, const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            if(b == Boundary::Left)
            {
                left_bc_[var] = false;
            }
            else if(b == Boundary::Right)
            {
                right_bc_[var] = false;
            }
        }
    }

    void AdvectionDiffusionReactionProblem::enable_spatial_terms(const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            spatial_terms_[var] = true;
        }
    }

    void AdvectionDiffusionReactionProblem::disable_spatial_terms(const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            spatial_terms_[var] = false;
        }
    }

    void AdvectionDiffusionReactionProblem::enable_output_time_column()
    {
        output_time_column_ = true;
    }

    void AdvectionDiffusionReactionProblem::disable_output_time_column()
    {
        output_time_column_ = false;
    }

    void AdvectionDiffusionReactionProblem::set_variable_names(const std::vector<std::string> &var_names)
    {
        assert(var_names.size() == n_var_);
        var_names_ = var_names;
    }

    double AdvectionDiffusionReactionProblem::integrate_solution(const unsigned var) const
    {
        assert(var >= 0 && var <= n_var_-1);

        return utilities::integrate_trapezium(u_[0].segment(var*n_node_, n_node_), dx_);
    }

    double AdvectionDiffusionReactionProblem::a() const
    {
        return a_;
    }

    double AdvectionDiffusionReactionProblem::b() const
    {
        return b_;
    }

    void AdvectionDiffusionReactionProblem::exact_solution(const double,
                                                           const double,
                                                           std::vector<double> &sol) const
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            sol[var] = 0.0;
        }
    }

    void AdvectionDiffusionReactionProblem::get_dbc_du(Boundary,
                    const unsigned,
                    std::vector<std::vector<double>> &da1_du,
                    std::vector<std::vector<double>> &da2_du,
                    std::vector<std::vector<double>> &da3_du) const
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

    void AdvectionDiffusionReactionProblem::get_dd_du(const unsigned,
                                                      const unsigned,
                                                      std::vector<std::vector<double>> &dd_du) const
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dd_du[var][var2] = 0.0;
            }
        }
    }

    void AdvectionDiffusionReactionProblem::get_dv_du(const unsigned,
                                                      const unsigned,
                                                      const unsigned,
                                                      const unsigned,
                                                      double &dv_du) const
    {
        dv_du = 0.0;
    }

    void AdvectionDiffusionReactionProblem::get_dr_du(const std::vector<double> &,
                                                      std::vector<std::vector<double>> &dr_du) const
    {
        for(unsigned var = 0; var < n_var_; ++var)
        {
            for(unsigned var2 = 0; var2 < n_var_; ++var2)
            {
                dr_du[var][var2] = 0.0;
            }
        }
    }

    const std::unordered_map<int, double>&
    AdvectionDiffusionReactionProblem::upwind_stencil_weights(const unsigned i,
                                                              const double v) const
    {
        assert(i >= 0 && i <= n_node_-1);

        //if(v > 2.0)
        //{
            //if(i >= 2)
                //return stencil::backward_1::weights;
            //else if(i == 1)
                //return stencil::central_1::weights;
            //else if(i == 0)
                //return stencil::forward_1::weights;
        //}
        //else if(v < -2.0)
        //{
            //if(i <= n_node_-3)
                //return stencil::forward_1::weights;
            //else if(i == n_node_-2)
                //return stencil::central_1::weights;
            //else if(i == n_node_-1)
                //return stencil::backward_1::weights;
        //}
        //else
        //{
            //if(i == 0)
                //return stencil::forward_1::weights;
            //else if(i == n_node_-1)
                //return stencil::backward_1::weights;
            //else
                //return stencil::central_1::weights;
        //}

        if(v > 0.0)
        {
            if(i >= 2)
                return stencil::backward_1::weights;
            else if(i == 1)
                return stencil::central_1::weights;
            else if(i == 0)
                return stencil::forward_1::weights;
        }
        else
        {
            if(i <= n_node_-3)
                return stencil::forward_1::weights;
            else if(i == n_node_-2)
                return stencil::central_1::weights;
            else if(i == n_node_-1)
                return stencil::backward_1::weights;
        }

        // This point is never reached but exit here just in case to remove a
        // compiler warning
        std::exit(1);
    }

    const std::unordered_map<int, double>&
    AdvectionDiffusionReactionProblem::central_1_stencil_weights(unsigned i) const
    {
        assert(i >= 0 && i <= n_node_-1);

        if(i == 0)
            return stencil::forward_1::weights;
        else if(i == n_node_-1)
            return stencil::backward_1::weights;
        else
            return stencil::central_1::weights;
    }

    double AdvectionDiffusionReactionProblem::x(const unsigned i) const
    {
        return a_ + static_cast<double>(i)*dx_;
    }
}
