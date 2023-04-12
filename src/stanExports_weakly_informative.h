// Generated by rstantools.  Do not edit by hand.

/*
    baldur is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    baldur is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with baldur3.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_weakly_informative_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_weakly_informative");
    reader.add_event(46, 44, "end", "model_weakly_informative");
    return reader;
}
#include <stan_meta_header.hpp>
class model_weakly_informative
  : public stan::model::model_base_crtp<model_weakly_informative> {
private:
        int N;
        int K;
        int C;
        matrix_d x;
        vector_d y;
        matrix_d c;
        double alpha;
        double beta;
        vector_d u;
        vector_d n_k;
        row_vector_d n_c;
        matrix_d abs_c;
public:
    model_weakly_informative(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_weakly_informative(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_weakly_informative_namespace::model_weakly_informative";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "C", "int", context__.to_vec());
            C = int(0);
            vals_i__ = context__.vals_i("C");
            pos__ = 0;
            C = vals_i__[pos__++];
            current_statement_begin__ = 5;
            validate_non_negative_index("x", "N", N);
            validate_non_negative_index("x", "K", K);
            context__.validate_dims("data initialization", "x", "matrix_d", context__.to_vec(N,K));
            x = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, K);
            vals_r__ = context__.vals_r("x");
            pos__ = 0;
            size_t x_j_2_max__ = K;
            size_t x_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < x_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < x_j_1_max__; ++j_1__) {
                    x(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(N));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("c", "K", K);
            validate_non_negative_index("c", "C", C);
            context__.validate_dims("data initialization", "c", "matrix_d", context__.to_vec(K,C));
            c = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(K, C);
            vals_r__ = context__.vals_r("c");
            pos__ = 0;
            size_t c_j_2_max__ = C;
            size_t c_j_1_max__ = K;
            for (size_t j_2__ = 0; j_2__ < c_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < c_j_1_max__; ++j_1__) {
                    c(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "alpha", "double", context__.to_vec());
            alpha = double(0);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            alpha = vals_r__[pos__++];
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "beta", "double", context__.to_vec());
            beta = double(0);
            vals_r__ = context__.vals_r("beta");
            pos__ = 0;
            beta = vals_r__[pos__++];
            current_statement_begin__ = 10;
            validate_non_negative_index("u", "N", N);
            context__.validate_dims("data initialization", "u", "vector_d", context__.to_vec(N));
            u = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("u");
            pos__ = 0;
            size_t u_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < u_j_1_max__; ++j_1__) {
                u(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            current_statement_begin__ = 14;
            validate_non_negative_index("n_k", "K", K);
            n_k = Eigen::Matrix<double, Eigen::Dynamic, 1>(K);
            stan::math::fill(n_k, DUMMY_VAR__);
            current_statement_begin__ = 15;
            validate_non_negative_index("n_c", "C", C);
            n_c = Eigen::Matrix<double, 1, Eigen::Dynamic>(C);
            stan::math::fill(n_c, DUMMY_VAR__);
            current_statement_begin__ = 16;
            validate_non_negative_index("abs_c", "K", K);
            validate_non_negative_index("abs_c", "C", C);
            abs_c = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(K, C);
            stan::math::fill(abs_c, DUMMY_VAR__);
            stan::math::assign(abs_c,stan::math::fabs(c));
            // execute transformed data statements
            current_statement_begin__ = 17;
            for (int i = 1; i <= K; ++i) {
                current_statement_begin__ = 18;
                stan::model::assign(n_k, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (1 / sum(stan::model::rvalue(x, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), "x"))), 
                            "assigning variable n_k");
            }
            current_statement_begin__ = 20;
            stan::math::assign(n_c, multiply(transpose(n_k), abs_c));
            current_statement_begin__ = 21;
            stan::math::assign(n_c, stan::math::sqrt(n_c));
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 25;
            validate_non_negative_index("mu", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 26;
            num_params_r__ += 1;
            current_statement_begin__ = 27;
            validate_non_negative_index("y_diff", "C", C);
            num_params_r__ += (1 * C);
            current_statement_begin__ = 28;
            validate_non_negative_index("eta", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 29;
            validate_non_negative_index("prior_mu_not", "K", K);
            num_params_r__ += K;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_weakly_informative() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 25;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "K", K);
        context__.validate_dims("parameter initialization", "mu", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu(K);
        size_t mu_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 26;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 27;
        if (!(context__.contains_r("y_diff")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable y_diff missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("y_diff");
        pos__ = 0U;
        validate_non_negative_index("y_diff", "C", C);
        context__.validate_dims("parameter initialization", "y_diff", "double", context__.to_vec(C));
        std::vector<double> y_diff(C, double(0));
        size_t y_diff_k_0_max__ = C;
        for (size_t k_0__ = 0; k_0__ < y_diff_k_0_max__; ++k_0__) {
            y_diff[k_0__] = vals_r__[pos__++];
        }
        size_t y_diff_i_0_max__ = C;
        for (size_t i_0__ = 0; i_0__ < y_diff_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(y_diff[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable y_diff: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 28;
        if (!(context__.contains_r("eta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable eta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("eta");
        pos__ = 0U;
        validate_non_negative_index("eta", "K", K);
        context__.validate_dims("parameter initialization", "eta", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta(K);
        size_t eta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            eta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(eta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable eta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 29;
        if (!(context__.contains_r("prior_mu_not")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable prior_mu_not missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("prior_mu_not");
        pos__ = 0U;
        validate_non_negative_index("prior_mu_not", "K", K);
        context__.validate_dims("parameter initialization", "prior_mu_not", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> prior_mu_not(K);
        size_t prior_mu_not_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < prior_mu_not_j_1_max__; ++j_1__) {
            prior_mu_not(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(prior_mu_not);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable prior_mu_not: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 25;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.vector_constrain(K, lp__);
            else
                mu = in__.vector_constrain(K);
            current_statement_begin__ = 26;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 27;
            std::vector<local_scalar_t__> y_diff;
            size_t y_diff_d_0_max__ = C;
            y_diff.reserve(y_diff_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < y_diff_d_0_max__; ++d_0__) {
                if (jacobian__)
                    y_diff.push_back(in__.scalar_constrain(lp__));
                else
                    y_diff.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 28;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta;
            (void) eta;  // dummy to suppress unused var warning
            if (jacobian__)
                eta = in__.vector_constrain(K, lp__);
            else
                eta = in__.vector_constrain(K);
            current_statement_begin__ = 29;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> prior_mu_not;
            (void) prior_mu_not;  // dummy to suppress unused var warning
            if (jacobian__)
                prior_mu_not = in__.vector_constrain(K, lp__);
            else
                prior_mu_not = in__.vector_constrain(K);
            // transformed parameters
            current_statement_begin__ = 33;
            validate_non_negative_index("mu_diff", "C", C);
            Eigen::Matrix<local_scalar_t__, 1, Eigen::Dynamic> mu_diff(C);
            stan::math::initialize(mu_diff, DUMMY_VAR__);
            stan::math::fill(mu_diff, DUMMY_VAR__);
            stan::math::assign(mu_diff,multiply(transpose(mu), c));
            current_statement_begin__ = 34;
            validate_non_negative_index("sigma_lfc", "C", C);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sigma_lfc(C);
            stan::math::initialize(sigma_lfc, DUMMY_VAR__);
            stan::math::fill(sigma_lfc, DUMMY_VAR__);
            stan::math::assign(sigma_lfc,multiply(sigma, transpose(n_c)));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 33;
            size_t mu_diff_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < mu_diff_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(mu_diff(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu_diff" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu_diff: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 34;
            size_t sigma_lfc_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < sigma_lfc_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(sigma_lfc(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: sigma_lfc" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable sigma_lfc: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 38;
            lp_accum__.add(gamma_log<propto__>(sigma, alpha, beta));
            current_statement_begin__ = 39;
            lp_accum__.add(normal_log<propto__>(eta, 0, 1));
            current_statement_begin__ = 40;
            lp_accum__.add(normal_log<propto__>(prior_mu_not, 0, 10));
            current_statement_begin__ = 41;
            lp_accum__.add(normal_log<propto__>(mu, add(prior_mu_not, multiply(sigma, eta)), sigma));
            current_statement_begin__ = 42;
            lp_accum__.add(normal_log<propto__>(y, multiply(x, mu), multiply(sigma, u)));
            current_statement_begin__ = 43;
            lp_accum__.add(normal_log<propto__>(y_diff, mu_diff, sigma_lfc));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu");
        names__.push_back("sigma");
        names__.push_back("y_diff");
        names__.push_back("eta");
        names__.push_back("prior_mu_not");
        names__.push_back("mu_diff");
        names__.push_back("sigma_lfc");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(C);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(C);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(C);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_weakly_informative_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu = in__.vector_constrain(K);
        size_t mu_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            vars__.push_back(mu(j_1__));
        }
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        std::vector<double> y_diff;
        size_t y_diff_d_0_max__ = C;
        y_diff.reserve(y_diff_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < y_diff_d_0_max__; ++d_0__) {
            y_diff.push_back(in__.scalar_constrain());
        }
        size_t y_diff_k_0_max__ = C;
        for (size_t k_0__ = 0; k_0__ < y_diff_k_0_max__; ++k_0__) {
            vars__.push_back(y_diff[k_0__]);
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta = in__.vector_constrain(K);
        size_t eta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            vars__.push_back(eta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> prior_mu_not = in__.vector_constrain(K);
        size_t prior_mu_not_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < prior_mu_not_j_1_max__; ++j_1__) {
            vars__.push_back(prior_mu_not(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 33;
            validate_non_negative_index("mu_diff", "C", C);
            Eigen::Matrix<double, 1, Eigen::Dynamic> mu_diff(C);
            stan::math::initialize(mu_diff, DUMMY_VAR__);
            stan::math::fill(mu_diff, DUMMY_VAR__);
            stan::math::assign(mu_diff,multiply(transpose(mu), c));
            current_statement_begin__ = 34;
            validate_non_negative_index("sigma_lfc", "C", C);
            Eigen::Matrix<double, Eigen::Dynamic, 1> sigma_lfc(C);
            stan::math::initialize(sigma_lfc, DUMMY_VAR__);
            stan::math::fill(sigma_lfc, DUMMY_VAR__);
            stan::math::assign(sigma_lfc,multiply(sigma, transpose(n_c)));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t mu_diff_j_1_max__ = C;
                for (size_t j_1__ = 0; j_1__ < mu_diff_j_1_max__; ++j_1__) {
                    vars__.push_back(mu_diff(j_1__));
                }
                size_t sigma_lfc_j_1_max__ = C;
                for (size_t j_1__ = 0; j_1__ < sigma_lfc_j_1_max__; ++j_1__) {
                    vars__.push_back(sigma_lfc(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_weakly_informative";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        size_t y_diff_k_0_max__ = C;
        for (size_t k_0__ = 0; k_0__ < y_diff_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_diff" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t eta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t prior_mu_not_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < prior_mu_not_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "prior_mu_not" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_diff_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < mu_diff_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_diff" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t sigma_lfc_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < sigma_lfc_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sigma_lfc" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t mu_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        size_t y_diff_k_0_max__ = C;
        for (size_t k_0__ = 0; k_0__ < y_diff_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_diff" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t eta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t prior_mu_not_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < prior_mu_not_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "prior_mu_not" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_diff_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < mu_diff_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_diff" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t sigma_lfc_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < sigma_lfc_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sigma_lfc" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_weakly_informative_namespace::model_weakly_informative stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
