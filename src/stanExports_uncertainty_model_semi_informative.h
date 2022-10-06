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
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-1-g67504470
#include <stan/model/model_header.hpp>
namespace model_uncertainty_model_semi_informative_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'string', line 15, column 2 to column 15)",
                                                      " (in 'string', line 16, column 2 to column 22)",
                                                      " (in 'string', line 17, column 2 to column 17)",
                                                      " (in 'string', line 18, column 2 to column 16)",
                                                      " (in 'string', line 19, column 2 to column 25)",
                                                      " (in 'string', line 22, column 2 to column 20)",
                                                      " (in 'string', line 23, column 2 to column 36)",
                                                      " (in 'string', line 26, column 2 to column 42)",
                                                      " (in 'string', line 27, column 2 to column 30)",
                                                      " (in 'string', line 28, column 2 to column 46)",
                                                      " (in 'string', line 29, column 2 to column 57)",
                                                      " (in 'string', line 30, column 2 to column 41)",
                                                      " (in 'string', line 31, column 2 to column 40)",
                                                      " (in 'string', line 2, column 2 to column 17)",
                                                      " (in 'string', line 3, column 2 to column 17)",
                                                      " (in 'string', line 4, column 2 to column 8)",
                                                      " (in 'string', line 5, column 9 to column 10)",
                                                      " (in 'string', line 5, column 12 to column 13)",
                                                      " (in 'string', line 5, column 2 to column 17)",
                                                      " (in 'string', line 6, column 9 to column 10)",
                                                      " (in 'string', line 6, column 2 to column 14)",
                                                      " (in 'string', line 7, column 8 to column 9)",
                                                      " (in 'string', line 7, column 2 to column 14)",
                                                      " (in 'string', line 8, column 2 to column 13)",
                                                      " (in 'string', line 9, column 2 to column 18)",
                                                      " (in 'string', line 10, column 9 to column 10)",
                                                      " (in 'string', line 10, column 2 to column 14)",
                                                      " (in 'string', line 11, column 9 to column 10)",
                                                      " (in 'string', line 11, column 2 to column 19)",
                                                      " (in 'string', line 12, column 2 to column 20)",
                                                      " (in 'string', line 15, column 9 to column 10)",
                                                      " (in 'string', line 17, column 14 to column 15)",
                                                      " (in 'string', line 18, column 9 to column 10)",
                                                      " (in 'string', line 19, column 9 to column 10)",
                                                      " (in 'string', line 22, column 9 to column 10)"};
#include <stan_meta_header.hpp>
class model_uncertainty_model_semi_informative final : public model_base_crtp<model_uncertainty_model_semi_informative> {
private:
  int N;
  int K;
  int C;
  Eigen::Matrix<double, -1, -1> x;
  Eigen::Matrix<double, -1, 1> y;
  std::vector<std::vector<int>> c;
  double alpha;
  double beta_gamma;
  Eigen::Matrix<double, -1, 1> u;
  Eigen::Matrix<double, -1, 1> mu_not;
  double sigma_mu_not;
 
public:
  ~model_uncertainty_model_semi_informative() { }
  
  inline std::string model_name() const final { return "model_uncertainty_model_semi_informative"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-1-g67504470", "stancflags = "};
  }
  
  
  model_uncertainty_model_semi_informative(stan::io::var_context& context__,
                                           unsigned int random_seed__ = 0,
                                           std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_uncertainty_model_semi_informative_namespace::model_uncertainty_model_semi_informative";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 14;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 14;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 14;
      current_statement__ = 14;
      check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 15;
      context__.validate_dims("data initialization","K","int",
          context__.to_vec());
      K = std::numeric_limits<int>::min();
      
      current_statement__ = 15;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 15;
      current_statement__ = 15;
      check_greater_or_equal(function__, "K", K, 0);
      current_statement__ = 16;
      context__.validate_dims("data initialization","C","int",
          context__.to_vec());
      C = std::numeric_limits<int>::min();
      
      current_statement__ = 16;
      C = context__.vals_i("C")[(1 - 1)];
      current_statement__ = 17;
      validate_non_negative_index("x", "N", N);
      current_statement__ = 18;
      validate_non_negative_index("x", "K", K);
      current_statement__ = 19;
      context__.validate_dims("data initialization","x","double",
          context__.to_vec(N, K));
      x = Eigen::Matrix<double, -1, -1>(N, K);
      stan::math::fill(x, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> x_flat__;
        current_statement__ = 19;
        assign(x_flat__, nil_index_list(), context__.vals_r("x"),
          "assigning variable x_flat__");
        current_statement__ = 19;
        pos__ = 1;
        current_statement__ = 19;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 19;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 19;
            assign(x,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              x_flat__[(pos__ - 1)], "assigning variable x");
            current_statement__ = 19;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 20;
      validate_non_negative_index("y", "N", N);
      current_statement__ = 21;
      context__.validate_dims("data initialization","y","double",
          context__.to_vec(N));
      y = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(y, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> y_flat__;
        current_statement__ = 21;
        assign(y_flat__, nil_index_list(), context__.vals_r("y"),
          "assigning variable y_flat__");
        current_statement__ = 21;
        pos__ = 1;
        current_statement__ = 21;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 21;
          assign(y, cons_list(index_uni(sym1__), nil_index_list()),
            y_flat__[(pos__ - 1)], "assigning variable y");
          current_statement__ = 21;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 22;
      validate_non_negative_index("c", "C", C);
      current_statement__ = 23;
      context__.validate_dims("data initialization","c","int",
          context__.to_vec(C, 2));
      c = std::vector<std::vector<int>>(C, std::vector<int>(2, std::numeric_limits<int>::min()));
      
      {
        std::vector<int> c_flat__;
        current_statement__ = 23;
        assign(c_flat__, nil_index_list(), context__.vals_i("c"),
          "assigning variable c_flat__");
        current_statement__ = 23;
        pos__ = 1;
        current_statement__ = 23;
        for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
          current_statement__ = 23;
          for (int sym2__ = 1; sym2__ <= C; ++sym2__) {
            current_statement__ = 23;
            assign(c,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              c_flat__[(pos__ - 1)], "assigning variable c");
            current_statement__ = 23;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 24;
      context__.validate_dims("data initialization","alpha","double",
          context__.to_vec());
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 24;
      alpha = context__.vals_r("alpha")[(1 - 1)];
      current_statement__ = 25;
      context__.validate_dims("data initialization","beta_gamma","double",
          context__.to_vec());
      beta_gamma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 25;
      beta_gamma = context__.vals_r("beta_gamma")[(1 - 1)];
      current_statement__ = 26;
      validate_non_negative_index("u", "N", N);
      current_statement__ = 27;
      context__.validate_dims("data initialization","u","double",
          context__.to_vec(N));
      u = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(u, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> u_flat__;
        current_statement__ = 27;
        assign(u_flat__, nil_index_list(), context__.vals_r("u"),
          "assigning variable u_flat__");
        current_statement__ = 27;
        pos__ = 1;
        current_statement__ = 27;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 27;
          assign(u, cons_list(index_uni(sym1__), nil_index_list()),
            u_flat__[(pos__ - 1)], "assigning variable u");
          current_statement__ = 27;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 28;
      validate_non_negative_index("mu_not", "K", K);
      current_statement__ = 29;
      context__.validate_dims("data initialization","mu_not","double",
          context__.to_vec(K));
      mu_not = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(mu_not, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> mu_not_flat__;
        current_statement__ = 29;
        assign(mu_not_flat__, nil_index_list(), context__.vals_r("mu_not"),
          "assigning variable mu_not_flat__");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 29;
          assign(mu_not, cons_list(index_uni(sym1__), nil_index_list()),
            mu_not_flat__[(pos__ - 1)], "assigning variable mu_not");
          current_statement__ = 29;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 30;
      context__.validate_dims("data initialization","sigma_mu_not","double",
          context__.to_vec());
      sigma_mu_not = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 30;
      sigma_mu_not = context__.vals_r("sigma_mu_not")[(1 - 1)];
      current_statement__ = 31;
      validate_non_negative_index("mu", "K", K);
      current_statement__ = 32;
      validate_non_negative_index("y_diff", "C", C);
      current_statement__ = 33;
      validate_non_negative_index("eta", "K", K);
      current_statement__ = 34;
      validate_non_negative_index("prior_mu_not", "K", K);
      current_statement__ = 35;
      validate_non_negative_index("mu_diff", "C", C);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += K;
      num_params_r__ += 1;
      num_params_r__ += C;
      num_params_r__ += K;
      num_params_r__ += K;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_uncertainty_model_semi_informative_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      Eigen::Matrix<local_scalar_t__, -1, 1> mu;
      mu = Eigen::Matrix<local_scalar_t__, -1, 1>(K);
      stan::math::fill(mu, DUMMY_VAR__);
      
      current_statement__ = 1;
      mu = in__.vector(K);
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0, lp__);
      } else {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0);
      }
      std::vector<local_scalar_t__> y_diff;
      y_diff = std::vector<local_scalar_t__>(C, DUMMY_VAR__);
      
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        current_statement__ = 3;
        assign(y_diff, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable y_diff");}
      Eigen::Matrix<local_scalar_t__, -1, 1> eta;
      eta = Eigen::Matrix<local_scalar_t__, -1, 1>(K);
      stan::math::fill(eta, DUMMY_VAR__);
      
      current_statement__ = 4;
      eta = in__.vector(K);
      Eigen::Matrix<local_scalar_t__, -1, 1> prior_mu_not;
      prior_mu_not = Eigen::Matrix<local_scalar_t__, -1, 1>(K);
      stan::math::fill(prior_mu_not, DUMMY_VAR__);
      
      current_statement__ = 5;
      prior_mu_not = in__.vector(K);
      Eigen::Matrix<local_scalar_t__, -1, 1> mu_diff;
      mu_diff = Eigen::Matrix<local_scalar_t__, -1, 1>(C);
      stan::math::fill(mu_diff, DUMMY_VAR__);
      
      current_statement__ = 7;
      assign(mu_diff, nil_index_list(),
        subtract(
          rvalue(mu,
            cons_list(
              index_multi(rvalue(c,
                            cons_list(index_omni(),
                              cons_list(index_uni(1), nil_index_list())),
                            "c")), nil_index_list()), "mu"),
          rvalue(mu,
            cons_list(
              index_multi(rvalue(c,
                            cons_list(index_omni(),
                              cons_list(index_uni(2), nil_index_list())),
                            "c")), nil_index_list()), "mu")),
        "assigning variable mu_diff");
      {
        current_statement__ = 8;
        lp_accum__.add(gamma_lpdf<propto__>(sigma, alpha, beta_gamma));
        current_statement__ = 9;
        lp_accum__.add(normal_lpdf<propto__>(eta, 0, 1));
        current_statement__ = 10;
        lp_accum__.add(
          normal_lpdf<propto__>(prior_mu_not, mu_not, sigma_mu_not));
        current_statement__ = 11;
        lp_accum__.add(
          normal_lpdf<propto__>(mu, add(prior_mu_not, multiply(sigma, eta)),
            sigma));
        current_statement__ = 12;
        lp_accum__.add(
          normal_lpdf<propto__>(y, multiply(x, mu), multiply(sigma, u)));
        current_statement__ = 13;
        lp_accum__.add(normal_lpdf<propto__>(y_diff, mu_diff, sigma));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_uncertainty_model_semi_informative_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      Eigen::Matrix<double, -1, 1> mu;
      mu = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(mu, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      mu = in__.vector(K);
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      sigma = stan::math::lb_constrain(sigma, 0);
      std::vector<double> y_diff;
      y_diff = std::vector<double>(C, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        current_statement__ = 3;
        assign(y_diff, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable y_diff");}
      Eigen::Matrix<double, -1, 1> eta;
      eta = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(eta, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 4;
      eta = in__.vector(K);
      Eigen::Matrix<double, -1, 1> prior_mu_not;
      prior_mu_not = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(prior_mu_not, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      prior_mu_not = in__.vector(K);
      Eigen::Matrix<double, -1, 1> mu_diff;
      mu_diff = Eigen::Matrix<double, -1, 1>(C);
      stan::math::fill(mu_diff, std::numeric_limits<double>::quiet_NaN());
      
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(mu[(sym1__ - 1)]);}
      vars__.emplace_back(sigma);
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        vars__.emplace_back(y_diff[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(eta[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(prior_mu_not[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 7;
      assign(mu_diff, nil_index_list(),
        subtract(
          rvalue(mu,
            cons_list(
              index_multi(rvalue(c,
                            cons_list(index_omni(),
                              cons_list(index_uni(1), nil_index_list())),
                            "c")), nil_index_list()), "mu"),
          rvalue(mu,
            cons_list(
              index_multi(rvalue(c,
                            cons_list(index_omni(),
                              cons_list(index_uni(2), nil_index_list())),
                            "c")), nil_index_list()), "mu")),
        "assigning variable mu_diff");
      if (emit_transformed_parameters__) {
        for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
          vars__.emplace_back(mu_diff[(sym1__ - 1)]);}
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      Eigen::Matrix<double, -1, 1> mu;
      mu = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(mu, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> mu_flat__;
        current_statement__ = 1;
        assign(mu_flat__, nil_index_list(), context__.vals_r("mu"),
          "assigning variable mu_flat__");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 1;
          assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
            mu_flat__[(pos__ - 1)], "assigning variable mu");
          current_statement__ = 1;
          pos__ = (pos__ + 1);}
      }
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      double sigma_free__;
      sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma_free__ = stan::math::lb_free(sigma, 0);
      std::vector<double> y_diff;
      y_diff = std::vector<double>(C, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 3;
      assign(y_diff, nil_index_list(), context__.vals_r("y_diff"),
        "assigning variable y_diff");
      Eigen::Matrix<double, -1, 1> eta;
      eta = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(eta, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> eta_flat__;
        current_statement__ = 4;
        assign(eta_flat__, nil_index_list(), context__.vals_r("eta"),
          "assigning variable eta_flat__");
        current_statement__ = 4;
        pos__ = 1;
        current_statement__ = 4;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 4;
          assign(eta, cons_list(index_uni(sym1__), nil_index_list()),
            eta_flat__[(pos__ - 1)], "assigning variable eta");
          current_statement__ = 4;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> prior_mu_not;
      prior_mu_not = Eigen::Matrix<double, -1, 1>(K);
      stan::math::fill(prior_mu_not, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> prior_mu_not_flat__;
        current_statement__ = 5;
        assign(prior_mu_not_flat__, nil_index_list(),
          context__.vals_r("prior_mu_not"),
          "assigning variable prior_mu_not_flat__");
        current_statement__ = 5;
        pos__ = 1;
        current_statement__ = 5;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 5;
          assign(prior_mu_not,
            cons_list(index_uni(sym1__), nil_index_list()),
            prior_mu_not_flat__[(pos__ - 1)],
            "assigning variable prior_mu_not");
          current_statement__ = 5;
          pos__ = (pos__ + 1);}
      }
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(mu[(sym1__ - 1)]);}
      vars__.emplace_back(sigma_free__);
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        vars__.emplace_back(y_diff[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(eta[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        vars__.emplace_back(prior_mu_not[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("mu");
    names__.emplace_back("sigma");
    names__.emplace_back("y_diff");
    names__.emplace_back("eta");
    names__.emplace_back("prior_mu_not");
    names__.emplace_back("mu_diff");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(K)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(C)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(K)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(K)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(C)});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "y_diff" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "eta" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "prior_mu_not" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mu_diff" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "y_diff" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "eta" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "prior_mu_not" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mu_diff" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"y_diff\",\"type\":{\"name\":\"array\",\"length\":" << C << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"eta\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"prior_mu_not\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"mu_diff\",\"type\":{\"name\":\"vector\",\"length\":" << C << "},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"y_diff\",\"type\":{\"name\":\"array\",\"length\":" << C << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"eta\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"prior_mu_not\",\"type\":{\"name\":\"vector\",\"length\":" << K << "},\"block\":\"parameters\"},{\"name\":\"mu_diff\",\"type\":{\"name\":\"vector\",\"length\":" << C << "},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_uncertainty_model_semi_informative_namespace::model_uncertainty_model_semi_informative;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_uncertainty_model_semi_informative_namespace::profiles__;
}
#endif
#endif
