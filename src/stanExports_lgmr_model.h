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
// Code generated by %%NAME%% %%VERSION%%
#include <stan/model/model_header.hpp>
namespace model_lgmr_model_namespace {
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
                                                      " (in 'string', line 18, column 2 to column 24)",
                                                      " (in 'string', line 19, column 2 to column 27)",
                                                      " (in 'string', line 20, column 2 to column 9)",
                                                      " (in 'string', line 21, column 2 to column 11)",
                                                      " (in 'string', line 22, column 2 to column 27)",
                                                      " (in 'string', line 23, column 2 to column 36)",
                                                      " (in 'string', line 26, column 2 to column 33)",
                                                      " (in 'string', line 27, column 2 to column 36)",
                                                      " (in 'string', line 42, column 2 to column 13)",
                                                      " (in 'string', line 44, column 11 to column 12)",
                                                      " (in 'string', line 44, column 4 to column 62)",
                                                      " (in 'string', line 45, column 4 to column 12)",
                                                      " (in 'string', line 46, column 4 to column 14)",
                                                      " (in 'string', line 47, column 4 to column 27)",
                                                      " (in 'string', line 43, column 2 to line 48, column 3)",
                                                      " (in 'string', line 49, column 2 to column 22)",
                                                      " (in 'string', line 30, column 2 to column 50)",
                                                      " (in 'string', line 31, column 2 to column 33)",
                                                      " (in 'string', line 32, column 2 to column 31)",
                                                      " (in 'string', line 33, column 2 to column 42)",
                                                      " (in 'string', line 34, column 2 to column 31)",
                                                      " (in 'string', line 35, column 2 to column 31)",
                                                      " (in 'string', line 37, column 11 to column 12)",
                                                      " (in 'string', line 37, column 4 to column 68)",
                                                      " (in 'string', line 38, column 4 to column 40)",
                                                      " (in 'string', line 36, column 2 to line 39, column 3)",
                                                      " (in 'string', line 9, column 2 to column 17)",
                                                      " (in 'string', line 10, column 18 to column 19)",
                                                      " (in 'string', line 10, column 2 to column 23)",
                                                      " (in 'string', line 11, column 9 to column 10)",
                                                      " (in 'string', line 11, column 2 to column 14)",
                                                      " (in 'string', line 14, column 2 to column 25)",
                                                      " (in 'string', line 15, column 9 to column 10)",
                                                      " (in 'string', line 15, column 2 to column 41)",
                                                      " (in 'string', line 23, column 31 to column 32)",
                                                      " (in 'string', line 3, column 11 to column 12)",
                                                      " (in 'string', line 3, column 4 to column 53)",
                                                      " (in 'string', line 4, column 14 to column 44)",
                                                      " (in 'string', line 5, column 4 to column 20)",
                                                      " (in 'string', line 2, column 85 to line 6, column 3)"};
template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T4__, typename T5__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
T2__, T3__,
T4__, stan::promote_args_t<T5__>>, -1, 1>
reg_function(const T0__& x_arg__, const T1__& p_arg__, const T2__& I,
             const T3__& I_L, const T4__& S, const T5__& S_L, const int& N,
             std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          stan::value_type_t<T1__>,
          T2__,
          T3__,
          T4__, stan::promote_args_t<T5__>>;
  const auto& x = to_ref(x_arg__);
  const auto& p = to_ref(p_arg__);
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    current_statement__ = 36;
    validate_non_negative_index("exp_beta", "N", N);
    Eigen::Matrix<local_scalar_t__, -1, 1> exp_beta;
    exp_beta = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
    stan::math::fill(exp_beta, DUMMY_VAR__);
    
    current_statement__ = 37;
    assign(exp_beta, nil_index_list(),
      multiply(.001,
        stan::math::exp(elt_multiply(p, subtract(I_L, multiply(S_L, x))))),
      "assigning variable exp_beta");
    current_statement__ = 38;
    assign(exp_beta, nil_index_list(),
      add(stan::model::deep_copy(exp_beta),
        stan::math::exp(subtract(I, multiply(S, x)))),
      "assigning variable exp_beta");
    current_statement__ = 39;
    return exp_beta;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct reg_function_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T4__, typename T5__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
T2__, T3__,
T4__, stan::promote_args_t<T5__>>, -1, 1>
operator()(const T0__& x, const T1__& p, const T2__& I, const T3__& I_L,
           const T4__& S, const T5__& S_L, const int& N,
           std::ostream* pstream__)  const 
{
return reg_function(x, p, I, I_L, S, S_L, N, pstream__);
}
};
#include <stan_meta_header.hpp>
class model_lgmr_model final : public model_base_crtp<model_lgmr_model> {
private:
  int N;
  Eigen::Matrix<double, -1, 1> y;
  Eigen::Matrix<double, -1, 1> x;
  double v_y;
  Eigen::Matrix<double, -1, 1> x_star;
 
public:
  ~model_lgmr_model() { }
  
  inline std::string model_name() const final { return "model_lgmr_model"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = %%NAME%%3 %%VERSION%%", "stancflags = "};
  }
  
  
  model_lgmr_model(stan::io::var_context& context__,
                   unsigned int random_seed__ = 0,
                   std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_lgmr_model_namespace::model_lgmr_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 27;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 27;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 27;
      current_statement__ = 27;
      check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 28;
      validate_non_negative_index("y", "N", N);
      current_statement__ = 29;
      context__.validate_dims("data initialization","y","double",
          context__.to_vec(N));
      y = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(y, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> y_flat__;
        current_statement__ = 29;
        assign(y_flat__, nil_index_list(), context__.vals_r("y"),
          "assigning variable y_flat__");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 29;
          assign(y, cons_list(index_uni(sym1__), nil_index_list()),
            y_flat__[(pos__ - 1)], "assigning variable y");
          current_statement__ = 29;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 29;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 29;
        current_statement__ = 29;
        check_greater_or_equal(function__, "y[sym1__]", y[(sym1__ - 1)], 0);}
      current_statement__ = 30;
      validate_non_negative_index("x", "N", N);
      current_statement__ = 31;
      context__.validate_dims("data initialization","x","double",
          context__.to_vec(N));
      x = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(x, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> x_flat__;
        current_statement__ = 31;
        assign(x_flat__, nil_index_list(), context__.vals_r("x"),
          "assigning variable x_flat__");
        current_statement__ = 31;
        pos__ = 1;
        current_statement__ = 31;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 31;
          assign(x, cons_list(index_uni(sym1__), nil_index_list()),
            x_flat__[(pos__ - 1)], "assigning variable x");
          current_statement__ = 31;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 32;
      v_y = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 32;
      v_y = variance(y);
      current_statement__ = 33;
      validate_non_negative_index("x_star", "N", N);
      current_statement__ = 34;
      x_star = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(x_star, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 34;
      assign(x_star, nil_index_list(), divide(subtract(x, mean(x)), sd(x)),
        "assigning variable x_star");
      current_statement__ = 35;
      validate_non_negative_index("p", "N", N);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 2;
      num_params_r__ += N;
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
    static const char* function__ = "model_lgmr_model_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ alpha;
      alpha = DUMMY_VAR__;
      
      current_statement__ = 1;
      alpha = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        alpha = stan::math::lb_constrain(alpha, 0, lp__);
      } else {
        current_statement__ = 1;
        alpha = stan::math::lb_constrain(alpha, 0);
      }
      local_scalar_t__ alpha_mu;
      alpha_mu = DUMMY_VAR__;
      
      current_statement__ = 2;
      alpha_mu = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        alpha_mu = stan::math::lb_constrain(alpha_mu, 0, lp__);
      } else {
        current_statement__ = 2;
        alpha_mu = stan::math::lb_constrain(alpha_mu, 0);
      }
      local_scalar_t__ I;
      I = DUMMY_VAR__;
      
      current_statement__ = 3;
      I = in__.scalar();
      local_scalar_t__ I_L;
      I_L = DUMMY_VAR__;
      
      current_statement__ = 4;
      I_L = in__.scalar();
      Eigen::Matrix<local_scalar_t__, -1, 1> eta;
      eta = Eigen::Matrix<local_scalar_t__, -1, 1>(2);
      stan::math::fill(eta, DUMMY_VAR__);
      
      current_statement__ = 5;
      eta = in__.vector(2);
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        current_statement__ = 5;
        if (jacobian__) {
          current_statement__ = 5;
          assign(eta, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(eta[(sym1__ - 1)], 0, lp__),
            "assigning variable eta");
        } else {
          current_statement__ = 5;
          assign(eta, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(eta[(sym1__ - 1)], 0),
            "assigning variable eta");
        }}
      Eigen::Matrix<local_scalar_t__, -1, 1> p;
      p = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
      stan::math::fill(p, DUMMY_VAR__);
      
      current_statement__ = 6;
      p = in__.vector(N);
      current_statement__ = 6;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 6;
        if (jacobian__) {
          current_statement__ = 6;
          assign(p, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lub_constrain(p[(sym1__ - 1)], 0, 1, lp__),
            "assigning variable p");
        } else {
          current_statement__ = 6;
          assign(p, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lub_constrain(p[(sym1__ - 1)], 0, 1),
            "assigning variable p");
        }}
      local_scalar_t__ S;
      S = DUMMY_VAR__;
      
      current_statement__ = 7;
      S = eta[(1 - 1)];
      local_scalar_t__ S_L;
      S_L = DUMMY_VAR__;
      
      current_statement__ = 8;
      S_L = (eta[(2 - 1)] * .1);
      current_statement__ = 7;
      current_statement__ = 7;
      check_greater_or_equal(function__, "S", S, 0);
      current_statement__ = 8;
      current_statement__ = 8;
      check_greater_or_equal(function__, "S_L", S_L, 0);
      {
        current_statement__ = 17;
        lp_accum__.add(exp_mod_normal_lpdf<propto__>(alpha, alpha_mu, 1, .1));
        current_statement__ = 18;
        lp_accum__.add(normal_lpdf<propto__>(alpha_mu, 50, 10));
        current_statement__ = 19;
        lp_accum__.add(std_normal_lpdf<propto__>(I));
        current_statement__ = 20;
        lp_accum__.add(skew_normal_lpdf<propto__>(I_L, 2, 15, 35));
        current_statement__ = 21;
        lp_accum__.add(std_normal_lpdf<propto__>(eta));
        current_statement__ = 22;
        lp_accum__.add(beta_lpdf<propto__>(p, .5, .5));
        {
          current_statement__ = 23;
          validate_non_negative_index("exp_beta", "N", N);
          Eigen::Matrix<local_scalar_t__, -1, 1> exp_beta;
          exp_beta = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
          stan::math::fill(exp_beta, DUMMY_VAR__);
          
          current_statement__ = 24;
          assign(exp_beta, nil_index_list(),
            reg_function(x_star, p, I, I_L, S, S_L, N, pstream__),
            "assigning variable exp_beta");
          current_statement__ = 25;
          lp_accum__.add(
            gamma_lpdf<propto__>(y, alpha, elt_divide(alpha, exp_beta)));
        }
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
    static const char* function__ = "model_lgmr_model_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha = in__.scalar();
      current_statement__ = 1;
      alpha = stan::math::lb_constrain(alpha, 0);
      double alpha_mu;
      alpha_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      alpha_mu = in__.scalar();
      current_statement__ = 2;
      alpha_mu = stan::math::lb_constrain(alpha_mu, 0);
      double I;
      I = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      I = in__.scalar();
      double I_L;
      I_L = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      I_L = in__.scalar();
      Eigen::Matrix<double, -1, 1> eta;
      eta = Eigen::Matrix<double, -1, 1>(2);
      stan::math::fill(eta, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      eta = in__.vector(2);
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        current_statement__ = 5;
        assign(eta, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(eta[(sym1__ - 1)], 0),
          "assigning variable eta");}
      Eigen::Matrix<double, -1, 1> p;
      p = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(p, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 6;
      p = in__.vector(N);
      current_statement__ = 6;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 6;
        assign(p, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lub_constrain(p[(sym1__ - 1)], 0, 1),
          "assigning variable p");}
      double S;
      S = std::numeric_limits<double>::quiet_NaN();
      
      double S_L;
      S_L = std::numeric_limits<double>::quiet_NaN();
      
      vars__.emplace_back(alpha);
      vars__.emplace_back(alpha_mu);
      vars__.emplace_back(I);
      vars__.emplace_back(I_L);
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        vars__.emplace_back(eta[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        vars__.emplace_back(p[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 7;
      S = eta[(1 - 1)];
      current_statement__ = 8;
      S_L = (eta[(2 - 1)] * .1);
      current_statement__ = 7;
      current_statement__ = 7;
      check_greater_or_equal(function__, "S", S, 0);
      current_statement__ = 8;
      current_statement__ = 8;
      check_greater_or_equal(function__, "S_L", S_L, 0);
      if (emit_transformed_parameters__) {
        vars__.emplace_back(S);
        vars__.emplace_back(S_L);
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      double nrmse;
      nrmse = std::numeric_limits<double>::quiet_NaN();
      
      {
        current_statement__ = 10;
        validate_non_negative_index("se", "N", N);
        Eigen::Matrix<double, -1, 1> se;
        se = Eigen::Matrix<double, -1, 1>(N);
        stan::math::fill(se, std::numeric_limits<double>::quiet_NaN());
        
        current_statement__ = 11;
        assign(se, nil_index_list(),
          reg_function(x_star, p, I, I_L, S, S_L, N, pstream__),
          "assigning variable se");
        current_statement__ = 12;
        assign(se, nil_index_list(), subtract(stan::model::deep_copy(se), y),
          "assigning variable se");
        current_statement__ = 13;
        assign(se, nil_index_list(), pow(stan::model::deep_copy(se), 2),
          "assigning variable se");
        current_statement__ = 14;
        nrmse = (mean(se) / v_y);
      }
      current_statement__ = 16;
      nrmse = stan::math::sqrt(nrmse);
      vars__.emplace_back(nrmse);
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
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha = context__.vals_r("alpha")[(1 - 1)];
      double alpha_free__;
      alpha_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha_free__ = stan::math::lb_free(alpha, 0);
      double alpha_mu;
      alpha_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      alpha_mu = context__.vals_r("alpha_mu")[(1 - 1)];
      double alpha_mu_free__;
      alpha_mu_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      alpha_mu_free__ = stan::math::lb_free(alpha_mu, 0);
      double I;
      I = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      I = context__.vals_r("I")[(1 - 1)];
      double I_L;
      I_L = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      I_L = context__.vals_r("I_L")[(1 - 1)];
      Eigen::Matrix<double, -1, 1> eta;
      eta = Eigen::Matrix<double, -1, 1>(2);
      stan::math::fill(eta, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> eta_flat__;
        current_statement__ = 5;
        assign(eta_flat__, nil_index_list(), context__.vals_r("eta"),
          "assigning variable eta_flat__");
        current_statement__ = 5;
        pos__ = 1;
        current_statement__ = 5;
        for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
          current_statement__ = 5;
          assign(eta, cons_list(index_uni(sym1__), nil_index_list()),
            eta_flat__[(pos__ - 1)], "assigning variable eta");
          current_statement__ = 5;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> eta_free__;
      eta_free__ = Eigen::Matrix<double, -1, 1>(2);
      stan::math::fill(eta_free__, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        current_statement__ = 5;
        assign(eta_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(eta[(sym1__ - 1)], 0),
          "assigning variable eta_free__");}
      Eigen::Matrix<double, -1, 1> p;
      p = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(p, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> p_flat__;
        current_statement__ = 6;
        assign(p_flat__, nil_index_list(), context__.vals_r("p"),
          "assigning variable p_flat__");
        current_statement__ = 6;
        pos__ = 1;
        current_statement__ = 6;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 6;
          assign(p, cons_list(index_uni(sym1__), nil_index_list()),
            p_flat__[(pos__ - 1)], "assigning variable p");
          current_statement__ = 6;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> p_free__;
      p_free__ = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(p_free__, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 6;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 6;
        assign(p_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lub_free(p[(sym1__ - 1)], 0, 1),
          "assigning variable p_free__");}
      vars__.emplace_back(alpha_free__);
      vars__.emplace_back(alpha_mu_free__);
      vars__.emplace_back(I);
      vars__.emplace_back(I_L);
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        vars__.emplace_back(eta_free__[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        vars__.emplace_back(p_free__[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("alpha");
    names__.emplace_back("alpha_mu");
    names__.emplace_back("I");
    names__.emplace_back("I_L");
    names__.emplace_back("eta");
    names__.emplace_back("p");
    names__.emplace_back("S");
    names__.emplace_back("S_L");
    names__.emplace_back("nrmse");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(2)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(N)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "alpha");
    param_names__.emplace_back(std::string() + "alpha_mu");
    param_names__.emplace_back(std::string() + "I");
    param_names__.emplace_back(std::string() + "I_L");
    for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "eta" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "p" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "S");
      param_names__.emplace_back(std::string() + "S_L");
    }
    
    if (emit_generated_quantities__) {
      param_names__.emplace_back(std::string() + "nrmse");
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "alpha");
    param_names__.emplace_back(std::string() + "alpha_mu");
    param_names__.emplace_back(std::string() + "I");
    param_names__.emplace_back(std::string() + "I_L");
    for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "eta" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "p" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "S");
      param_names__.emplace_back(std::string() + "S_L");
    }
    
    if (emit_generated_quantities__) {
      param_names__.emplace_back(std::string() + "nrmse");
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"alpha_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"I\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"I_L\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"eta\",\"type\":{\"name\":\"vector\",\"length\":" << 2 << "},\"block\":\"parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"vector\",\"length\":" << N << "},\"block\":\"parameters\"},{\"name\":\"S\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"S_L\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"nrmse\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"alpha_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"I\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"I_L\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"eta\",\"type\":{\"name\":\"vector\",\"length\":" << 2 << "},\"block\":\"parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"vector\",\"length\":" << N << "},\"block\":\"parameters\"},{\"name\":\"S\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"S_L\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"nrmse\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]";
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
using stan_model = model_lgmr_model_namespace::model_lgmr_model;
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
  return model_lgmr_model_namespace::profiles__;
}
#endif
#endif
