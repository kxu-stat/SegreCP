#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// Compiled kernels of SegreCP. Tensors are R arrays in column-major order
// (first mode fastest); a CP component is list(lambda, U) with unit factors.
//
//   tangent_proj_cpp                     projection onto a tangent space
//   rgd_step_fused_cpp                   one RGD update of one component
//   make_rank1_vec_cpp                   dense rank-one tensor
//   update_rank1_cache_inplace_cpp       replace a component in the residual
//   initialize_decomposition_state_cpp   residual and caches at the start
//   rgd_sweep_gs_batch_cpp               one RGD iteration over all components
//   modewise_design_fibres_cpp           regression design contracted with
//                                        all factors but one

// Projection of the tensor Tvec onto the tangent space at u_1 x ... x u_d:
//   sum_k u_1 x ... x g_k x ... x u_d - (d - 1) <T, u> u,
// where g_k is T contracted with all factors but u_k and u = u_1 x ... x u_d.
// [[Rcpp::export]]
NumericVector tangent_proj_cpp(NumericVector Tvec, List Ulist) {
  IntegerVector dims = Tvec.attr("dim");
  int d = dims.size();
  if (Ulist.size() != d) stop("Length(U) must match tensor order.");

  std::vector< std::vector<double> > U(d);
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    if (uk.size() != dims[k])
      stop("U[[%d]] has wrong length (expected %d).", k+1, dims[k]);
    U[k] = std::vector<double>(uk.begin(), uk.end());
  }

  R_xlen_t N = Tvec.size();

  // Column-major strides for decoding linear indices
  std::vector<R_xlen_t> stride(d, 1);
  for (int k = 1; k < d; ++k) stride[k] = stride[k-1] * (R_xlen_t)dims[k-1];

  std::vector< std::vector<double> > g(d);
  for (int k = 0; k < d; ++k) g[k].assign(dims[k], 0.0);
  double lambda_all = 0.0;

  // Prefix/suffix buffers (reused each entry): pre[k] = prod_{l<k} u_l[i_l]
  //                                            suf[k] = prod_{l>=k} u_l[i_l]
  std::vector<double> pre(d + 1), suf(d + 1);
  std::vector<int>    idx(d);

  // ---- Pass 1: compute lambda_all and g_k with O(d) work per entry ----
  for (R_xlen_t lin = 0; lin < N; ++lin) {
    R_xlen_t tmp = lin;
    for (int k = d-1; k >= 0; --k) {
      idx[k] = (int)(tmp / stride[k]);
      tmp   -= (R_xlen_t)idx[k] * stride[k];
    }

    // Build prefix and suffix products
    pre[0] = 1.0;
    for (int k = 0; k < d; ++k) pre[k+1] = pre[k] * U[k][idx[k]];
    suf[d] = 1.0;
    for (int k = d-1; k >= 0; --k) suf[k] = suf[k+1] * U[k][idx[k]];

    double Tij = Tvec[lin];
    lambda_all += Tij * pre[d];           // pre[d] = prod_all

    // g_k[i_k] += Tij * pre[k] * suf[k+1]  (= Tij * prod_{l != k} u_l[i_l])
    for (int k = 0; k < d; ++k)
      g[k][idx[k]] += Tij * pre[k] * suf[k+1];
  }

  // ---- Pass 2: assemble projected tensor with O(d) work per entry ----
  NumericVector out(N);
  out.attr("dim") = clone(dims);

  for (R_xlen_t lin = 0; lin < N; ++lin) {
    R_xlen_t tmp = lin;
    for (int k = d-1; k >= 0; --k) {
      idx[k] = (int)(tmp / stride[k]);
      tmp   -= (R_xlen_t)idx[k] * stride[k];
    }

    pre[0] = 1.0;
    for (int k = 0; k < d; ++k) pre[k+1] = pre[k] * U[k][idx[k]];
    suf[d] = 1.0;
    for (int k = d-1; k >= 0; --k) suf[k] = suf[k+1] * U[k][idx[k]];

    double val = 0.0;
    for (int k = 0; k < d; ++k)
      val += g[k][idx[k]] * pre[k] * suf[k+1];

    // Subtract (d-1)*lambda_all*prod_all (tangent-space correction)
    val -= (double)(d - 1) * lambda_all * pre[d];
    out[lin] = val;
  }

  return out;
}


// Rank-(1, ..., 1) T-HOSVD of the RGD proposal
//   Z = lambda_i u - alpha * P(R) = a u + sum_k u_1 x ... x w_k x ... x u_d,
// with u = u_1 x ... x u_d, P the tangent projection at u, g_k the
// contractions of the residual R with all factors but u_k, lambda_full =
// <R, u>, and w_k = -alpha (g_k - lambda_full u_k), which is orthogonal to
// u_k. The mode-k unfolding of Z has its column space in span(u_k, w_k), so
// its leading left singular vector comes from a 2 x 2 eigenproblem, solved
// here in closed form.
static void exact_tangent_hosvd(
    const std::vector< std::vector<double> >& U0,
    const std::vector< std::vector<double> >& g,
    double lambda_i, double alpha, double lambda_full,
    std::vector< std::vector<double> >& Unew,
    double& lambda_new) {
  const int d = U0.size();
  std::vector< std::vector<double> > w(d);
  std::vector<double> b(d, 0.0), ccoef(d, 1.0), scoef(d, 0.0);

  // In exact arithmetic g_k - lambda_full u_k is perpendicular to u_k.
  // Re-orthogonalise to remove contraction roundoff and transfer the tiny
  // parallel part back into the shared scale coefficient, preserving Z.
  double a = lambda_i - alpha * lambda_full;
  for (int k = 0; k < d; ++k) {
    const int pk = U0[k].size();
    w[k].resize(pk);
    double parallel = 0.0;
    for (int j = 0; j < pk; ++j) {
      w[k][j] = -alpha * (g[k][j] - lambda_full * U0[k][j]);
      parallel += U0[k][j] * w[k][j];
    }
    a += parallel;
    double w_nrm = 0.0;
    for (int j = 0; j < pk; ++j) {
      w[k][j] -= parallel * U0[k][j];
      w_nrm = std::hypot(w_nrm, w[k][j]);
    }
    b[k] = w_nrm;
  }

  double total_b2 = 0.0;
  for (int k = 0; k < d; ++k) total_b2 += b[k] * b[k];
  Unew.resize(d);

  for (int k = 0; k < d; ++k) {
    const double A = a * a + total_b2 - b[k] * b[k];
    const double B = a * b[k];
    const double C = b[k] * b[k];
    const double theta = 0.5 * std::atan2(2.0 * B, A - C);
    ccoef[k] = std::cos(theta);
    scoef[k] = std::sin(theta);

    const int pk = U0[k].size();
    Unew[k].resize(pk);
    double u_nrm = 0.0;
    if (b[k] > 0.0) {
      for (int j = 0; j < pk; ++j) {
        Unew[k][j] = ccoef[k] * U0[k][j] + scoef[k] * w[k][j] / b[k];
        u_nrm = std::hypot(u_nrm, Unew[k][j]);
      }
    } else {
      Unew[k] = U0[k];
      ccoef[k] = 1.0;
      scoef[k] = 0.0;
      u_nrm = 1.0;
    }
    if (u_nrm > 0.0) {
      for (int j = 0; j < pk; ++j) Unew[k][j] /= u_nrm;
      ccoef[k] /= u_nrm;
      scoef[k] /= u_nrm;
    }
  }

  std::vector<double> pre(d + 1, 1.0), suf(d + 1, 1.0);
  for (int k = 0; k < d; ++k) pre[k + 1] = pre[k] * ccoef[k];
  for (int k = d - 1; k >= 0; --k) suf[k] = suf[k + 1] * ccoef[k];
  lambda_new = a * pre[d];
  for (int k = 0; k < d; ++k)
    lambda_new += b[k] * scoef[k] * pre[k] * suf[k + 1];
}

// One RGD update of the component lambda_i u_1 x ... x u_d, with residual
// Rvec = (sum of components) - observed tensor and step alpha: one pass over
// Rvec for the contractions, then exact_tangent_hosvd().
// [[Rcpp::export]]
List rgd_step_fused_cpp(NumericVector Rvec, List Ulist,
                        double lambda_i, double alpha) {
  IntegerVector dims = Rvec.attr("dim");
  int d = dims.size();
  R_xlen_t N = Rvec.size();

  // Factors and column-major strides.
  std::vector< std::vector<double> > U0(d), Ucur(d);
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    U0[k]   = std::vector<double>(uk.begin(), uk.end());
    Ucur[k] = U0[k];
  }
  std::vector<R_xlen_t> stride(d, 1);
  for (int k = 1; k < d; ++k)
    stride[k] = stride[k-1] * (R_xlen_t)dims[k-1];

  // Contractions, in one pass over the residual:
  //   g_k[i_k]    = sum over all indices but i_k of R[i] prod_{l != k} u_l[i_l]
  //   lambda_full = <R, u_1 x ... x u_d>.
  // Orders 3 and 4 have explicit loops. The general case keeps the
  // multi-index in an odometer idx[] that is advanced with carries.
  std::vector< std::vector<double> > g(d);
  for (int k = 0; k < d; ++k) g[k].assign(dims[k], 0.0);
  double lambda_full = 0.0;

  if (d == 3) {
    // The contraction of each mode-0 fibre with u_0 is reused for g[1], g[2]
    // and lambda_full.
    R_xlen_t lin = 0;
    for (int i2 = 0; i2 < dims[2]; ++i2) {
      const double u2 = U0[2][i2];
      for (int i1 = 0; i1 < dims[1]; ++i1) {
        const double u1  = U0[1][i1];
        const double u21 = u2 * u1;
        double contracted0 = 0.0;
        for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
          const double Rij = Rvec[lin];
          g[0][i0] += Rij * u21;
          contracted0 += Rij * U0[0][i0];
        }
        g[1][i1] += contracted0 * u2;
        g[2][i2] += contracted0 * u1;
        lambda_full += contracted0 * u21;
      }
    }
  } else if (d == 4) {
    // As for order 3: the mode-0 fibre contraction is reused for modes 1 to 3
    // and lambda_full.
    R_xlen_t lin = 0;
    for (int i3 = 0; i3 < dims[3]; ++i3) {
      const double u3 = U0[3][i3];
      for (int i2 = 0; i2 < dims[2]; ++i2) {
        const double u2  = U0[2][i2];
        const double u32 = u3 * u2;
        for (int i1 = 0; i1 < dims[1]; ++i1) {
          const double u1   = U0[1][i1];
          const double u321 = u32 * u1;
          double contracted0 = 0.0;
          for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
            const double Rij = Rvec[lin];
            g[0][i0] += Rij * u321;
            contracted0 += Rij * U0[0][i0];
          }
          g[1][i1] += contracted0 * u32;
          g[2][i2] += contracted0 * u1 * u3;
          g[3][i3] += contracted0 * u1 * u2;
          lambda_full += contracted0 * u321;
        }
      }
    }
  } else {
    std::vector<double> pre(d+1), suf(d+1);
    std::vector<int>    idx(d, 0);

    for (R_xlen_t lin = 0; lin < N; ++lin) {
      pre[0] = 1.0;
      for (int k = 0; k < d; ++k)
        pre[k+1] = pre[k] * U0[k][idx[k]];
      suf[d] = 1.0;
      for (int k = d-1; k >= 0; --k)
        suf[k] = suf[k+1] * U0[k][idx[k]];

      const double Rij = Rvec[lin];
      lambda_full += Rij * pre[d];
      for (int k = 0; k < d; ++k)
        g[k][idx[k]] += Rij * pre[k] * suf[k+1];

      if (lin < N - 1) {
        int carry = 0;
        idx[carry]++;
        while (carry < d - 1 && idx[carry] == dims[carry]) {
          idx[carry] = 0;
          carry++;
          idx[carry]++;
        }
      }
    }
  }

  // Retraction.
  double lambda_new = 0.0;
  exact_tangent_hosvd(
    U0, g, lambda_i, alpha, lambda_full, Ucur, lambda_new
  );

  List U_out(d);
  for (int k = 0; k < d; ++k) {
    NumericVector uk(dims[k]);
    for (int i = 0; i < dims[k]; ++i) uk[i] = Ucur[k][i];
    U_out[k] = uk;
  }
  return List::create(Named("U") = U_out, Named("lambda") = lambda_new);
}

// Dense tensor lambda * u_1 x ... x u_d. Each entry is computed as
// (((lambda * u_1[i_1]) * u_2[i_2]) * ...) * u_d[i_d].
// [[Rcpp::export]]
NumericVector make_rank1_vec_cpp(double lambda, List Ulist) {
  int d = Ulist.size();
  std::vector< std::vector<double> > U(d);
  std::vector<int> dims(d);
  R_xlen_t N = 1;
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    U[k]    = std::vector<double>(uk.begin(), uk.end());
    dims[k] = (int)uk.size();
    N      *= dims[k];
  }

  NumericVector out(N);
  std::vector<int> idx(d, 0);

  for (R_xlen_t lin = 0; lin < N; ++lin) {
    double val = lambda;
    for (int k = 0; k < d; ++k) val *= U[k][idx[k]];
    out[lin] = val;

    // Advance odometer
    if (lin < N - 1) {
      int carry = 0;
      idx[carry]++;
      while (carry < d - 1 && idx[carry] == dims[carry]) {
        idx[carry] = 0;
        carry++;
        idx[carry]++;
      }
    }
  }

  // The dim attribute makes the result an R array without a copy.
  IntegerVector rdims(d);
  for (int k = 0; k < d; ++k) rdims[k] = dims[k];
  out.attr("dim") = rdims;
  return out;
}

// Replaces a component in place: residual += new - cache, cache = new, where
// new = lambda * u_1 x ... x u_d.
// [[Rcpp::export]]
void update_rank1_cache_inplace_cpp(NumericVector residual,
                                    NumericVector component_cache,
                                    double lambda, List Ulist) {
  int d = Ulist.size();
  std::vector< std::vector<double> > U(d);
  std::vector<int> dims(d);
  R_xlen_t N = 1;
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    U[k]    = std::vector<double>(uk.begin(), uk.end());
    dims[k] = (int)uk.size();
    N      *= dims[k];
  }
  if (residual.size() != N || component_cache.size() != N)
    stop("Residual/cache length does not match the factor dimensions.");

  std::vector<int> idx(d, 0);
  for (R_xlen_t lin = 0; lin < N; ++lin) {
    double value = lambda;
    for (int k = 0; k < d; ++k) value *= U[k][idx[k]];

    residual[lin] += value - component_cache[lin];
    component_cache[lin] = value;

    if (lin < N - 1) {
      int carry = 0;
      idx[carry]++;
      while (carry < d - 1 && idx[carry] == dims[carry]) {
        idx[carry] = 0;
        carry++;
        idx[carry]++;
      }
    }
  }
}

// Dense components (caches), residual = sum of components - observed, and the
// objective 0.5 ||residual||^2 at the start, in one pass. Entries are
// multiplied as in make_rank1_vec_cpp and summed over components in order.
// [[Rcpp::export]]
List initialize_decomposition_state_cpp(NumericVector observed,
                                        List components) {
  IntegerVector dims = observed.attr("dim");
  const int d = dims.size(), r = components.size();
  const R_xlen_t N = observed.size();
  if (d < 2 || r < 1) stop("Initialization requires tensor modes and components.");
  R_xlen_t expected = 1;
  for (int k = 0; k < d; ++k) {
    if (dims[k] < 1 || expected > N / dims[k])
      stop("Observed tensor dimensions do not match its length.");
    expected *= dims[k];
  }
  if (expected != N) stop("Observed tensor dimensions do not match its length.");
  typedef std::vector< std::vector<double> > FactorVectors;
  std::vector<FactorVectors> U(r, FactorVectors(d));
  std::vector<double> lambda(r);
  List caches(r);
  std::vector<double*> cache_ptr(r);
  for (int i = 0; i < r; ++i) {
    List component = components[i];
    lambda[i] = as<double>(component["lambda"]);
    if (!std::isfinite(lambda[i]) || lambda[i] <= 0.0)
      stop("Components must have positive finite weights.");
    List factors = component["U"];
    if (factors.size() != d) stop("Component order does not match observed tensor.");
    for (int k = 0; k < d; ++k) {
      NumericVector u = factors[k];
      if (u.size() != dims[k]) stop("Component factors have wrong dimensions.");
      U[i][k] = std::vector<double>(u.begin(), u.end());
    }
    NumericVector component_cache(N);
    caches[i] = component_cache;
    cache_ptr[i] = component_cache.begin();
  }
  NumericVector residual(N);
  residual.attr("dim") = clone(dims);
  const double* observed_ptr = observed.begin();
  double* residual_ptr = residual.begin();
  std::vector<int> idx(d, 0);
  long double residual_sq = 0.0L;
  for (R_xlen_t lin = 0; lin < N; ++lin) {
    double estimate = 0.0;
    for (int i = 0; i < r; ++i) {
      double value = lambda[i];
      for (int k = 0; k < d; ++k) value *= U[i][k][idx[k]];
      cache_ptr[i][lin] = value;
      if (i == 0) estimate = value; else estimate += value;
    }
    const double value = estimate - observed_ptr[lin];
    residual_ptr[lin] = value;
    residual_sq += static_cast<long double>(value * value);
    int k = 0;
    while (k < d && ++idx[k] == dims[k]) { idx[k] = 0; ++k; }
  }
  return List::create(Named("residual") = residual, Named("caches") = caches,
    Named("objective") = static_cast<double>(0.5L * residual_sq));
}

// One RGD iteration (step alpha; RGN for alpha = 1) over all components,
// updated one at a time in the order order_i (one-based), each with the
// residual that includes the components already updated in this iteration.
//
// 1. One pass over the residual computes the contractions of every
//    component (as in rgd_step_fused_cpp).
// 2. Before component i is updated, the change of each earlier component j,
//    from T_j to T_j', is added to these contractions exactly:
//    <T_j' - T_j, ...> needs only inner products of factors, O(r^2 sum p_k)
//    operations in total.
// 3. One pass updates the residual and the caches.
//
// The result equals updating the components one at a time with
// rgd_step_fused_cpp up to rounding. When the residual is so small that the
// corrections in step 2 could lose accuracy, the components are updated one
// at a time instead. previous_objective, if given, is 0.5 ||residual||^2 and
// saves one pass. residual and cache_list are modified in place; the returned
// factors are new vectors.
// [[Rcpp::export]]
List rgd_sweep_gs_batch_cpp(NumericVector residual, List cache_list,
                           List components, IntegerVector order_i,
                           double alpha,
                           Nullable<NumericVector> previous_objective = R_NilValue) {
  IntegerVector dims = residual.attr("dim");
  const int d = dims.size(), r = components.size();
  const R_xlen_t N = residual.size();
  if (d < 2 || r < 1 || cache_list.size() != r || order_i.size() != r)
    stop("Invalid dimensions, caches, or update order for a GS sweep.");
  if (!std::isfinite(alpha) || alpha <= 0.0)
    stop("The GS step size must be positive and finite.");
  R_xlen_t expected = 1;
  for (int k = 0; k < d; ++k) {
    if (dims[k] < 1) stop("Tensor mode dimensions must be positive.");
    expected *= dims[k];
  }
  if (expected != N) stop("Residual length does not match its dimensions.");
  std::vector<int> order(r), seen(r, 0);
  for (int s = 0; s < r; ++s) {
    if (order_i[s] == NA_INTEGER || order_i[s] < 1 || order_i[s] > r)
      stop("order_i must be a one-based component permutation.");
    const int i = order_i[s] - 1;
    if (seen[i])
      stop("order_i must be a one-based component permutation.");
    order[s] = i; seen[i] = 1;
  }

  typedef std::vector< std::vector<double> > FactorVectors;
  std::vector<FactorVectors> U(r, FactorVectors(d)), g(r, FactorVectors(d));
  std::vector<double> lambda(r), full(r, 0.0);
  std::vector<NumericVector> cache;
  cache.reserve(r);
  double weight_sum = 0.0;
  for (int i = 0; i < r; ++i) {
    List comp = components[i];
    lambda[i] = as<double>(comp["lambda"]);
    if (!std::isfinite(lambda[i]) || lambda[i] <= 0.0)
      stop("GS components must have positive finite weights.");
    weight_sum += lambda[i];
    List factors = comp["U"];
    if (factors.size() != d) stop("Component order does not match residual.");
    for (int k = 0; k < d; ++k) {
      NumericVector u = factors[k];
      if (u.size() != dims[k]) stop("Component factors have wrong dimensions.");
      U[i][k] = std::vector<double>(u.begin(), u.end());
      g[i][k].assign(dims[k], 0.0);
    }
    NumericVector ci = cache_list[i];
    if (ci.size() != N) stop("Component cache has wrong length.");
    cache.push_back(ci);
  }

  double* residual_ptr = residual.begin();
  // ||residual||^2, from previous_objective if given.
  long double initial_residual_sq = 0.0L;
  if (previous_objective.isNotNull()) {
    NumericVector previous(previous_objective.get());
    if (previous.size() != 1 || !std::isfinite(previous[0]) || previous[0] < 0.0)
      stop("previous_objective must be one finite nonnegative scalar.");
    initial_residual_sq = 2.0L * previous[0];
  } else {
    for (R_xlen_t j = 0; j < N; ++j)
      initial_residual_sq += static_cast<long double>(residual_ptr[j] * residual_ptr[j]);
  }
  // Small residual: update the components one at a time (see above).
  const double precision_scale = 64.0 * std::sqrt(std::numeric_limits<double>::epsilon()) * weight_sum;
  long double residual_norm = std::sqrt(initial_residual_sq);
  if (previous_objective.isNotNull() &&
      std::abs(residual_norm - precision_scale) <=
        8.0L * std::numeric_limits<double>::epsilon() * precision_scale) {
    // At the threshold, recompute the norm from the residual so that the
    // choice does not depend on how previous_objective was rounded.
    initial_residual_sq = 0.0L;
    for (R_xlen_t j = 0; j < N; ++j)
      initial_residual_sq += static_cast<long double>(residual_ptr[j] * residual_ptr[j]);
    residual_norm = std::sqrt(initial_residual_sq);
  }
  if (residual_norm <= precision_scale) {
    List out = clone(components);
    for (int s = 0; s < r; ++s) {
      const int i = order[s];
      List old = out[i];
      List proposal = rgd_step_fused_cpp(residual, as<List>(old["U"]),
        as<double>(old["lambda"]), alpha);
      double next_lambda = as<double>(proposal["lambda"]);
      List next_U = proposal["U"];
      if (next_lambda < 0.0) {
        next_lambda = -next_lambda;
        NumericVector first = clone(as<NumericVector>(next_U[0]));
        for (R_xlen_t j = 0; j < first.size(); ++j) first[j] = -first[j];
        next_U[0] = first;
      }
      if (!std::isfinite(next_lambda) || next_lambda <= std::numeric_limits<double>::epsilon())
        stop("Component %d retracted to a zero/nonfinite tensor.", i + 1);
      update_rank1_cache_inplace_cpp(residual, cache[i], next_lambda, next_U);
      out[i] = List::create(Named("lambda") = next_lambda, Named("U") = next_U);
    }
    long double residual_sq = 0.0L;
    for (R_xlen_t j = 0; j < N; ++j)
      residual_sq += static_cast<long double>(residual_ptr[j] * residual_ptr[j]);
    return List::create(Named("components") = out,
      Named("objective") = static_cast<double>(0.5L * residual_sq),
      Named("sequential_fallback") = true);
  }

  std::vector<double*> g0(r);
  std::vector<const double*> u0(r);
  for (int i = 0; i < r; ++i) {
    g0[i] = g[i][0].data();
    u0[i] = U[i][0].data();
  }
  R_xlen_t lin = 0;
  if (d == 4 && r == 3) {
    // Step 1 for rank 3 and order 4, unrolled over the components.
    for (int i3 = 0; i3 < dims[3]; ++i3)
      for (int i2 = 0; i2 < dims[2]; ++i2)
        for (int i1 = 0; i1 < dims[1]; ++i1) {
          double c0 = 0.0, c1 = 0.0, c2 = 0.0;
          const double p0 = U[0][1][i1] * U[0][2][i2] * U[0][3][i3];
          const double p1 = U[1][1][i1] * U[1][2][i2] * U[1][3][i3];
          const double p2 = U[2][1][i1] * U[2][2][i2] * U[2][3][i3];
          for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
            const double value = residual_ptr[lin];
            g0[0][i0] += value * p0; c0 += value * u0[0][i0];
            g0[1][i0] += value * p1; c1 += value * u0[1][i0];
            g0[2][i0] += value * p2; c2 += value * u0[2][i0];
          }
          const double contracted[3] = {c0, c1, c2};
          for (int i = 0; i < 3; ++i) {
            const double u1 = U[i][1][i1], u2 = U[i][2][i2], u3 = U[i][3][i3];
            g[i][1][i1] += contracted[i] * u2 * u3;
            g[i][2][i2] += contracted[i] * u1 * u3;
            g[i][3][i3] += contracted[i] * u1 * u2;
            full[i] += contracted[i] * u1 * u2 * u3;
          }
        }
  } else {
    // Step 1, general order and rank: the contraction of each mode-0 fibre
    // is reused for the other modes.
    std::vector<int> idx(d, 0);
    std::vector<double> product(r), contracted(r), pre(d + 1), suf(d + 1);
    const R_xlen_t fibres = N / dims[0];
    for (R_xlen_t f = 0; f < fibres; ++f) {
      for (int i = 0; i < r; ++i) {
        product[i] = 1.0; contracted[i] = 0.0;
        for (int k = 1; k < d; ++k) product[i] *= U[i][k][idx[k]];
      }
      for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
        const double value = residual_ptr[lin];
        for (int i = 0; i < r; ++i) {
          g0[i][i0] += value * product[i];
          contracted[i] += value * u0[i][i0];
        }
      }
      for (int i = 0; i < r; ++i) {
        pre[1] = 1.0; suf[d] = 1.0;
        for (int k = 1; k < d; ++k) pre[k + 1] = pre[k] * U[i][k][idx[k]];
        for (int k = d - 1; k >= 1; --k) suf[k] = suf[k + 1] * U[i][k][idx[k]];
        full[i] += contracted[i] * product[i];
        for (int k = 1; k < d; ++k)
          g[i][k][idx[k]] += contracted[i] * pre[k] * suf[k + 1];
      }
      int k = 1;
      while (k < d && ++idx[k] == dims[k]) { idx[k] = 0; ++k; }
    }
  }

  // Step 2 and the retractions, in the order order_i.
  std::vector<FactorVectors> updated = U;
  std::vector<double> next_lambda = lambda;
  std::vector<long double> pre_old(d + 1), pre_new(d + 1),
    suf_old(d + 1), suf_new(d + 1), dot_old(d), dot_new(d);
  for (int s = 0; s < r; ++s) {
    const int i = order[s];
    for (int previous = 0; previous < s; ++previous) {
      const int j = order[previous];
      for (int k = 0; k < d; ++k) {
        dot_old[k] = 0.0L; dot_new[k] = 0.0L;
        for (int a = 0; a < dims[k]; ++a) {
          dot_old[k] += static_cast<long double>(U[j][k][a]) * U[i][k][a];
          dot_new[k] += static_cast<long double>(updated[j][k][a]) * U[i][k][a];
        }
      }
      pre_old[0] = pre_new[0] = suf_old[d] = suf_new[d] = 1.0L;
      for (int k = 0; k < d; ++k) {
        pre_old[k + 1] = pre_old[k] * dot_old[k];
        pre_new[k + 1] = pre_new[k] * dot_new[k];
      }
      for (int k = d - 1; k >= 0; --k) {
        suf_old[k] = suf_old[k + 1] * dot_old[k];
        suf_new[k] = suf_new[k + 1] * dot_new[k];
      }
      full[i] += static_cast<double>(next_lambda[j] * pre_new[d] - lambda[j] * pre_old[d]);
      for (int k = 0; k < d; ++k) {
        const long double old_scale = lambda[j] * pre_old[k] * suf_old[k + 1];
        const long double new_scale = next_lambda[j] * pre_new[k] * suf_new[k + 1];
        for (int a = 0; a < dims[k]; ++a)
          g[i][k][a] += static_cast<double>(new_scale * updated[j][k][a] - old_scale * U[j][k][a]);
      }
    }
    exact_tangent_hosvd(U[i], g[i], lambda[i], alpha, full[i], updated[i], next_lambda[i]);
    if (next_lambda[i] < 0.0) {
      next_lambda[i] = -next_lambda[i];
      for (int a = 0; a < dims[0]; ++a) updated[i][0][a] = -updated[i][0][a];
    }
    if (!std::isfinite(next_lambda[i]) || next_lambda[i] <= std::numeric_limits<double>::epsilon())
      stop("Component %d retracted to a zero/nonfinite tensor.", i + 1);
  }

  // Step 3: residual and caches in one pass, adding the component changes in
  // the order order_i.
  std::vector<int> idx(d, 0);
  std::vector<double> scale(r);
  std::vector<double*> ordered_cache(r);
  std::vector<const double*> ordered_u0(r);
  for (int s = 0; s < r; ++s) {
    const int i = order[s];
    ordered_cache[s] = cache[i].begin();
    ordered_u0[s] = updated[i][0].data();
  }
  long double residual_sq = 0.0L;
  lin = 0;
  for (R_xlen_t f = 0; f < N / dims[0]; ++f) {
    for (int s = 0; s < r; ++s) {
      const int i = order[s];
      scale[s] = next_lambda[i];
      for (int k = 1; k < d; ++k) scale[s] *= updated[i][k][idx[k]];
    }
    if (r == 3) {
      double* cache0 = ordered_cache[0];
      double* cache1 = ordered_cache[1];
      double* cache2 = ordered_cache[2];
      const double* factor0 = ordered_u0[0];
      const double* factor1 = ordered_u0[1];
      const double* factor2 = ordered_u0[2];
      const double scale0 = scale[0], scale1 = scale[1], scale2 = scale[2];
      for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
        const double value0 = scale0 * factor0[i0];
        const double value1 = scale1 * factor1[i0];
        const double value2 = scale2 * factor2[i0];
        double next_residual = residual_ptr[lin];
        next_residual += value0 - cache0[lin];
        next_residual += value1 - cache1[lin];
        next_residual += value2 - cache2[lin];
        cache0[lin] = value0; cache1[lin] = value1; cache2[lin] = value2;
        residual_ptr[lin] = next_residual;
        residual_sq += static_cast<long double>(next_residual * next_residual);
      }
    } else {
      // General rank.
      for (int i0 = 0; i0 < dims[0]; ++i0, ++lin) {
        double next_residual = residual_ptr[lin];
        for (int s = 0; s < r; ++s) {
          const double value = scale[s] * ordered_u0[s][i0];
          next_residual += value - ordered_cache[s][lin];
          ordered_cache[s][lin] = value;
        }
        residual_ptr[lin] = next_residual;
        residual_sq += static_cast<long double>(next_residual * next_residual);
      }
    }
    int k = 1;
    while (k < d && ++idx[k] == dims[k]) { idx[k] = 0; ++k; }
  }
  List out(r);
  for (int i = 0; i < r; ++i) {
    List factors(d);
    for (int k = 0; k < d; ++k) factors[k] = wrap(updated[i][k]);
    out[i] = List::create(Named("lambda") = next_lambda[i], Named("U") = factors);
  }
  return List::create(Named("components") = out,
    Named("objective") = static_cast<double>(0.5L * residual_sq),
    Named("sequential_fallback") = false);
}


// For a design X (n x prod(p_k), rows vec(X_m)) and factors u_1, ..., u_d,
// returns d matrices: column j of matrix k holds, for each observation, X_m
// contracted with all factors but u_k, at index j of mode k.
// [[Rcpp::export]]
List modewise_design_fibres_cpp(NumericMatrix X, List Ulist) {
  const int n = X.nrow();
  const int d = Ulist.size();
  if (d < 1) stop("Ulist must contain at least one mode vector.");

  std::vector<int> dims(d);
  std::vector< std::vector<double> > U(d);
  R_xlen_t p_tot = 1;
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    dims[k] = uk.size();
    if (dims[k] < 1) stop("Every mode vector must be non-empty.");
    U[k] = std::vector<double>(uk.begin(), uk.end());
    p_tot *= (R_xlen_t)dims[k];
  }
  if ((R_xlen_t)X.ncol() != p_tot)
    stop("ncol(X) does not equal prod(lengths(Ulist)).");

  std::vector<R_xlen_t> stride(d, 1);
  for (int k = 1; k < d; ++k)
    stride[k] = stride[k - 1] * (R_xlen_t)dims[k - 1];

  std::vector<NumericMatrix> fibres;
  fibres.reserve(d);
  for (int k = 0; k < d; ++k) fibres.emplace_back(n, dims[k]);

  std::vector<int> idx(d, 0);
  std::vector<double> uval(d, 0.0), pre(d + 1, 1.0), suf(d + 1, 1.0);
  const double* xptr = REAL(X);

  for (R_xlen_t lin = 0; lin < p_tot; ++lin) {
    R_xlen_t tmp = lin;
    for (int k = d - 1; k >= 0; --k) {
      idx[k] = (int)(tmp / stride[k]);
      tmp -= (R_xlen_t)idx[k] * stride[k];
      uval[k] = U[k][idx[k]];
    }
    pre[0] = 1.0;
    for (int k = 0; k < d; ++k) pre[k + 1] = pre[k] * uval[k];
    suf[d] = 1.0;
    for (int k = d - 1; k >= 0; --k) suf[k] = uval[k] * suf[k + 1];

    const double* xcol = xptr + (R_xlen_t)n * lin;
    for (int k = 0; k < d; ++k) {
      const double weight = pre[k] * suf[k + 1];
      double* out = REAL(fibres[k]) + (R_xlen_t)n * idx[k];
      for (int row = 0; row < n; ++row) out[row] += xcol[row] * weight;
    }
  }

  List ans(d);
  for (int k = 0; k < d; ++k) ans[k] = fibres[k];
  return ans;
}
