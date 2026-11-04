#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// The kernels below are extracted from the active decomposition and regression
// implementations used by the manuscript experiments. They operate on R's
// column-major array layout.

// [[Rcpp::export]]
NumericVector tangent_proj_cpp(NumericVector Tvec, List Ulist) {
  // Tvec : tensor entries as a numeric vector with dim attribute
  // Ulist: list of mode vectors u_k (each length dims[k])

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


// Exact rank-1 T-HOSVD of a rank-1 tensor plus a tangent tensor.
// The proposal has the orthogonal representation
//   Z = a U + sum_k w_k tensor U_{-k},  with w_k perpendicular to u_k.
// Every mode Gramian is therefore supported on span(u_k, w_k), so its leading
// eigenvector is obtained from a 2 x 2 symmetric eigendecomposition.
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

// ============================================================
// Kernel 1 — Fused RGD step (stride-counter, no divisions)
// ============================================================
// [[Rcpp::export]]
List rgd_step_fused_cpp(NumericVector Rvec, List Ulist,
                        double lambda_i, double alpha, int n_iter) {
  // n_iter is retained in the signature for API compatibility but unused.
  IntegerVector dims = Rvec.attr("dim");
  int d = dims.size();
  R_xlen_t N = Rvec.size();

  // ----- Extract U0 and column-major strides ----------------------------------
  std::vector< std::vector<double> > U0(d), Ucur(d);
  for (int k = 0; k < d; ++k) {
    NumericVector uk = Ulist[k];
    U0[k]   = std::vector<double>(uk.begin(), uk.end());
    Ucur[k] = U0[k];
  }
  std::vector<R_xlen_t> stride(d, 1);
  for (int k = 1; k < d; ++k)
    stride[k] = stride[k-1] * (R_xlen_t)dims[k-1];

  // ----- Pass 1: g_k and lambda_full via prefix-suffix + odometer counter ----
  //   g_k[i_k] = Σ_{l≠k} R[i] · ∏_{l≠k} u_l^{(0)}[i_l]
  //   lambda_full = ⟨R_arr, ⊗_l u_l^{(0)}⟩
  //
  //   Index decode: instead of d integer divisions per entry we maintain an
  //   "odometer" idx[] that is incremented by carries — O(1) amortised per
  //   entry vs O(d) divisions in the original implementation.
  std::vector< std::vector<double> > g(d);
  for (int k = 0; k < d; ++k) g[k].assign(dims[k], 0.0);
  double lambda_full = 0.0;

  if (d == 3) {
    // Explicit loops avoid dynamic indexing for three-mode tensors.
    // Contract mode 0 once inside each (i1,i2) fibre and reuse that scalar for
    // g[1], g[2], and lambda_full.  Compared with updating all three g arrays
    // at every tensor entry, this removes two N-length accumulation streams
    // while evaluating the same contractions.
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
    // Fourth-order counterpart of the staged third-order contraction: the
    // mode-0 fibre contraction is reused for modes 1--3 and lambda_full.
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

  // ----- Exact structure-exploiting T-HOSVD retraction -----------------------
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

// ============================================================
// Kernel 2 — Rank-1 array formation (odometer, no outer())
// ============================================================
// [[Rcpp::export]]
NumericVector make_rank1_vec_cpp(double lambda, List Ulist) {
  // Ulist: list of d vectors; result = lambda * u_1 ⊗ u_2 ⊗ ... ⊗ u_d
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

  // Attach dim attribute so the R wrapper can set it without a copy
  IntegerVector rdims(d);
  for (int k = 0; k < d; ++k) rdims[k] = dims[k];
  out.attr("dim") = rdims;
  return out;
}

// ============================================================
// Kernel 4 — Fused rank-1 cache replacement and residual update
// ============================================================
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

// Build all initial rank-one caches, their residual, and its objective in one
// tensor traversal.  The caller supplies validated observed data and unit
// factors.  Multiplication within each rank-one entry and addition across
// components retain the order of make_rank1_vec_cpp followed by Reduce(`+`).
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

// One exact-arithmetic Gauss--Seidel sweep.  The initial residual is streamed
// once for all component contractions.  Before component i is updated, add
// the contractions of every preceding new-minus-old rank-one component.
// Component i's factors have not yet changed, so these are exactly the
// contractions of its current sequential residual.  No Jacobi substitution
// or approximation to the existing T-HOSVD retraction is made.
//
// residual = sum_i T_i - observed; order_i is a one-based permutation fixed
// for this sweep.  Only residual and cache_list are mutated.  All factor
// vectors returned below own fresh storage, preserving stored CP histories.
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
  // The public wrapper already holds 0.5 * ||residual||^2 from initialization
  // or the preceding sweep.  Reuse it; the optional scan preserves internal
  // five-argument calls without recomputing norms in the contraction loop.
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
  // Numerical fallback, not a stopping or acceptance condition.  It is
  // selected before computing any batch contractions when the residual is
  // small enough that new-minus-old low-rank corrections could cancel.
  const double precision_scale = 64.0 * std::sqrt(std::numeric_limits<double>::epsilon()) * weight_sum;
  long double residual_norm = std::sqrt(initial_residual_sq);
  if (previous_objective.isNotNull() &&
      std::abs(residual_norm - precision_scale) <=
        8.0L * std::numeric_limits<double>::epsilon() * precision_scale) {
    // Only at the branch boundary, recover the directly accumulated norm:
    // rounding a long-double sum to the stored double objective must not
    // choose a different numerical path from the five-argument interface.
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
        as<double>(old["lambda"]), alpha, 0);
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
    // Unroll the three components for four-mode tensors. The residual is
    // held fixed only during this preliminary contraction pass.
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
    // Generic order/rank: reuse the contiguous mode-zero fibre contraction.
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

  // Refresh caches and the dense residual in one pass, adding component
  // changes in the requested GS order.  Products constant along a contiguous
  // fibre are hoisted, including the rank-three/four-mode case.
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
      // Ordered pointers remove component lookup and Rcpp access from the
      // innermost loop for every rank, not only the rank-three dispatch.
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
