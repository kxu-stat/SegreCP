#include <Rcpp.h>
#include <algorithm>
#include <cmath>
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
    // Explicit loops remove dynamic indexing in the main benchmark regime.
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
