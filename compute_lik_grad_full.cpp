// compute_lik_grad_full.cpp
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// 0-based findInterval: idx with grid[idx] <= x < grid[idx+1]; -1 / n-1 for out-of-range
inline int find_interval(const NumericVector &grid, double x) {
  int n = grid.size();
  if (n < 2) return -1;
  if (x < grid[0]) return -1;
  if (x >= grid[n-1]) return n-1;
  int lo = 0, hi = n-1;
  while (hi - lo > 1) {
    int mid = (lo + hi) >> 1;
    if (grid[mid] <= x) lo = mid; else hi = mid;
  }
  return lo;
}

// Bilinear gradient sampler.
// Input: locs (n x 2), covlist: list of covs each with x (vec), y (vec), z (matrix nrow=length(x), ncol=length(y)).
// Output: NumericVector with dim = c(n_cov, n, 2) storing (df/dx, df/dy).
NumericVector bilinearGradVec_cpp(const NumericMatrix &locs, const List &covlist) {
  int n_cov = covlist.size();
  int n_obs = locs.nrow();
  NumericVector out((size_t)n_cov * n_obs * 2);
  out.attr("dim") = IntegerVector::create(n_cov, n_obs, 2);
  
  for (int j = 0; j < n_cov; ++j) {
    List cov = covlist[j];
    NumericVector x_grid = cov["x"];
    NumericVector y_grid = cov["y"];
    NumericMatrix z      = cov["z"]; // dims: length(x) x length(y)
    
    int nx = x_grid.size();
    int ny = y_grid.size();
    for (int ii = 0; ii < n_obs; ++ii) {
      double lx = locs(ii, 0);
      double ly = locs(ii, 1);
      
      int ix = find_interval(x_grid, lx);
      int iy = find_interval(y_grid, ly);
      
      // invalid cell -> set grads to 0 (safer than NA for downstream algebra)
      if (ix < 0 || iy < 0 || ix >= nx-1 || iy >= ny-1) {
        out[j * n_obs * 2 + ii * 2 + 0] = 0.0;
        out[j * n_obs * 2 + ii * 2 + 1] = 0.0;
        continue;
      }
      
      double x1 = x_grid[ix],     x2 = x_grid[ix + 1];
      double y1 = y_grid[iy],     y2 = y_grid[iy + 1];
      double dx = x2 - x1,        dy = y2 - y1;
      
      double f11 = z(ix,     iy);
      double f21 = z(ix + 1, iy);
      double f12 = z(ix,     iy + 1);
      double f22 = z(ix + 1, iy + 1);
      
      double dfdx = ((y2 - ly) * (f21 - f11) + (ly - y1) * (f22 - f12)) / (dy * dx);
      double dfdy = ((x2 - lx) * (f12 - f11) + (lx - x1) * (f22 - f21)) / (dy * dx);
      
      out[j * n_obs * 2 + ii * 2 + 0] = R_finite(dfdx) ? dfdx : 0.0;
      out[j * n_obs * 2 + ii * 2 + 1] = R_finite(dfdy) ? dfdy : 0.0;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_lik_grad_full_cpp(
    const Rcpp::NumericMatrix &full_x,   // M x (N+2)
    const Rcpp::NumericMatrix &full_y,   // M x (N+2)
    const Rcpp::NumericVector &L_k,      // length M
    const Rcpp::NumericMatrix &X,        // n x 2  (x,y at segment start)
    int i_index,                         // 1-based segment index i
    const Rcpp::NumericVector &par,      // [beta_1 .. beta_p, s] with p >= 1
    double delta,
    int N,
    const Rcpp::List &covlist
) {
  using namespace Rcpp;
  
  const int M    = full_x.nrow();
  const int cols = full_x.ncol();
  if (cols != N + 2) stop("full_x must have N+2 columns");
  if (full_y.nrow() != M || full_y.ncol() != cols) stop("full_y dims mismatch");
  if ((int)L_k.size() != M) stop("L_k length mismatch");
  
  const int i = i_index - 1;
  if (i < 0 || i >= X.nrow() - 1) stop("i out of bounds relative to X");
  
  // grad at start point (s=0)
  NumericMatrix loc0(1, 2);
  loc0(0, 0) = X(i, 0);
  loc0(0, 1) = X(i, 1);
  NumericVector grad0_raw = bilinearGradVec_cpp(loc0, covlist); // (n_cov, 1, 2)
  IntegerVector gdim = grad0_raw.attr("dim");
  if (gdim.size() != 3 || gdim[2] != 2) stop("grad0_raw has unexpected dim");
  const int n_cov = gdim[0];
  
  // number of covariates actually used (from par)
  if ((int)par.size() < 2) stop("par must have at least 2 entries (beta_1, s)");
  const int p = std::min<int>(n_cov, (int)par.size() - 1);
  
  // betas and variance parameter s = gamma^2
  std::vector<double> beta(p);
  for (int k = 0; k < p; ++k) beta[k] = par[k];
  const double s = par[p];
  if (!(R_finite(s) && s > 0.0)) stop("s must be positive");
  
  // grad_0: (n_cov x 2) -> slice first p rows
  NumericMatrix grad_0(p, 2);
  for (int k = 0; k < p; ++k) {
    grad_0(k, 0) = grad0_raw[k * 2 + 0];
    grad_0(k, 1) = grad0_raw[k * 2 + 1];
  }
  
  // constants
  const int Tsteps = N + 1;
  const double factor_u      = (delta * s) / (2.0 * Tsteps);
  const double var_step      = (delta * s) / double(Tsteps);
  const double log_norm_const= std::log(2.0 * M_PI * var_step);
  
  // output rows: 1 (w) + p (S1) + 1 (score_s) + p^2 (S2) + 1 (sumsq)
  const int out_rows = 1 + p + 1 + p*p + 1;
  NumericMatrix out(out_rows, M);
  
  // loop over paths
  for (int j = 0; j < M; ++j) {
    // locations (N x 2) for interior points
    NumericMatrix locs(N, 2);
    for (int sidx = 0; sidx < N; ++sidx) {
      locs(sidx, 0) = full_x(j, sidx + 1);
      locs(sidx, 1) = full_y(j, sidx + 1);
    }
    NumericVector grads_raw = bilinearGradVec_cpp(locs, covlist); // (n_cov, N, 2)
    IntegerVector gdim2 = grads_raw.attr("dim");
    if (gdim2.size() != 3 || gdim2[0] < p || gdim2[1] != N || gdim2[2] != 2)
      stop("grads_raw has unexpected dim");
    
    // u_s for s=0..N
    std::vector<double> usx(Tsteps), usy(Tsteps);
    double u0x = 0.0, u0y = 0.0;
    for (int k = 0; k < p; ++k) {
      u0x += beta[k] * grad_0(k, 0);
      u0y += beta[k] * grad_0(k, 1);
    }
    u0x *= factor_u; u0y *= factor_u;
    usx[0] = u0x;    usy[0] = u0y;
    
    for (int sidx = 0; sidx < N; ++sidx) {
      double sumx = 0.0, sumy = 0.0;
      for (int k = 0; k < p; ++k) {
        // grads_raw(k, sidx, {0,1})
        const double gx = grads_raw[k * N * 2 + sidx * 2 + 0];
        const double gy = grads_raw[k * N * 2 + sidx * 2 + 1];
        sumx += beta[k] * (R_finite(gx) ? gx : 0.0);
        sumy += beta[k] * (R_finite(gy) ? gy : 0.0);
      }
      usx[sidx + 1] = factor_u * sumx;
      usy[sidx + 1] = factor_u * sumy;
    }
    
    // residuals and sumsq
    const int Dlen = 2 * Tsteps;
    std::vector<double> D(Dlen);
    double sumsq = 0.0;
    for (int sidx = 0; sidx < Tsteps; ++sidx) {
      const double dx = full_x(j, sidx + 1) - full_x(j, sidx) - usx[sidx];
      const double dy = full_y(j, sidx + 1) - full_y(j, sidx) - usy[sidx];
      D[2 * sidx + 0] = dx;
      D[2 * sidx + 1] = dy;
      sumsq += dx * dx + dy * dy;
    }
    
    // log-density product across steps
    double logdens = 0.0;
    for (int sidx = 0; sidx < Tsteps; ++sidx) {
      const double dx = D[2 * sidx + 0];
      const double dy = D[2 * sidx + 1];
      const double q  = (dx * dx + dy * dy) / (2.0 * var_step);
      logdens += -q - log_norm_const;
    }
    const double prod_density = std::exp(logdens);
    
    // G (p x Dlen) flattened by cov then time
    std::vector<double> gvec(p * Dlen, 0.0);
    for (int k = 0; k < p; ++k) {              // step 0 from grad_0
      gvec[k * Dlen + 0] = grad_0(k, 0);
      gvec[k * Dlen + 1] = grad_0(k, 1);
    }
    for (int sidx = 0; sidx < N; ++sidx) {     // steps 1..N from grads_raw
      const int idx_x = 2 * (sidx + 1);
      const int idx_y = idx_x + 1;
      for (int k = 0; k < p; ++k) {
        const double gx = grads_raw[k * N * 2 + sidx * 2 + 0];
        const double gy = grads_raw[k * N * 2 + sidx * 2 + 1];
        gvec[k * Dlen + idx_x] = R_finite(gx) ? gx : 0.0;
        gvec[k * Dlen + idx_y] = R_finite(gy) ? gy : 0.0;
      }
    }
    
    // S1 = G %*% D  (length p)
    std::vector<double> gdotD(p, 0.0);
    for (int k = 0; k < p; ++k) {
      double acc = 0.0;
      const double* gk = &gvec[k * Dlen];
      for (int t = 0; t < Dlen; ++t) acc += gk[t] * D[t];
      gdotD[k] = acc;
    }
    
    // inner = (G beta)^T D
    double inner = 0.0;
    for (int t = 0; t < Dlen; ++t) {
      double v = 0.0;
      for (int k = 0; k < p; ++k) v += gvec[k * Dlen + t] * beta[k];
      inner += v * D[t];
    }
    
    // score for s = gamma^2
    const double score_s =
      -double(Tsteps) / s
      + double(Tsteps) / (2.0 * delta * s * s) * sumsq
      + (1.0 / (2.0 * s)) * inner;
      
      // S2 = delta * G G^T  (p x p), flattened row-major
      std::vector<double> S2(p * p, 0.0);
      for (int a = 0; a < p; ++a) {
        for (int b = 0; b < p; ++b) {
          double acc = 0.0;
          const double* ga = &gvec[a * Dlen];
          const double* gb = &gvec[b * Dlen];
          for (int t = 0; t < Dlen; ++t) acc += ga[t] * gb[t];
          S2[a * p + b] = delta * acc;
        }
      }
      
      // fill output
      int r = 0;
      out(r++, j) = prod_density * L_k[j];       // w
      for (int k = 0; k < p; ++k) out(r++, j) = gdotD[k]; // S1 (p rows)
      out(r++, j) = score_s;                     // score for s
      for (int idx = 0; idx < p * p; ++idx) out(r++, j) = S2[idx]; // S2 (row-major)
      out(r++, j) = sumsq;                       // sumsq
  }
  
  return out;
}



// compute_log_lik_grad_full_cpp_p3.cpp

// ---- column-major index helpers for arrays with dim = (n_cov, N, 2) ----
inline int idx_grad(int k, int sidx, int comp, int n_cov, int N) {
  // comp in {0:x, 1:y}
  return k + n_cov * (sidx + N * comp);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_log_lik_grad_full_cpp(
    const Rcpp::NumericMatrix &full_x,   // M x (N+2)
    const Rcpp::NumericMatrix &full_y,   // M x (N+2)
    const Rcpp::NumericVector &log_L_k,  // length M, *log* proposal weights
    const Rcpp::NumericMatrix &X,        // n x 2  (x,y at segment start)
    int i_index,                         // 1-based segment index i
    const Rcpp::NumericVector &par,      // [beta_1 .. beta_p, s] with p >= 1
    double delta,
    int N,
    const Rcpp::List &covlist
) {
  using namespace Rcpp;
  
  const int M    = full_x.nrow();
  const int cols = full_x.ncol();
  if (cols != N + 2) stop("full_x must have N+2 columns");
  if (full_y.nrow() != M || full_y.ncol() != cols) stop("full_y dims mismatch");
  if ((int)log_L_k.size() != M) stop("log_L_k length mismatch");
  
  const int i = i_index - 1;
  if (i < 0 || i >= X.nrow() - 1) stop("i out of bounds relative to X");
  
  // grad at start point (s=0)
  NumericMatrix loc0(1, 2);
  loc0(0, 0) = X(i, 0);
  loc0(0, 1) = X(i, 1);
  NumericVector grad0_raw = bilinearGradVec_cpp(loc0, covlist); // (n_cov, 1, 2)
  IntegerVector gdim = grad0_raw.attr("dim");
  if (gdim.size() != 3 || gdim[2] != 2) stop("grad0_raw has unexpected dim");
  const int n_cov = gdim[0];
  
  // number of covariates actually used (from par)
  if ((int)par.size() < 2) stop("par must have at least 2 entries (beta_1, s)");
  const int p = std::min<int>(n_cov, (int)par.size() - 1);
  
  // betas and variance parameter s = gamma^2
  std::vector<double> beta(p);
  for (int k = 0; k < p; ++k) beta[k] = par[k];
  const double s = par[p];
  if (!(R_finite(s) && s > 0.0)) stop("s must be positive");
  
  // grad_0: (n_cov x 2) -> slice first p rows
  NumericMatrix grad_0(p, 2);
  for (int k = 0; k < p; ++k) {
    grad_0(k, 0) = grad0_raw[k * 2 + 0];
    grad_0(k, 1) = grad0_raw[k * 2 + 1];
  }
  
  // constants
  const int Tsteps = N + 1;
  const double factor_u      = (delta * s) / (2.0 * Tsteps);
  const double var_step      = (delta * s) / double(Tsteps);
  const double log_norm_const= std::log(2.0 * M_PI * var_step);
  
  // output rows: 1 (log_w) + p (S1) + 1 (score_s)
  const int out_rows = 1 + p + 1;
  NumericMatrix out(out_rows, M);
  
  // loop over paths
  for (int j = 0; j < M; ++j) {
    // locations (N x 2) for interior points
    NumericMatrix locs(N, 2);
    for (int sidx = 0; sidx < N; ++sidx) {
      locs(sidx, 0) = full_x(j, sidx + 1);
      locs(sidx, 1) = full_y(j, sidx + 1);
    }
    NumericVector grads_raw = bilinearGradVec_cpp(locs, covlist); // (n_cov, N, 2)
    IntegerVector gdim2 = grads_raw.attr("dim");
    if (gdim2.size() != 3 || gdim2[0] < p || gdim2[1] != N || gdim2[2] != 2)
      stop("grads_raw has unexpected dim");
    
    // u_s for s=0..N
    std::vector<double> usx(Tsteps), usy(Tsteps);
    double u0x = 0.0, u0y = 0.0;
    for (int k = 0; k < p; ++k) {
      u0x += beta[k] * grad_0(k, 0);
      u0y += beta[k] * grad_0(k, 1);
    }
    u0x *= factor_u; u0y *= factor_u;
    usx[0] = u0x;    usy[0] = u0y;
    
    for (int sidx = 0; sidx < N; ++sidx) {
      double sumx = 0.0, sumy = 0.0;
      for (int k = 0; k < p; ++k) {
        const double gx = grads_raw[k * N * 2 + sidx * 2 + 0];
        const double gy = grads_raw[k * N * 2 + sidx * 2 + 1];
        sumx += beta[k] * (R_finite(gx) ? gx : 0.0);
        sumy += beta[k] * (R_finite(gy) ? gy : 0.0);
      }
      usx[sidx + 1] = factor_u * sumx;
      usy[sidx + 1] = factor_u * sumy;
    }
    
    // residuals and sumsq
    const int Dlen = 2 * Tsteps;
    std::vector<double> D(Dlen);
    double sumsq = 0.0;
    for (int sidx = 0; sidx < Tsteps; ++sidx) {
      const double dx = full_x(j, sidx + 1) - full_x(j, sidx) - usx[sidx];
      const double dy = full_y(j, sidx + 1) - full_y(j, sidx) - usy[sidx];
      D[2 * sidx + 0] = dx;
      D[2 * sidx + 1] = dy;
      sumsq += dx * dx + dy * dy;
    }
    
    // log-density sum across steps
    double logdens = 0.0;
    for (int sidx = 0; sidx < Tsteps; ++sidx) {
      const double dx = D[2 * sidx + 0];
      const double dy = D[2 * sidx + 1];
      const double q  = (dx * dx + dy * dy) / (2.0 * var_step);
      logdens += -q - log_norm_const;
    }
    // log-weight: log p(path | par) + log proposal correction
    double logw = logdens + log_L_k[j];
    
    // G (p x Dlen) flattened by cov then time
    std::vector<double> gvec(p * Dlen, 0.0);
    for (int k = 0; k < p; ++k) {              // step 0 from grad_0
      gvec[k * Dlen + 0] = grad_0(k, 0);
      gvec[k * Dlen + 1] = grad_0(k, 1);
    }
    for (int sidx = 0; sidx < N; ++sidx) {     // steps 1..N from grads_raw
      const int idx_x = 2 * (sidx + 1);
      const int idx_y = idx_x + 1;
      for (int k = 0; k < p; ++k) {
        const double gx = grads_raw[k * N * 2 + sidx * 2 + 0];
        const double gy = grads_raw[k * N * 2 + sidx * 2 + 1];
        gvec[k * Dlen + idx_x] = R_finite(gx) ? gx : 0.0;
        gvec[k * Dlen + idx_y] = R_finite(gy) ? gy : 0.0;
      }
    }
    
    // S1 = G %*% D  (length p)
    std::vector<double> gdotD(p, 0.0);
    for (int k = 0; k < p; ++k) {
      double acc = 0.0;
      const double* gk = &gvec[k * Dlen];
      for (int t = 0; t < Dlen; ++t) acc += gk[t] * D[t];
      gdotD[k] = acc;
    }
    
    // inner = (G beta)^T D
    double inner = 0.0;
    for (int t = 0; t < Dlen; ++t) {
      double v = 0.0;
      for (int k = 0; k < p; ++k) v += gvec[k * Dlen + t] * beta[k];
      inner += v * D[t];
    }
    
    // score for s = gamma^2
    const double score_s =
      -double(Tsteps) / s
      + double(Tsteps) / (2.0 * delta * s * s) * sumsq
      + (1.0 / (2.0 * s)) * inner;
      
      // fill output
      int r = 0;
      out(r++, j) = logw;              // row 1: log weight
      for (int k = 0; k < p; ++k) {
        out(r++, j) = gdotD[k];        // rows 2..(p+1): S1
      }
      out(r++, j) = score_s;           // row p+2: score for s
  }
  
  return out;
}












