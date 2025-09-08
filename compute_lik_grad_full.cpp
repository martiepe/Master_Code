// compute_lik_grad_full.cpp
#include <Rcpp.h>
using namespace Rcpp;

// Helper: find interval index like R::findInterval but returns 0-based index
inline int find_interval(const NumericVector &grid, double x) {
  int n = grid.size();
  if (x < grid[0]) return -1;
  if (x >= grid[n-1]) return n-1; // mark as out-of-right-bound
  // binary search: want idx such that grid[idx] <= x < grid[idx+1]
  int lo = 0, hi = n-1;
  while (hi - lo > 1) {
    int mid = (lo + hi) >> 1;
    if (grid[mid] <= x) lo = mid;
    else hi = mid;
  }
  return lo; // 0-based
}

// bilinearGradVec implemented in C++
// locs: n x 2 matrix of locations
// covlist: list of covs where each cov is a list with x (vec), y (vec), z (matrix with nrow = length(x))
//
// returns NumericVector with dim = c(n_cov, n, 2) (same ordering as your R function)
NumericVector bilinearGradVec_cpp(const NumericMatrix &locs, const List &covlist) {
  int n_cov = covlist.size();
  int n_obs = locs.nrow();
  NumericVector out((size_t)n_cov * n_obs * 2);
  IntegerVector dims = IntegerVector::create(n_cov, n_obs, 2);
  out.attr("dim") = dims;
  
  for (int j = 0; j < n_cov; ++j) {
    List cov = covlist[j];
    NumericVector x_grid = cov["x"];
    NumericVector y_grid = cov["y"];
    NumericMatrix z = cov["z"]; // in R z has nrow = length(x_grid)
    int nx = x_grid.size();
    int ny = y_grid.size();
    
    for (int ii = 0; ii < n_obs; ++ii) {
      double lx = locs(ii, 0);
      double ly = locs(ii, 1);
      
      int ix = find_interval(x_grid, lx);
      int iy = find_interval(y_grid, ly);
      
      // valid if 0 <= ix < nx-1 and 0 <= iy < ny-1
      if (ix <= -1 || iy <= -1 || ix >= nx - 1 || iy >= ny - 1) {
        // leave as NA-like: we'll set to NaN which mirrors NA_real_
        out[j * n_obs * 2 + ii * 2 + 0] = NA_REAL;
        out[j * n_obs * 2 + ii * 2 + 1] = NA_REAL;
        continue;
      }
      
      double x1 = x_grid[ix];
      double x2 = x_grid[ix + 1];
      double y1 = y_grid[iy];
      double y2 = y_grid[iy + 1];
      
      double dx = x2 - x1;
      double dy = y2 - y1;
      
      // z is indexed: z(row = ix+1 in R, col = iy+1 in R), but in C++ it's 0-based: z(ix, iy)
      double f11 = z(ix,     iy);
      double f21 = z(ix + 1, iy);
      double f12 = z(ix,     iy + 1);
      double f22 = z(ix + 1, iy + 1);
      
      double dfdx = ((y2 - ly) * (f21 - f11) + (ly - y1) * (f22 - f12)) / (dy * dx);
      double dfdy = ((x2 - lx) * (f12 - f11) + (lx - x1) * (f22 - f21)) / (dy * dx);
      
      out[j * n_obs * 2 + ii * 2 + 0] = dfdx;
      out[j * n_obs * 2 + ii * 2 + 1] = dfdy;
    } // ii
  } // j
  
  return out;
}


// [[Rcpp::export]]
NumericMatrix compute_lik_grad_full_cpp(
    const NumericMatrix &full_x,    // M x (N+2)
    const NumericMatrix &full_y,    // M x (N+2)
    const NumericVector &L_k,       // length M  (importance weight as in your R)
    const NumericMatrix &X,         // full X path used to compute grad_0 (nrow >= i+1)
    int i_index,                    // 1-based i index (we'll convert inside)
    const NumericVector &par,       // length >= 4
    double delta,
    int N,
    const List &covlist
) {
  // Basic checks
  int M = full_x.nrow();
  int cols = full_x.ncol();
  if (cols != N + 2) stop("full_x must have N+2 columns");
  if (full_y.nrow() != M || full_y.ncol() != cols) stop("full_y dims mismatch");
  if ((int)L_k.size() != M) stop("L_k length mismatch");
  if (par.size() < 4) stop("par must have length >= 4");
  
  // Convert i_index from 1-based (R) to 0-based for indexing into X
  int i = i_index - 1;
  if (i < 0 || i >= X.nrow() - 1) stop("i out of bounds relative to X");
  
  // Precompute grad_0 (1 x 2) but bilinearGradVec returns 3 x 1 x 2 style structure;
  NumericMatrix loc0(1, 2);
  loc0(0,0) = X(i, 0);
  loc0(0,1) = X(i, 1);
  NumericVector grad0_raw = bilinearGradVec_cpp(loc0, covlist); // dim (n_cov, 1, 2)
  IntegerVector gdim = grad0_raw.attr("dim");
  int n_cov = gdim[0];
  
  // Extract grad_0 as a n_cov x 2 matrix
  NumericMatrix grad_0(n_cov, 2);
  for (int j = 0; j < n_cov; ++j) {
    grad_0(j, 0) = grad0_raw[j * 1 * 2 + 0];
    grad_0(j, 1) = grad0_raw[j * 1 * 2 + 1];
  }
  
  // Prepare output 5 x M (rows: prod_density*L_k, then 3 components g%*%D, then scalar term)
  NumericMatrix out(5, M);
  
  double par1 = par[0], par2 = par[1], par3 = par[2], par4 = par[3];
  double factor_u = (delta * par4 / ((N + 1) * 2.0));
  double var = delta * par4 / double(N + 1);
  double log_norm_const = std::log(2.0 * M_PI * var);
  
  // For each sample j:
  for (int j = 0; j < M; ++j) {
    // Build grads for N points: we need locs = N x 2 matrix of (x_samples[j,], y_samples[j,])
    NumericMatrix locs(N, 2);
    for (int s = 0; s < N; ++s) {
      locs(s, 0) = full_x(j, s+1); // note full_x row contains [X[i,1], x_samples..., X[i+1,1]]
      locs(s, 1) = full_y(j, s+1);
    }
    NumericVector grads_raw = bilinearGradVec_cpp(locs, covlist); // dim (n_cov, N, 2)
    IntegerVector gdim2 = grads_raw.attr("dim"); // {n_cov, N, 2}
    
    // Build 'us' for s=0..N: us is (N+1) x 2 where s=0 uses u_0
    std::vector<double> usx(N+1), usy(N+1);
    // u_0 from grad_0: it's a 1 x 2 vector for each cov combined by par1..par3
    // compute u_0_x = factor_u * sum_k par_k * grad_0[k, 0], etc.
    double u0x = 0.0, u0y = 0.0;
    for (int k = 0; k < n_cov; ++k) {
      double g0x = grad_0(k, 0);
      double g0y = grad_0(k, 1);
      u0x += par1 * g0x + par2 * g0x * 0.0 + par3 * g0x * 0.0; // keep formula but only par1..par3 apply per cov (like R)
      // The R code sums par[1]*grad1 + par[2]*grad2 + par[3]*grad3 by cov index,
      // Here we assume covlist[[1..3]] correspond to the three covariates and par1..par3 map to them.
      // BUT the indexing below (see filling us for s>=1) uses k separately, so above we should do:
    }
    // Correction: compute u0 properly:
    u0x = 0.0; u0y = 0.0;
    for (int k = 0; k < n_cov; ++k) {
      double g0x_k = grad_0(k, 0);
      double g0y_k = grad_0(k, 1);
      double pk = 0.0;
      if (k == 0) pk = par1;
      else if (k == 1) pk = par2;
      else if (k == 2) pk = par3;
      else pk = 0.0;
      u0x += pk * g0x_k;
      u0y += pk * g0y_k;
    }
    u0x *= factor_u; u0y *= factor_u;
    usx[0] = u0x; usy[0] = u0y;
    
    // For steps s = 1..N compute from grads
    for (int s = 0; s < N; ++s) {
      double sumx = 0.0, sumy = 0.0;
      for (int k = 0; k < n_cov; ++k) {
        // fetch grads_raw element for (k, s, :)
        double gx = grads_raw[k * N * 2 + s * 2 + 0];
        double gy = grads_raw[k * N * 2 + s * 2 + 1];
        double pk = 0.0;
        if (k == 0) pk = par1;
        else if (k == 1) pk = par2;
        else if (k == 2) pk = par3;
        else pk = 0.0;
        // add contribution
        sumx += pk * gx;
        sumy += pk * gy;
      }
      usx[s+1] = factor_u * sumx;
      usy[s+1] = factor_u * sumy;
    }
    
    // Build D vector and subtract us
    int Dlen = 2 * (N + 1);
    std::vector<double> D(Dlen);
    double sumsq = 0.0;
    for (int s = 0; s <= N; ++s) {
      double dx = full_x(j, s+1) - full_x(j, s);
      double dy = full_y(j, s+1) - full_y(j, s);
      dx -= usx[s];
      dy -= usy[s];
      D[2*s + 0] = dx;
      D[2*s + 1] = dy;
      sumsq += dx*dx + dy*dy;
    }
    
    // compute log product density for the N+1 bivariate normals
    double logdens = 0.0;
    for (int s = 0; s <= N; ++s) {
      double dx = D[2*s + 0];
      double dy = D[2*s + 1];
      double q = (dx*dx + dy*dy) / (2.0 * var);
      logdens += - q - log_norm_const;
    }
    double prod_density = std::exp(logdens);
    
    // Build g (3 x Dlen): first two cols come from grad_0; then N blocks from grads_raw
    // g(k, idx) where k=0..2 (cov index)
    std::vector<double> gvec(3 * Dlen, 0.0);
    // first step s=0 uses grad_0(k, :)
    for (int k = 0; k < n_cov && k < 3; ++k) {
      gvec[k * Dlen + 0] = grad_0(k, 0); // x
      gvec[k * Dlen + 1] = grad_0(k, 1); // y
    }
    // subsequent steps s = 1..N take grads_raw
    for (int s = 0; s < N; ++s) {
      for (int k = 0; k < n_cov && k < 3; ++k) {
        double gx = grads_raw[k * N * 2 + s * 2 + 0];
        double gy = grads_raw[k * N * 2 + s * 2 + 1];
        int idx_x = 2 * (s + 1);
        int idx_y = idx_x + 1;
        gvec[k * Dlen + idx_x] = gx;
        gvec[k * Dlen + idx_y] = gy;
      }
    }
    
    // compute g %*% D -> yields 3 entries
    NumericVector gdotD(3);
    for (int k = 0; k < 3; ++k) {
      double acc = 0.0;
      for (int t = 0; t < Dlen; ++t) acc += gvec[k * Dlen + t] * D[t];
      gdotD[k] = acc;
    }
    
    // compute temp = t(g) %*% par[1:3] -> Dlen vector where temp[t] = sum_k g[k,t] * par[k]
    std::vector<double> temp(Dlen, 0.0);
    for (int t = 0; t < Dlen; ++t) {
      double v = 0.0;
      if (n_cov > 0) v += gvec[0 * Dlen + t] * par1;
      if (n_cov > 1) v += gvec[1 * Dlen + t] * par2;
      if (n_cov > 2) v += gvec[2 * Dlen + t] * par3;
      temp[t] = v;
    }
    // scalar inner = temp' * D
    double inner = 0.0;
    for (int t = 0; t < Dlen; ++t) inner += temp[t] * D[t];
    
    double term3 = -double(N + 1) / par4 + double(N + 1) / (2.0 * delta * par4 * par4) * sumsq
    + (1.0 / (2.0 * par4)) * inner;
    
    // Fill out the columns
    out(0, j) = prod_density * L_k[j];
    out(1, j) = gdotD[0];
    out(2, j) = gdotD[1];
    out(3, j) = gdotD[2];
    out(4, j) = term3;
  } // j
  
  return out;
}
