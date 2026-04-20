// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// reshape_configs
//   Rebuilds a list of coordinate matrices from a flat parameter vector.
//   vec  : flat vector of all coordinates (column-major per time slice)
//   P    : integer vector, P[t] = number of nodes at time t
//   ndim : number of embedding dimensions
//   Tmax : number of time points
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List reshape_configs(const NumericVector& vec,
		     const IntegerVector& P,
		     int ndim,
		     int Tmax) {

  List configs(Tmax);

  // cumulative sum of node counts (0-indexed: cumP[t] = sum of P[0..t-1])
  IntegerVector cumP(Tmax + 1);
  cumP[0] = 0;
  for (int t = 0; t < Tmax; t++) {
    cumP[t + 1] = cumP[t] + P[t];
  }

  for (int t = 0; t < Tmax; t++) {
    int npt   = P[t];
    int start = cumP[t] * ndim;   // 0-based index into vec
    int len   = npt * ndim;

    // Build P[t] x ndim matrix in column-major order (matches R's matrix())
    NumericMatrix mat(npt, ndim);
    for (int j = 0; j < ndim; j++) {
      for (int i = 0; i < npt; i++) {
        mat(i, j) = vec[start + j * npt + i];
      }
    }
    configs[t] = mat;
  }
  return configs;
}

// ---------------------------------------------------------------------------
// stress_DynMDS
//   Kruskal stress + temporal-penalty stress for dynamic MDS.
//
//   vec      : flat coordinate vector (input to optimizer)
//   diss_list: List of numeric matrices (dissimilarity matrices, may have
//              colnames; dim 0 x 0 means "empty / missing" time slice)
//   P        : integer vector of node counts per time point
//   ndim     : embedding dimensionality
//   Tmax     : number of time points
//   lambda   : temporal penalty weight
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double stress_DynMDS(const NumericVector& vec,
		     const List&         diss_list,
		     const IntegerVector& P,
		     int    ndim,
		     int    Tmax,
		     double lambda) {

  // ---- Rebuild coordinate matrices ----------------------------------------
  List configs = reshape_configs(vec, P, ndim, Tmax);

  double stress = 0.0;

  // ---- Kruskal stress (sum over time slices) --------------------------------
  for (int t = 0; t < Tmax; t++) {
    NumericMatrix delta = diss_list[t];
    int n = delta.nrow();
    if (n == 0) continue;               // empty slice → skip

    NumericMatrix coords = configs[t];  // n x ndim

    // Pairwise Euclidean distances (full n x n, including i==j = 0)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double d2 = 0.0;
        for (int k = 0; k < ndim; k++) {
          double diff = coords(i, k) - coords(j, k);
          d2 += diff * diff;
        }
        double d    = std::sqrt(d2);
        double res  = d - delta(i, j);
        stress += res * res;
      }
    }
  }

  // ---- Temporal penalty ---------------------------------------------------
  //
  // Mirrors the R logic exactly:
  //   1. Build node_list / config_list as working copies.
  //   2. Forward-fill NULL (empty colname) slices.
  //   3. For each t >= 1 with a non-empty dissimilarity matrix, find nodes
  //      common to t and t-1, extract their coordinates in column-major order
  //      (matching R's matrix[logical] recycled indexing), and accumulate
  //      lambda * sum(diff^2).

  // Collect column names for each slice (empty vector = NULL in R)
  std::vector<CharacterVector> node_list(Tmax);
  std::vector<NumericMatrix>   config_list(Tmax);

  for (int t = 0; t < Tmax; t++) {
    NumericMatrix delta = diss_list[t];
    // colnames() returns R_NilValue if absent → CharacterVector of length 0
    RObject cn_obj = colnames(delta);
    if (cn_obj.isNULL()) {
      node_list[t] = CharacterVector(0);
    } else {
      node_list[t] = as<CharacterVector>(cn_obj);
    }
    config_list[t] = as<NumericMatrix>(configs[t]);
  }

  // Forward-fill: if node_list[i] is empty, copy from i-1
  for (int i = 1; i < Tmax; i++) {
    if (node_list[i].size() == 0) {
      node_list[i]   = node_list[i - 1];
      config_list[i] = config_list[i - 1];
    }
  }

  // Compute penalty for t = 1 .. Tmax-1
  for (int t = 1; t < Tmax; t++) {
    NumericMatrix delta_t = diss_list[t];
    if (delta_t.nrow() == 0) continue;     // empty slice → skip

    const CharacterVector& nt   = node_list[t];
    const CharacterVector& nt_1 = node_list[t - 1];

    int sz_t   = nt.size();
    int sz_t_1 = nt_1.size();

    // rid_t[i]   = true iff nt[i]   is found in nt_1  (R: node_t   %in% node_t_1)
    // rid_t_1[j] = true iff nt_1[j] is found in nt    (R: node_t_1 %in% node_t)
    std::vector<bool> rid_t(sz_t, false);
    std::vector<bool> rid_t_1(sz_t_1, false);

    for (int i = 0; i < sz_t; i++) {
      for (int j = 0; j < sz_t_1; j++) {
        if (nt[i] == nt_1[j]) { rid_t[i] = true; break; }
      }
    }
    for (int j = 0; j < sz_t_1; j++) {
      for (int i = 0; i < sz_t; i++) {
        if (nt_1[j] == nt[i]) { rid_t_1[j] = true; break; }
      }
    }

    // Replicate R's matrix[logical_vector] column-major recycled extraction:
    //   config_list[[t]][rid_t] → iterate columns, then rows within each col,
    //   keeping elements where rid_t (recycled) is TRUE.
    // Because rid_t has exactly nrow(config_list[[t]]) elements, it recycles
    // perfectly once per column → we get all matched rows across all columns.
    const NumericMatrix& cfg_t   = config_list[t];
    const NumericMatrix& cfg_t_1 = config_list[t - 1];

    // Collect values column-major for matched rows at time t
    std::vector<double> vals_t, vals_t_1;
    vals_t.reserve(sz_t * ndim);
    vals_t_1.reserve(sz_t_1 * ndim);

    for (int d = 0; d < ndim; d++) {
      for (int i = 0; i < sz_t; i++) {
        if (rid_t[i]) vals_t.push_back(cfg_t(i, d));
      }
    }
    for (int d = 0; d < ndim; d++) {
      for (int j = 0; j < sz_t_1; j++) {
        if (rid_t_1[j]) vals_t_1.push_back(cfg_t_1(j, d));
      }
    }

    // Both vectors must be the same length (one entry per matched node × ndim)
    int m = (int)std::min(vals_t.size(), vals_t_1.size());
    for (int k = 0; k < m; k++) {
      double diff = vals_t[k] - vals_t_1[k];
      stress += lambda * diff * diff;
    }
  }

  return stress;
}
