// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// splmds_stress_from_coords_cpp  (internal, not exported)
//
//   Normalized Kruskal stress between lower-triangle dissimilarities and
//   pairwise Euclidean distances computed from two coordinate columns.
//
//   c1, c2   : length-P coordinate vectors (the two embedding dimensions)
//   diss_vec : length P*(P-1)/2, lower-triangle dissimilarities stored in
//              column-major order (j < i), matching R's lower.tri / dist()
// ---------------------------------------------------------------------------
static double splmds_stress_from_coords_cpp(const vec& c1,
                                            const vec& c2,
                                            const vec& diss_vec) {
  int P = (int)c1.n_elem;
  int m = P * (P - 1) / 2;

  // Pairwise Euclidean distances, column-major lower triangle (j outer, i inner)
  // This replicates R: dist(cbind(c1, c2), method="euclidean", upper=FALSE)
  vec dist_t(m);
  int idx = 0;
  for (int j = 0; j < P; j++) {
    for (int i = j + 1; i < P; i++) {
      double d1 = c1(i) - c1(j);
      double d2 = c2(i) - c2(j);
      dist_t(idx++) = std::sqrt(d1 * d1 + d2 * d2);
    }
  }

  // sqrt( sum((diss - dist)^2) / sum(diss^2) )
  double denom = dot(diss_vec, diss_vec);
  if (denom == 0.0) return 0.0;
  vec diff = diss_vec - dist_t;
  return std::sqrt(dot(diff, diff) / denom);
}

// ---------------------------------------------------------------------------
// splmds_stress_from_coords  (exported R-callable wrapper)
//
//   Direct translation of the R function of the same name.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double splmds_stress_from_coords(const NumericVector& c1_r,
				 const NumericVector& c2_r,
				 const NumericVector& diss_vec_r) {
  return splmds_stress_from_coords_cpp(as<vec>(c1_r),
                                       as<vec>(c2_r),
                                       as<vec>(diss_vec_r));
}

// ---------------------------------------------------------------------------
// stress_SplMDS
//
//   Full spline-MDS stress summed over all time points.
//
//   xi_vec       : flat vector [xi1 (P*K), xi2 (P*K)] — column-major
//   tid_vec      : 1-based integer time indices into Xmat / Xmat2dev rows
//   dissim_list  : List of P×P dissimilarity matrices (one per time point)
//   P            : number of variables (nodes)
//   init_coord   : P×2 matrix of initialisation coordinates
//   lambda       : temporal smoothness penalty weight
//   Xmat         : basis-function matrix (nrow × K)
//   Xmat2dev     : second-derivative basis matrix (nrow × K)
//   diss_vec_list: optional pre-extracted lower-triangle dissimilarities
//                  (list of length n_t; pass R_NilValue / NULL to compute
//                  them on the fly from dissim_list)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double stress_SplMDS(const NumericVector& xi_vec,
		     const IntegerVector& tid_vec,
		     const List&          dissim_list,
		     int                  P,
		     const NumericMatrix& init_coord,
		     double               lambda,
		     const NumericMatrix& Xmat,
		     const NumericMatrix& Xmat2dev,
		     Nullable<List>       diss_vec_list_r = R_NilValue) {

  int K   = Xmat.ncol();
  int n_t = tid_vec.size();

  // ---- Parse xi_vec → xi1 (P×K) and xi2 (P×K), column-major --------------
  // R: xi1 <- matrix(xi_vec[1:(P*K)],        nrow = P)
  //    xi2 <- matrix(xi_vec[(P*K+1):(2*P*K)], nrow = P)
  mat xi1(P, K), xi2(P, K);
  for (int j = 0; j < K; j++) {
    for (int i = 0; i < P; i++) {
      xi1(i, j) = xi_vec[j * P + i];
      xi2(i, j) = xi_vec[P * K + j * P + i];
    }
  }

  // ---- Initialisation coordinate columns ----------------------------------
  vec ic1(P), ic2(P);
  for (int i = 0; i < P; i++) {
    ic1(i) = init_coord(i, 0);
    ic2(i) = init_coord(i, 1);
  }

  // ---- Subset Xmat / Xmat2dev to the requested time rows (1-based) --------
  mat Xsub(n_t, K), X2sub(n_t, K);
  for (int i = 0; i < n_t; i++) {
    int row = tid_vec[i] - 1;               // convert to 0-based
    for (int j = 0; j < K; j++) {
      Xsub(i, j)  = Xmat(row, j);
      X2sub(i, j) = Xmat2dev(row, j);
    }
  }

  // ---- Coordinate matrices (P×n_t) ----------------------------------------
  // c1_all = xi1 %*% t(Xmat[tid_vec,]) + init_coord[,1]   (broadcast per row)
  // c2_all = xi2 %*% t(Xmat[tid_vec,]) + init_coord[,2]
  mat c1_all = xi1 * Xsub.t();             // P×n_t
  mat c2_all = xi2 * Xsub.t();
  c1_all.each_col() += ic1;
  c2_all.each_col() += ic2;

  // ---- Penalty vector (length n_t) ----------------------------------------
  // z1 = xi1 %*% t(Xmat2dev[tid_vec,]),  z2 = xi2 %*% t(Xmat2dev[tid_vec,])
  // penal_vec = lambda * (colSums(z1^2) + colSums(z2^2))
  mat z1 = xi1 * X2sub.t();               // P×n_t
  mat z2 = xi2 * X2sub.t();
  rowvec penal_vec = lambda * (sum(square(z1), 0) + sum(square(z2), 0));

  // ---- Dissimilarity vectors (lower triangle, column-major) ---------------
  // Pre-extract once (or use caller-supplied cache)
  int m = P * (P - 1) / 2;
  std::vector<vec> diss_vecs(n_t);

  bool has_precomp = diss_vec_list_r.isNotNull();
  List diss_vec_list;
  if (has_precomp) diss_vec_list = as<List>(diss_vec_list_r);

  for (int i = 0; i < n_t; i++) {
    if (has_precomp) {
      diss_vecs[i] = as<vec>(as<NumericVector>(diss_vec_list[i]));
    } else {
      NumericMatrix dm = dissim_list[i];   // P×P
      vec dv(m);
      int idx = 0;
      // column-major lower triangle: j outer, i inner (i > j)
      for (int j = 0; j < P; j++) {
        for (int ii = j + 1; ii < P; ii++) {
          dv(idx++) = dm(ii, j);
        }
      }
      diss_vecs[i] = dv;
    }
  }

  // ---- Sum stress over all time points (na.rm = TRUE) ----------------------
  double total = 0.0;
  for (int i = 0; i < n_t; i++) {
    double s = splmds_stress_from_coords_cpp(c1_all.col(i),
                                             c2_all.col(i),
                                             diss_vecs[i]);
    double contribution = s + penal_vec(i);
    if (!std::isnan(contribution)) total += contribution;
  }

  return total;
}
