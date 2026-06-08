#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Double-center a matrix (subtract row means, then column means)
//'
//' Matches the Stata demeanmat() function.
//' @param mat n x n matrix
//' @return Double-centered matrix
// [[Rcpp::export]]
arma::mat double_center(const arma::mat& mat) {
    int n = mat.n_rows;
    arma::mat result = mat;

    // Subtract row means
    result = result - (arma::sum(result, 1) / n) * arma::ones<arma::rowvec>(n);

    // Subtract column means
    result = result - arma::ones<arma::vec>(n) * (arma::sum(result, 0) / n);

    return result;
}

//' Compute LBM covariance matrix from distance matrix
//'
//' sigma_lbm = 0.5 * (J * dist[,1]' + dist[1,] * J' - dist)
//' Uses the first location as the origin.
//' @param distmat n x n distance matrix
//' @return n x n LBM covariance matrix
// [[Rcpp::export]]
arma::mat sigma_lbm(const arma::mat& distmat) {
    int n = distmat.n_rows;
    arma::mat sig(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sig(i, j) = 0.5 * (distmat(i, 0) + distmat(0, j) - distmat(i, j));
        }
    }

    return sig;
}

//' Compute demeaned LBM covariance matrix
//'
//' sigma_lbm_dm = demeanmat(sigma_lbm(distmat))
//' @param distmat n x n distance matrix
//' @return Double-centered LBM covariance matrix
// [[Rcpp::export]]
arma::mat sigma_lbm_dm(const arma::mat& distmat) {
    arma::mat sig = sigma_lbm(distmat);
    return double_center(sig);
}

//' Compute demeaned spatial covariance with exponential decay
//'
//' sigma = exp(-c * distmat), then double-center.
//' @param distmat n x n distance matrix
//' @param c Decay parameter
//' @return Double-centered exponential covariance matrix
// [[Rcpp::export]]
arma::mat sigma_dm(const arma::mat& distmat, double c) {
    arma::mat sigma = arma::exp(-c * distmat);
    return double_center(sigma);
}

//' Compute residual-projected spatial covariance
//'
//' sigma = exp(-c * distmat), then sigma_dm = M * sigma * M'.
//' @param distmat n x n distance matrix
//' @param c Decay parameter
//' @param M n x n projection matrix (I - X(X'X)^{-1}X')
//' @return Residual-projected covariance matrix
// [[Rcpp::export]]
arma::mat sigma_residual(const arma::mat& distmat, double c,
                          const arma::mat& M) {
    arma::mat sigma = arma::exp(-c * distmat);
    return M * sigma * M.t();
}

//' Get top q eigenvectors of a symmetric matrix
//'
//' Computes eigendecomposition, sorts eigenvalues descending,
//' returns first q eigenvectors.
//' Matches the Stata get_R() function.
//' @param sigma n x n symmetric matrix
//' @param q Number of eigenvectors to extract
//' @return n x q matrix of top eigenvectors
// [[Rcpp::export]]
arma::mat get_R(const arma::mat& sigma, int q) {
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, sigma);

    // Sort by descending eigenvalues
    arma::uvec idx = arma::sort_index(eigval, "descend");
    arma::vec sorted_eigval = eigval(idx);
    arma::mat sorted_eigvec = eigvec.cols(idx);

    // Return first q eigenvectors
    return sorted_eigvec.cols(0, q - 1);
}

//' Compute LBM-GLS transformation matrix
//'
//' From sigma_lbm_dm eigendecomposition, keep eigenvalues > 1e-10,
//' return V * diag(1/sqrt(eval)) * V'.
//' @param distmat n x n distance matrix
//' @return n x n LBM-GLS transformation matrix
// [[Rcpp::export]]
arma::mat lbm_gls_matrix(const arma::mat& distmat) {
    double small = 1.0e-10;

    arma::mat sig = sigma_lbm_dm(distmat);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sig);

    // Sort descending
    arma::uvec idx = arma::sort_index(eigval, "descend");
    eigval = eigval(idx);
    eigvec = eigvec.cols(idx);

    // Keep eigenvalues > small
    arma::uvec keep = arma::find(eigval > small);
    arma::vec evals_keep = eigval(keep);
    arma::mat evecs_keep = eigvec.cols(keep);

    // Dsi = diag(1 / sqrt(eval))
    arma::vec dsi = 1.0 / arma::sqrt(evals_keep);
    arma::mat result = evecs_keep * arma::diagmat(dsi) * evecs_keep.t();

    return result;
}

//' Compute nearest-neighbor transformation matrix
//'
//' NN_mat = I - rowmins_mat, where each row has 1/k for the k nearest neighbors.
//' @param s n x d matrix of locations
//' @param latlong Whether coordinates are latitude/longitude
//' @return n x n NN transformation matrix
// [[Rcpp::export]]
arma::mat nn_matrix(const arma::mat& s, bool latlong) {
    int n = s.n_rows;
    arma::mat distmat;

    if (latlong) {
        double cc = M_PI / 180.0;
        distmat.set_size(n, n);
        arma::vec lat = s.col(0);
        arma::vec lon = s.col(1);
        for (int i = 0; i < n; i++) {
            arma::vec dlat = 0.5 * cc * (lat - lat(i));
            arma::vec dlon = 0.5 * cc * (lon - lon(i));
            double cos_lat_i = std::cos(cc * lat(i));
            arma::vec term = arma::square(arma::sin(dlat)) +
                             cos_lat_i * arma::cos(cc * lat) %
                             arma::square(arma::sin(dlon));
            distmat.col(i) = arma::asin(arma::sqrt(term)) / M_PI;
        }
    } else {
        distmat = arma::mat(n, n);
        for (int i = 0; i < n; i++) {
            arma::rowvec si = s.row(i);
            arma::mat diff = s.each_row() - si;
            distmat.col(i) = arma::sqrt(arma::sum(arma::square(diff), 1));
        }
    }

    // Normalize
    double maxdist = distmat.max();
    if (maxdist > 0) distmat = distmat / maxdist;

    // Replace zeros with large values for finding minimum non-zero distance
    arma::mat distmat_nozeros = distmat;
    distmat_nozeros.replace(0.0, 1e10);

    // Find row minima
    arma::vec rowmins = arma::min(distmat_nozeros, 1);

    // Build rowmins_mat: 1/k for each nearest neighbor
    arma::mat rowmins_mat(n, n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (std::abs(distmat_nozeros(i, j) - rowmins(i)) < 1e-12) {
                rowmins_mat(i, j) = 1.0;
            }
        }
    }

    // Normalize rows
    arma::vec rowsums = arma::sum(rowmins_mat, 1);
    for (int i = 0; i < n; i++) {
        if (rowsums(i) > 0) {
            rowmins_mat.row(i) = rowmins_mat.row(i) / rowsums(i);
        }
    }

    return arma::eye<arma::mat>(n, n) - rowmins_mat;
}

//' Compute isotropic transformation matrix
//'
//' iso_mat = I - dist_below_b / rowsum(dist_below_b)
//' where dist_below_b indicates distances <= radius
//' @param s n x d matrix of locations
//' @param radius Radius for differencing
//' @param latlong Whether coordinates are latitude/longitude
//' @return n x n isotropic transformation matrix
// [[Rcpp::export]]
arma::mat iso_matrix(const arma::mat& s, double radius, bool latlong) {
    int n = s.n_rows;
    arma::mat distmat;

    if (latlong) {
        double cc = M_PI / 180.0;
        distmat.set_size(n, n);
        arma::vec lat = s.col(0);
        arma::vec lon = s.col(1);
        for (int i = 0; i < n; i++) {
            arma::vec dlat = 0.5 * cc * (lat - lat(i));
            arma::vec dlon = 0.5 * cc * (lon - lon(i));
            double cos_lat_i = std::cos(cc * lat(i));
            arma::vec term = arma::square(arma::sin(dlat)) +
                             cos_lat_i * arma::cos(cc * lat) %
                             arma::square(arma::sin(dlon));
            distmat.col(i) = arma::asin(arma::sqrt(term)) / M_PI;
        }
        // Convert to metres: multiply by pi * 6371000.009 * 2
        distmat = distmat * M_PI * 6371000.009 * 2;
    } else {
        distmat = arma::mat(n, n);
        for (int i = 0; i < n; i++) {
            arma::rowvec si = s.row(i);
            arma::mat diff = s.each_row() - si;
            distmat.col(i) = arma::sqrt(arma::sum(arma::square(diff), 1));
        }
    }

    // Replace zeros with large values
    arma::mat distmat_nozeros = distmat;
    distmat_nozeros.replace(0.0, 1e10);

    // Find distances below or equal to radius
    arma::mat dist_below_b(n, n, arma::fill::zeros);
    int num_no_neighbors = 0;
    for (int i = 0; i < n; i++) {
        bool has_neighbor = false;
        for (int j = 0; j < n; j++) {
            if (distmat_nozeros(i, j) <= radius) {
                dist_below_b(i, j) = 1.0;
                has_neighbor = true;
            }
        }
        if (!has_neighbor) num_no_neighbors++;
    }

    // Normalize rows
    arma::vec rowsums = arma::sum(dist_below_b, 1);
    for (int i = 0; i < n; i++) {
        if (rowsums(i) > 0) {
            dist_below_b.row(i) = dist_below_b.row(i) / rowsums(i);
        }
    }

    return arma::eye<arma::mat>(n, n) - dist_below_b;
}

//' Compute cluster transformation matrix
//'
//' clust_mat = I - (clust_id == clust_id') / rowsum(clust_id == clust_id')
//' @param cluster n x 1 vector of integer cluster IDs
//' @return n x n cluster transformation matrix
// [[Rcpp::export]]
arma::mat cluster_matrix(const arma::vec& cluster) {
    int n = cluster.n_elem;
    arma::mat clust(n, n);

    // Build cluster membership indicator matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            clust(i, j) = (cluster(i) == cluster(j)) ? 1.0 : 0.0;
        }
    }

    // Normalize rows by cluster size
    arma::vec rowsums = arma::sum(clust, 1);
    for (int i = 0; i < n; i++) {
        if (rowsums(i) > 0) {
            clust.row(i) = clust.row(i) / rowsums(i);
        }
    }

    return arma::eye<arma::mat>(n, n) - clust;
}

//' Apply transformation matrix to a variable
//'
//' hy = H * y. For lbmgls, additionally mean-center the result.
//' @param y n x 1 vector
//' @param H n x n transformation matrix
//' @param is_lbmgls Whether to mean-center the result
//' @return n x 1 transformed vector
// [[Rcpp::export]]
arma::vec apply_transform(const arma::vec& y, const arma::mat& H,
                           bool is_lbmgls) {
    arma::vec hy = H * y;
    if (is_lbmgls) {
        hy = hy - arma::mean(hy);
    }
    return hy;
}

//' Half-vectorization: extract lower triangular elements (excluding diagonal)
//'
//' Matches the Stata lvech() function.
//' @param S n x n matrix
//' @return Column vector of lower-triangular elements
// [[Rcpp::export]]
arma::vec lvech(const arma::mat& S) {
    int n = S.n_rows;
    int nlower = (n * (n - 1)) / 2;
    arma::vec v(nlower);
    int idx = 0;
    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            v(idx) = S(i, j);
            idx++;
        }
    }
    return v;
}
