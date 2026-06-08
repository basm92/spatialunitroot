#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Compute Euclidean distance matrix
//'
//' @param coords n x d matrix of coordinates
//' @return n x n Euclidean distance matrix
// [[Rcpp::export]]
arma::mat euclidean_distances(const arma::mat& coords) {
    int n = coords.n_rows;
    arma::mat distmat(n, n);

    for (int i = 0; i < n; i++) {
        arma::rowvec si = coords.row(i);
        arma::mat diff = coords.each_row() - si;
        distmat.col(i) = arma::sqrt(arma::sum(arma::square(diff), 1));
    }

    return distmat;
}

//' Compute Haversine great-circle distance matrix
//'
//' Returns asin(sqrt(haversine_term)) / pi for each pair, i.e. the central
//' angle scaled by 1 / (2 * pi). This scaling is arbitrary because the tests
//' normalize distances (dividing by the maximum), so any constant factor
//' cancels. The great-circle distance in metres is recovered downstream by
//' multiplying the raw value by 2 * pi * earth_radius (see spurhalflife()).
//'
//' @param coords n x 2 matrix with latitude in column 0, longitude in column 1
//' @return n x n matrix of central angles scaled by 1 / (2 * pi)
// [[Rcpp::export]]
arma::mat haversine_distances(const arma::mat& coords) {
    int n = coords.n_rows;
    double c = M_PI / 180.0;
    arma::mat distmat(n, n);
    arma::vec lat = coords.col(0);
    arma::vec lon = coords.col(1);

    for (int i = 0; i < n; i++) {
        arma::vec dlat = 0.5 * c * (lat - lat(i));
        arma::vec dlon = 0.5 * c * (lon - lon(i));

        arma::vec sin_dlat = arma::sin(dlat);
        arma::vec sin_dlon = arma::sin(dlon);
        double cos_lat_i = std::cos(c * lat(i));
        arma::vec cos_lat = arma::cos(c * lat);

        arma::vec term = arma::square(sin_dlat) +
                         cos_lat_i * cos_lat % arma::square(sin_dlon);

        // distance = asin(sqrt(term)) / pi
        distmat.col(i) = arma::asin(arma::sqrt(term)) / M_PI;
    }

    return distmat;
}

//' Normalize distance matrix so maximum distance is 1
//'
//' @param distmat n x n distance matrix
//' @return Normalized distance matrix
// [[Rcpp::export]]
arma::mat normalize_distances(const arma::mat& distmat) {
    double maxdist = distmat.max();
    if (maxdist > 0) {
        return distmat / maxdist;
    }
    return distmat;
}
