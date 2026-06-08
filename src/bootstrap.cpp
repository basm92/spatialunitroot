#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Eigen-based robust alternatives to chol() and inv_sympd()
inline arma::mat eigen_chol(const arma::mat& X_in) {
    arma::mat X = 0.5 * (X_in + X_in.t());
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, X);
    arma::vec d = arma::sqrt(eigval);
    return eigvec * arma::diagmat(d) * eigvec.t();
}

inline arma::mat eigen_inv(const arma::mat& X_in) {
    arma::mat X = 0.5 * (X_in + X_in.t());
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, X);
    return eigvec * arma::diagmat(1.0 / eigval) * eigvec.t();
}

//' Find c such that mean(exp(-c * lower_tri(distmat))) equals rhobar
//'
//' Uses exponential bracketing followed by geometric bisection.
//' Matches the Stata getcbar() function exactly.
//' @param rhobar Target correlation
//' @param distmat n x n distance matrix
//' @return The constant cbar
// [[Rcpp::export]]
double getcbar(double rhobar, const arma::mat& distmat) {
    int n = distmat.n_rows;

    // Extract lower triangular elements (excluding diagonal)
    int nlower = (n * (n - 1)) / 2;
    arma::vec vd(nlower);
    int idx = 0;
    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            vd(idx) = distmat(i, j);
            idx++;
        }
    }

    double c0 = 10.0;
    double c1 = 10.0;
    double v, cm, cbar;

    // Find c0 such that v > rhobar
    int i1 = 0;
    int jj = 0;
    while (i1 == 0) {
        v = arma::mean(arma::exp(-c0 * vd));
        i1 = (v > rhobar) ? 1 : 0;
        if (i1 == 0) {
            c1 = c0;
            c0 = c0 / 2.0;
            jj++;
        }
        if (jj > 500) {
            Rcpp::stop("getcbar: rhobar too large");
        }
    }

    // Find c1 such that v < rhobar
    i1 = 0;
    jj = 0;
    while (i1 == 0) {
        v = arma::mean(arma::exp(-c1 * vd));
        i1 = (v < rhobar) ? 1 : 0;
        if (i1 == 0) {
            c0 = c1;
            c1 = 2.0 * c1;
            jj++;
        }
        if (c1 > 10000.0) {
            i1 = 1;
        }
        if (jj > 500) {
            Rcpp::stop("getcbar: rhobar too small");
        }
    }

    // Geometric bisection
    while ((c1 - c0) > 0.001) {
        cm = std::sqrt(c0 * c1);
        v = arma::mean(arma::exp(-cm * vd));
        if (v < rhobar) {
            c1 = cm;
        } else {
            c0 = cm;
        }
    }

    cbar = std::sqrt(c0 * c1);
    return cbar;
}

//' Compute power of quadratic form test
//'
//' Compares quadratic forms under null (om0) vs alternative (om1) covariance.
//' Matches the Stata getpow_qf() function.
//' @param om0 Null hypothesis covariance matrix (q x q)
//' @param om1 Alternative hypothesis covariance matrix (q x q)
//' @param e q x nrep matrix of standard normal draws
//' @return Power (scalar between 0 and 1)
// [[Rcpp::export]]
double getpow_qf(const arma::mat& om0_in, const arma::mat& om1_in,
                  const arma::mat& e) {
    // Force symmetry
    arma::mat om0 = 0.5 * (om0_in + om0_in.t());
    arma::mat om1 = 0.5 * (om1_in + om1_in.t());

    arma::mat om0i = eigen_inv(om0);
    arma::mat om1i = eigen_inv(om1);

    // Eigen-based Cholesky decompositions (tolerates near-symmetry)
    arma::mat ch_om0 = eigen_chol(om0);
    arma::mat ch_om1 = eigen_chol(om1);
    arma::mat ch_om0i = eigen_chol(om0i);
    arma::mat ch_om1i = eigen_chol(om1i);

    // ho = ch_om1i * ch_om0'
    arma::mat ho = ch_om1i * ch_om0.t();
    // ha = ch_om0i * ch_om1'
    arma::mat ha = ch_om0i * ch_om1.t();

    // Quadratic forms
    arma::rowvec qe = arma::sum(arma::square(e), 0);

    arma::mat ya_o = ho * e;
    arma::mat yo_a = ha * e;

    arma::rowvec qa_o = arma::sum(arma::square(ya_o), 0);
    arma::rowvec qo_a = arma::sum(arma::square(yo_a), 0);

    // LR statistics
    arma::rowvec lr_o = qe / qa_o;
    arma::rowvec lr_a = qo_a / qe;

    // Critical value (95th percentile under H0)
    arma::vec prob = {0.95};
    double cv = arma::as_scalar(arma::quantile(lr_o, prob));

    // Power = proportion of H1 draws exceeding CV
    double pow = arma::mean(lr_a > cv);

    return pow;
}
