#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Forward declarations of functions in other compilation units
arma::mat sigma_dm(const arma::mat& distmat, double c);
arma::mat sigma_residual(const arma::mat& distmat, double c, const arma::mat& M);
arma::mat get_R(const arma::mat& sigma, int q);
arma::mat sigma_lbm_dm(const arma::mat& distmat);
double getcbar(double rhobar, const arma::mat& distmat);
double getpow_qf(const arma::mat& om0, const arma::mat& om1, const arma::mat& e);

// Eigen-based robust alternatives to chol() and inv_sympd()
// These use eig_sym which tolerates near-symmetric matrices better

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

// ---- I(1) Test ----

//' Spatial I(1) unit root test
//'
//' Tests H0: spatial unit root (I(1)) vs H1: not I(1).
//' @param Y n x 1 vector of observations
//' @param distmat n x n normalized distance matrix
//' @param emat q x nrep matrix of standard normal draws
//' @param q Number of low-frequency components
//' @return List with LR, pvalue, ha_parm, cv_vec
// [[Rcpp::export]]
Rcpp::List spatial_i1_test(const arma::vec& Y, const arma::mat& distmat,
                            const arma::mat& emat, int q) {
    int n = distmat.n_rows;

    // BM covariance matrix (demeaned)
    arma::mat sigdm_bm = sigma_lbm_dm(distmat);

    // Low-frequency eigenvectors
    arma::mat R = get_R(sigdm_bm, q);

    // Null hypothesis covariance
    arma::mat om_ho = R.t() * sigdm_bm * R;

    // Find ha_parm with ~50% power
    double pow50 = 0.5;
    double pow = 1.0;
    double ctry = getcbar(0.95, distmat);
    double c, c1, c2;
    arma::mat sigdm_c, om_c;
    int maxiter = 20;

    // Step 1: decrease c until power < 0.5
    while (pow > pow50) {
        c = ctry;
        sigdm_c = sigma_dm(distmat, c);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);
        ctry = ctry / 2.0;
    }
    c1 = c;

    // Step 2: increase c until power > 0.5
    pow = 0.0;
    ctry = getcbar(0.01, distmat);
    while (pow < pow50) {
        c = ctry;
        sigdm_c = sigma_dm(distmat, c);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);
        ctry = 2.0 * ctry;
    }
    c2 = c;

    // Step 3: bisection
    int iter = 0;
    while (std::abs(pow - pow50) > 0.01) {
        c = (c1 + c2) / 2.0;
        sigdm_c = sigma_dm(distmat, c);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);

        if (pow > pow50) {
            c2 = c;
        } else if (pow < pow50) {
            c1 = c;
        }
        iter++;
        if (iter > maxiter) break;
    }

    double ha_parm = c;

    // Alternative hypothesis covariance
    arma::mat sigdm_ha = sigma_dm(distmat, ha_parm);
    arma::mat om_ha = R.t() * sigdm_ha * R;

    // Use eigen helper functions (handles symmetry internally)
    arma::mat ch_om_ho  = eigen_chol(om_ho);
    arma::mat omi_ho    = eigen_inv(om_ho);
    arma::mat omi_ha    = eigen_inv(om_ha);
    arma::mat ch_omi_ho = eigen_chol(omi_ho);
    arma::mat ch_omi_ha = eigen_chol(omi_ha);

    // Null distribution of LR
    int nrep = emat.n_cols;
    arma::mat y_ho = ch_om_ho.t() * emat;  // q x nrep
    arma::mat y_ho_ho = ch_omi_ho * y_ho;   // q x nrep
    arma::mat y_ho_ha = ch_omi_ha * y_ho;   // q x nrep
    arma::rowvec q_ho_ho = arma::sum(arma::square(y_ho_ho), 0);
    arma::rowvec q_ho_ha = arma::sum(arma::square(y_ho_ha), 0);
    arma::rowvec lr_ho = q_ho_ho / q_ho_ha;

    // Critical values (1%, 5%, 10%)
    arma::vec sz_vec = {0.01, 0.05, 0.10};
    arma::vec cv_vec = arma::quantile(lr_ho, 1.0 - sz_vec);

    // Test statistic for data
    arma::vec X = Y - arma::mean(Y);
    arma::vec P = R.t() * X;
    double LR = arma::as_scalar((P.t() * omi_ho * P) / (P.t() * omi_ha * P));
    double pvalue = arma::mean(lr_ho > LR);

    return Rcpp::List::create(
        Rcpp::Named("LR") = LR,
        Rcpp::Named("pvalue") = pvalue,
        Rcpp::Named("ha_parm") = ha_parm,
        Rcpp::Named("cv_vec") = cv_vec
    );
}

// ---- I(0) Test ----

//' Spatial I(0) unit root test
//'
//' Tests H0: no spatial unit root (I(0)) vs H1: spatial unit root.
//' Uses a grid of rho values for the null hypothesis.
//' @param Y n x 1 vector of observations
//' @param distmat n x n normalized distance matrix
//' @param emat q x nrep matrix of standard normal draws
//' @param q Number of low-frequency components
//' @return List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
// [[Rcpp::export]]
Rcpp::List spatial_i0_test(const arma::vec& Y, const arma::mat& distmat,
                            const arma::mat& emat, int q) {
    int n = distmat.n_rows;
    int nrep = emat.n_cols;

    // BM covariance for weighting
    double rho_bm = 0.999;
    double c_bm = getcbar(rho_bm, distmat);
    arma::mat sigdm_bm = sigma_dm(distmat, c_bm);
    arma::mat R = get_R(sigdm_bm, q);

    // om_rho with rho = 0.001 (white noise approximation)
    double rho_init = 0.001;
    double c_init = getcbar(rho_init, distmat);
    arma::mat sigdm_rho = sigma_dm(distmat, c_init);
    arma::mat om_rho = R.t() * sigdm_rho * R;
    arma::mat om_bm = R.t() * sigdm_bm * R;

    // Find ha_parm g (~50% power)
    arma::mat om_i0 = om_rho;
    arma::mat om_ho = om_rho;

    double pow50 = 0.5;
    double pow = 1.0;
    double gtry = 1.0;
    double g, g1, g2;
    int maxiter = 20;

    // Step 1: decrease g until power < 0.5
    while (pow > pow50) {
        g = gtry;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        gtry = g / 2.0;
    }
    g1 = g;

    // Step 2: increase g until power > 0.5
    pow = 0.0;
    gtry = 30.0;
    while (pow < pow50) {
        g = gtry;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        gtry = g * 2.0;
    }
    g2 = g;

    // Step 3: bisection
    int iter = 0;
    while (std::abs(pow - pow50) > 0.01) {
        g = (g1 + g2) / 2.0;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        if (pow > pow50) {
            g2 = g;
        } else if (pow < pow50) {
            g1 = g;
        }
        iter++;
        if (iter > maxiter) break;
    }

    double ha_parm = g;
    arma::mat om_ha = om_i0 + ha_parm * om_bm;

    // Force symmetry
    om_ho = 0.5 * (om_ho + om_ho.t());
    om_ha = 0.5 * (om_ha + om_ha.t());

    // Cholesky for null and alternative
    arma::mat ch_omi_ho = eigen_chol(eigen_inv(om_ho));
    arma::mat ch_omi_ha = eigen_chol(eigen_inv(om_ha));

    // LR for data
    arma::vec X = Y - arma::mean(Y);
    arma::vec P = R.t() * X;
    arma::vec y_P_ho = ch_omi_ho * P;
    arma::vec y_P_ha = ch_omi_ha * P;
    double q_P_ho = arma::as_scalar(arma::sum(arma::square(y_P_ho)));
    double q_P_ha = arma::as_scalar(arma::sum(arma::square(y_P_ha)));
    double LR = q_P_ho / q_P_ha;

    // Grid of rho values (log-spaced, 30 points)
    double rho_min = 0.0001;
    double rho_max = 0.03;
    int n_rho = 30;
    arma::vec rho_grid(n_rho);
    for (int i = 0; i < n_rho; i++) {
        rho_grid(i) = std::exp(std::log(rho_min) +
                               (double)i / (n_rho - 1) *
                               (std::log(rho_max) - std::log(rho_min)));
    }

    // Store Cholesky of om_ho for each rho
    std::vector<arma::mat> ch_om_ho_list(n_rho);
    for (int i = 0; i < n_rho; i++) {
        double rho = rho_grid(i);
        arma::mat om_ho_rho;
        if (rho > 0) {
            double c_rho = getcbar(rho, distmat);
            arma::mat sigdm_ho = sigma_dm(distmat, c_rho);
            om_ho_rho = R.t() * sigdm_ho * R;
        } else {
            om_ho_rho = arma::eye<arma::mat>(q, q);
        }
        om_ho_rho = 0.5 * (om_ho_rho + om_ho_rho.t());
        ch_om_ho_list[i] = eigen_chol(om_ho_rho);
    }

    // Compute p-values and critical values across grid
    arma::vec sz_vec = {0.01, 0.05, 0.10};
    arma::mat pvalue_mat(n_rho, 1);
    arma::mat cvalue_mat(n_rho, 3);

    for (int ir = 0; ir < n_rho; ir++) {
        arma::mat ch_om_ho_grid = ch_om_ho_list[ir];

        arma::mat y_ho = ch_om_ho_grid.t() * emat;
        arma::mat y_ho_ho = ch_omi_ho * y_ho;
        arma::mat y_ho_ha = ch_omi_ha * y_ho;
        arma::rowvec q_ho_ho = arma::sum(arma::square(y_ho_ho), 0);
        arma::rowvec q_ho_ha = arma::sum(arma::square(y_ho_ha), 0);
        arma::rowvec lr_ho = q_ho_ho / q_ho_ha;

        arma::vec cv_vec = arma::quantile(lr_ho, 1.0 - sz_vec);
        cvalue_mat.row(ir) = cv_vec.t();
        pvalue_mat(ir, 0) = arma::mean(lr_ho > LR);
    }

    // colmax over rho grid
    double pvalue = pvalue_mat.max();
    arma::vec cvalue = arma::max(cvalue_mat, 0).t();

    return Rcpp::List::create(
        Rcpp::Named("LR") = LR,
        Rcpp::Named("pvalue") = pvalue,
        Rcpp::Named("cvalue") = cvalue,
        Rcpp::Named("ha_parm") = ha_parm,
        Rcpp::Named("rho_grid") = rho_grid,
        Rcpp::Named("pvalue_mat") = pvalue_mat,
        Rcpp::Named("cvalue_mat") = cvalue_mat
    );
}

// ---- I(1) Residual Test ----

//' Spatial I(1) residual unit root test
//'
//' Tests residuals from regression Y on X for spatial unit root.
//' @param Y n x 1 response vector
//' @param X_in n x k design matrix (including intercept)
//' @param distmat n x n normalized distance matrix
//' @param emat q x nrep matrix of standard normal draws
//' @param q Number of low-frequency components
//' @return List with LR, pvalue, ha_parm, cv_vec
// [[Rcpp::export]]
Rcpp::List spatial_i1_test_residual(const arma::vec& Y,
                                      const arma::mat& X_in,
                                      const arma::mat& distmat,
                                      const arma::mat& emat, int q) {
    int n = distmat.n_rows;

    // Projection matrix M = I - X*(X'X)^{-1}*X'
    arma::mat XtX_inv = eigen_inv(X_in.t() * X_in);
    arma::mat M = arma::eye<arma::mat>(n, n) - X_in * XtX_inv * X_in.t();

    // BM covariance with residual projection
    double rho_bm = 0.999;
    double c_bm = getcbar(rho_bm, distmat);
    arma::mat sigdm_bm = sigma_residual(distmat, c_bm, M);
    arma::mat R = get_R(sigdm_bm, q);

    arma::mat om_ho = R.t() * sigdm_bm * R;

    // Find ha_parm
    double pow50 = 0.5;
    double pow = 1.0;
    double ctry = getcbar(0.95, distmat);
    double c, c1, c2;
    arma::mat sigdm_c, om_c;
    int maxiter = 20;

    while (pow > pow50) {
        c = ctry;
        sigdm_c = sigma_residual(distmat, c, M);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);
        ctry = ctry / 2.0;
    }
    c1 = c;

    pow = 0.0;
    ctry = getcbar(0.01, distmat);
    while (pow < pow50) {
        c = ctry;
        sigdm_c = sigma_residual(distmat, c, M);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);
        ctry = 2.0 * ctry;
    }
    c2 = c;

    int iter = 0;
    while (std::abs(pow - pow50) > 0.01) {
        c = (c1 + c2) / 2.0;
        sigdm_c = sigma_residual(distmat, c, M);
        om_c = R.t() * sigdm_c * R;
        pow = getpow_qf(om_ho, om_c, emat);
        if (pow > pow50) { c2 = c; }
        else if (pow < pow50) { c1 = c; }
        iter++;
        if (iter > maxiter) break;
    }

    double ha_parm = c;
    arma::mat sigdm_ha = sigma_residual(distmat, ha_parm, M);
    arma::mat om_ha = R.t() * sigdm_ha * R;

    // Force symmetry
    om_ho = 0.5 * (om_ho + om_ho.t());
    om_ha = 0.5 * (om_ha + om_ha.t());

    // Test statistic
    arma::mat ch_om_ho = eigen_chol(om_ho);
    arma::mat omi_ho = eigen_inv(om_ho);
    arma::mat omi_ha = eigen_inv(om_ha);
    arma::mat ch_omi_ho = eigen_chol(omi_ho);
    arma::mat ch_omi_ha = eigen_chol(omi_ha);

    int nrep = emat.n_cols;
    arma::mat y_ho = ch_om_ho.t() * emat;
    arma::mat y_ho_ho = ch_omi_ho * y_ho;
    arma::mat y_ho_ha = ch_omi_ha * y_ho;
    arma::rowvec q_ho_ho = arma::sum(arma::square(y_ho_ho), 0);
    arma::rowvec q_ho_ha = arma::sum(arma::square(y_ho_ha), 0);
    arma::rowvec lr_ho = q_ho_ho / q_ho_ha;

    arma::vec sz_vec = {0.01, 0.05, 0.10};
    arma::vec cv_vec = arma::quantile(lr_ho, 1.0 - sz_vec);

    arma::vec Xd = Y - arma::mean(Y);
    arma::vec P = R.t() * Xd;
    double LR = arma::as_scalar((P.t() * omi_ho * P) / (P.t() * omi_ha * P));
    double pvalue = arma::mean(lr_ho > LR);

    return Rcpp::List::create(
        Rcpp::Named("LR") = LR,
        Rcpp::Named("pvalue") = pvalue,
        Rcpp::Named("ha_parm") = ha_parm,
        Rcpp::Named("cv_vec") = cv_vec
    );
}

// ---- I(0) Residual Test ----

//' Spatial I(0) residual unit root test
//'
//' Tests residuals from regression Y on X for no spatial unit root (I(0)).
//' @param Y n x 1 response vector
//' @param X_in n x k design matrix (including intercept)
//' @param distmat n x n normalized distance matrix
//' @param emat q x nrep matrix of standard normal draws
//' @param q Number of low-frequency components
//' @return List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
// [[Rcpp::export]]
Rcpp::List spatial_i0_test_residual(const arma::vec& Y,
                                      const arma::mat& X_in,
                                      const arma::mat& distmat,
                                      const arma::mat& emat, int q) {
    int n = distmat.n_rows;
    int nrep = emat.n_cols;

    // Projection matrix
    arma::mat XtX_inv = eigen_inv(X_in.t() * X_in);
    arma::mat M = arma::eye<arma::mat>(n, n) - X_in * XtX_inv * X_in.t();

    // BM covariance for weighting
    double rho_bm = 0.999;
    double c_bm = getcbar(rho_bm, distmat);
    arma::mat sigdm_bm = sigma_residual(distmat, c_bm, M);
    arma::mat R = get_R(sigdm_bm, q);

    // om_rho with rho=0.001
    double rho_init = 0.001;
    double c_init = getcbar(rho_init, distmat);
    arma::mat sigdm_rho = sigma_residual(distmat, c_init, M);
    arma::mat om_rho = R.t() * sigdm_rho * R;
    arma::mat om_bm = R.t() * sigdm_bm * R;

    arma::mat om_i0 = om_rho;
    arma::mat om_ho = om_rho;

    // Find ha_parm g
    double pow50 = 0.5;
    double pow = 1.0;
    double gtry = 1.0;
    double g, g1, g2;
    int maxiter = 20;

    while (pow > pow50) {
        g = gtry;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        gtry = g / 2.0;
    }
    g1 = g;

    pow = 0.0;
    gtry = 30.0;
    while (pow < pow50) {
        g = gtry;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        gtry = g * 2.0;
    }
    g2 = g;

    int iter = 0;
    while (std::abs(pow - pow50) > 0.01) {
        g = (g1 + g2) / 2.0;
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, emat);
        if (pow > pow50) { g2 = g; }
        else if (pow < pow50) { g1 = g; }
        iter++;
        if (iter > maxiter) break;
    }

    double ha_parm = g;
    arma::mat om_ha = om_i0 + ha_parm * om_bm;

    om_ho = 0.5 * (om_ho + om_ho.t());
    om_ha = 0.5 * (om_ha + om_ha.t());

    arma::mat ch_omi_ho = eigen_chol(eigen_inv(om_ho));
    arma::mat ch_omi_ha = eigen_chol(eigen_inv(om_ha));

    // LR for data
    arma::vec X = Y - arma::mean(Y);
    arma::vec P = R.t() * X;
    arma::vec y_P_ho = ch_omi_ho * P;
    arma::vec y_P_ha = ch_omi_ha * P;
    double q_P_ho = arma::as_scalar(arma::sum(arma::square(y_P_ho)));
    double q_P_ha = arma::as_scalar(arma::sum(arma::square(y_P_ha)));
    double LR = q_P_ho / q_P_ha;

    // Grid of rho values
    double rho_min = 0.0001;
    double rho_max = 0.03;
    int n_rho = 30;
    arma::vec rho_grid(n_rho);
    for (int i = 0; i < n_rho; i++) {
        rho_grid(i) = std::exp(std::log(rho_min) +
                               (double)i / (n_rho - 1) *
                               (std::log(rho_max) - std::log(rho_min)));
    }

    std::vector<arma::mat> ch_om_ho_list(n_rho);
    for (int i = 0; i < n_rho; i++) {
        double rho = rho_grid(i);
        arma::mat om_ho_rho;
        if (rho > 0) {
            double c_rho = getcbar(rho, distmat);
            arma::mat sigdm_ho = sigma_residual(distmat, c_rho, M);
            om_ho_rho = R.t() * sigdm_ho * R;
        } else {
            om_ho_rho = arma::eye<arma::mat>(q, q);
        }
        om_ho_rho = 0.5 * (om_ho_rho + om_ho_rho.t());
        ch_om_ho_list[i] = eigen_chol(om_ho_rho);
    }

    arma::vec sz_vec = {0.01, 0.05, 0.10};
    arma::mat pvalue_mat(n_rho, 1);
    arma::mat cvalue_mat(n_rho, 3);

    for (int ir = 0; ir < n_rho; ir++) {
        arma::mat ch_om_ho_grid = ch_om_ho_list[ir];
        arma::mat y_ho = ch_om_ho_grid.t() * emat;
        arma::mat y_ho_ho = ch_omi_ho * y_ho;
        arma::mat y_ho_ha = ch_omi_ha * y_ho;
        arma::rowvec q_ho_ho = arma::sum(arma::square(y_ho_ho), 0);
        arma::rowvec q_ho_ha = arma::sum(arma::square(y_ho_ha), 0);
        arma::rowvec lr_ho = q_ho_ho / q_ho_ha;

        arma::vec cv_vec = arma::quantile(lr_ho, 1.0 - sz_vec);
        cvalue_mat.row(ir) = cv_vec.t();
        pvalue_mat(ir, 0) = arma::mean(lr_ho > LR);
    }

    double pvalue = pvalue_mat.max();
    arma::vec cvalue = arma::max(cvalue_mat, 0).t();

    return Rcpp::List::create(
        Rcpp::Named("LR") = LR,
        Rcpp::Named("pvalue") = pvalue,
        Rcpp::Named("cvalue") = cvalue,
        Rcpp::Named("ha_parm") = ha_parm,
        Rcpp::Named("rho_grid") = rho_grid,
        Rcpp::Named("pvalue_mat") = pvalue_mat,
        Rcpp::Named("cvalue_mat") = cvalue_mat
    );
}

// ---- Half-life / Spatial Persistence ----

//' Compute spatial half-life confidence interval
//'
//' Implements the Bayesian posterior-based confidence set for the
//' spatial half-life parameter. Matches Stata spatial_persistence().
//' @param Z n x 1 vector of observations
//' @param distmat n x n normalized distance matrix
//' @param emat q x nrep matrix of standard normal draws
//' @param level Confidence level (e.g., 0.95)
//' @return List with ci_lower, ci_upper
// [[Rcpp::export]]
Rcpp::List spatial_persistence(const arma::vec& Z, const arma::mat& distmat,
                                const arma::mat& emat, double level) {
    int n = distmat.n_rows;
    int nrep = emat.n_cols;
    int q = emat.n_rows;

    // Half-life grids
    int n_hl_ho1 = 100;
    arma::vec hl_grid_ho1 = arma::linspace(0.001, 1.0, n_hl_ho1);

    int n_hl_ho2 = 30;
    arma::vec hl_grid_ho2 = arma::linspace(1.01, 3.0, n_hl_ho2);

    int n_hl_total = n_hl_ho1 + n_hl_ho2 + 1;
    arma::vec hl_grid_ho = arma::join_cols(
        arma::join_cols(hl_grid_ho1, hl_grid_ho2),
        arma::vec({100.0})
    );

    int n_hl_ha = 50;
    arma::vec hl_grid_ha = arma::linspace(0.001, 1.0, n_hl_ha);

    // Map half-life to decay parameter: c = -log(0.5) / half_life
    double log2 = std::log(0.5);
    arma::vec c_grid_ho = -log2 / hl_grid_ho;
    arma::vec c_grid_ha = -log2 / hl_grid_ha;

    // BM covariance for low-frequency extraction
    double rho_bm = 0.999;
    double c_bm = getcbar(rho_bm, distmat);
    arma::mat sigdm_bm = sigma_dm(distmat, c_bm);
    arma::mat R = get_R(sigdm_bm, q);

    // Precompute Cholesky and constants for null grid
    std::vector<arma::mat> ch_om_ho_vec(n_hl_total);
    std::vector<double> const_den_vec(n_hl_total);

    for (int i = 0; i < n_hl_total; i++) {
        double c = c_grid_ho(i);
        arma::mat sigdm = sigma_dm(distmat, c);
        arma::mat om = R.t() * sigdm * R;
        om = 0.5 * (om + om.t());
        ch_om_ho_vec[i] = eigen_chol(om);

        arma::mat omi = eigen_inv(om);
        double logdet = 0.0;
        double det_val = 0.0;
        arma::log_det(logdet, det_val, omi);
        double const_den = std::sqrt(std::abs(det_val)) * 0.5 *
                           std::tgamma((double)q / 2.0) /
                           std::pow(M_PI, (double)q / 2.0);
        const_den_vec[i] = const_den;
    }

    // Precompute Cholesky constants for alternative grid
    std::vector<arma::mat> ch_omi_ha_vec(n_hl_ha);
    std::vector<double> const_den_ha_vec(n_hl_ha);

    for (int i = 0; i < n_hl_ha; i++) {
        double c = c_grid_ha(i);
        arma::mat sigdm = sigma_dm(distmat, c);
        arma::mat om = R.t() * sigdm * R;
        om = 0.5 * (om + om.t());
        arma::mat omi = eigen_inv(om);
        ch_omi_ha_vec[i] = eigen_chol(omi);

        double logdet = 0.0;
        double det_val = 0.0;
        arma::log_det(logdet, det_val, omi);
        double const_den = std::sqrt(std::abs(det_val)) * 0.5 *
                           std::tgamma((double)q / 2.0) /
                           std::pow(M_PI, (double)q / 2.0);
        const_den_ha_vec[i] = const_den;
    }

    // Data projection
    arma::vec X = R.t() * Z;

    // Compute p-values over the null grid
    arma::vec pv_vec(n_hl_total);

    for (int i = 0; i < n_hl_total; i++) {
        arma::mat ch_null = ch_om_ho_vec[i];
        double const_den = const_den_vec[i];
        arma::mat ch_omi = eigen_chol(eigen_inv(
            ch_null * ch_null.t()));

        // Density under null
        arma::vec Xc = ch_omi * X;
        double sum_sq_X = arma::as_scalar(arma::sum(arma::square(Xc)));
        double den_ho_X = const_den * std::pow(sum_sq_X, -(double)q / 2.0);

        // Average density under alternative
        double den_ha_sum = 0.0;
        for (int j = 0; j < n_hl_ha; j++) {
            arma::mat ch_omi_ha = ch_omi_ha_vec[j];
            double const_den_ha = const_den_ha_vec[j];
            arma::vec Xc_ha = ch_omi_ha * X;
            double sum_sq_ha = arma::as_scalar(arma::sum(arma::square(Xc_ha)));
            double den_ha = const_den_ha *
                            std::pow(sum_sq_ha, -(double)q / 2.0);
            den_ha_sum += den_ha;
        }
        double den_ha_avg = den_ha_sum / n_hl_ha;

        // Null distribution and p-value
        arma::mat e_scaled = ch_null.t() * emat;
        arma::mat e_scaled2 = ch_omi * e_scaled;
        arma::rowvec sum_sq_e = arma::sum(arma::square(e_scaled2), 0);
        arma::rowvec den_ho_e = const_den *
            arma::pow(sum_sq_e, -(double)q / 2.0);

        arma::rowvec den_ha_e_mat(nrep);
        for (int j = 0; j < n_hl_ha; j++) {
            arma::mat ch_omi_ha = ch_omi_ha_vec[j];
            double const_den_ha = const_den_ha_vec[j];
            arma::mat e_ha = ch_omi_ha * e_scaled;
            arma::rowvec sum_sq_ha_e = arma::sum(arma::square(e_ha), 0);
            arma::rowvec den_j = const_den_ha *
                arma::pow(sum_sq_ha_e, -(double)q / 2.0);
            den_ha_e_mat += den_j;
        }
        arma::rowvec den_ha_avg_e = den_ha_e_mat / n_hl_ha;
        arma::rowvec lr_e = den_ha_avg_e / den_ho_e;
        double lr_X = den_ha_avg / den_ho_X;

        pv_vec(i) = arma::mean(lr_e > lr_X);
    }

    // Find confidence bounds: half-life values where p > 1 - level
    arma::uvec ci_idx = arma::find(pv_vec > 1.0 - level);
    double ci_lower, ci_upper;

    if (ci_idx.n_elem > 0) {
        ci_lower = arma::min(hl_grid_ho(ci_idx));
        ci_upper = arma::max(hl_grid_ho(ci_idx));
    } else {
        ci_lower = R_NaN;
        ci_upper = R_NaN;
    }

    return Rcpp::List::create(
        Rcpp::Named("ci_lower") = ci_lower,
        Rcpp::Named("ci_upper") = ci_upper
    );
}
