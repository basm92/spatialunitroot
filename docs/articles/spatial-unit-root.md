# Spatial Unit Root Diagnostics with spatialunitroot

## Introduction

This vignette demonstrates the `spatialunitroot` R package, which
implements the spatial unit root diagnostic tests and transformations of
Müller & Watson (2024, *Econometrica*). The package is an R port of the
SPUR Stata package by Becker, Boll & Voth (2025).

Spatial data often exhibit strong dependence across locations. Just as
time series with unit roots can produce spurious regressions, spatial
processes with (near) unit roots can lead to misleading inference. This
package provides tools to:

1.  **Diagnose** spatial unit roots using I(1) and I(0) tests
2.  **Remove** spatial unit roots via spatial differencing
    transformations
3.  **Characterize** spatial persistence via half-life confidence
    intervals

## The Chetty et al. (2014) Commuting Zone Data

We use the dataset from Chetty, Hendren, Kline & Saez (2014), which
contains measures of intergenerational mobility and its correlates
across 722 US commuting zones. This is the same data used in the Stata
package’s supplementary materials.

``` r
data(chetty)
str(chetty[, 1:10], vec.len = 2)
```

    ## Classes 'tbl_df', 'tbl' and 'data.frame':    722 obs. of  10 variables:
    ##  $ state    : chr  "TN" "TN" ...
    ##  $ cz       : num  100 200 301 302 401 ...
    ##  $ czname   : chr  "Johnson City" "Morristown" ...
    ##  $ s_1      : num  36.5 36 ...
    ##  $ s_2      : num  -82.4 -83.4 ...
    ##  $ am       : num  38.4 37.8 ...
    ##  $ fracblack: num  0.0208 0.0198 ...
    ##  $ racseg   : num  0.0904 0.0932 ...
    ##  $ segpov25 : num  0.0302 0.0279 ...
    ##  $ fraccom15: num  0.325 0.276 ...

Spatial coordinates are stored in `s_1` (latitude) and `s_2`
(longitude), following the convention of Müller & Watson’s `scpc`
package.

## Diagnosing Spatial Unit Roots

We illustrate the procedure using absolute upward mobility (`am`), a
variable known to be highly spatially persistent.

### I(0) Test: Is the variable spatially stationary?

The I(0) test has the null hypothesis of **no** spatial unit root. A low
p-value indicates strong evidence *against* spatial stationarity.

``` r
set.seed(42)
t_i0 <- spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
                 type = "i0", latlong = TRUE, q = 15, nrep = 10000)
summary(t_i0)
```

    ## 
    ## Spatial I(0) Test Summary
    ## =======================================
    ## 
    ## Test Statistic :  3.04364 
    ## P-value        :  0.0011 
    ## HA Parameter   :  41.56 
    ## 
    ## Critical Values:
    ##   1%  :   2.5225
    ##   5%  :   2.1020
    ##   10% :   1.9393
    ## 
    ## Test Configuration:
    ##   Type        :  i0 
    ##   q (components) :  15 
    ##   MC replications:  10000

The p-value is well below 0.01, so we **strongly reject** the null that
`am` is I(0). This is consistent with the Stata results reported in
Becker, Boll & Voth (2025, §4.2), where the I(0) test on a persistent
spatial variable yields p = 0.0004.

### I(1) Test: Is there a spatial unit root?

The I(1) test has the null hypothesis of a spatial unit root. A high
p-value means we cannot reject the presence of a unit root.

``` r
set.seed(42)
t_i1 <- spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
                 type = "i1", latlong = TRUE, q = 15, nrep = 10000)
summary(t_i1)
```

    ## 
    ## Spatial I(1) Test Summary
    ## =======================================
    ## 
    ## Test Statistic :  5.75884 
    ## P-value        :  0.3874 
    ## HA Parameter   :  9.238 
    ## 
    ## Critical Values:
    ##   1%  :   9.7841
    ##   5%  :   8.5872
    ##   10% :   7.8594
    ## 
    ## Test Configuration:
    ##   Type        :  i1 
    ##   q (components) :  15 
    ##   MC replications:  10000

With a p-value of 0.387, we **fail to reject** the spatial unit root
null. Again, this matches the Stata article’s finding (p = 0.574 for
their simulated persistent variable). Together, the I(0) and I(1) tests
provide strong evidence that `am` has a spatial unit root.

### Spatial Half-Life

The half-life measures the distance at which spatial correlation drops
to one-half.

``` r
set.seed(42)
hl <- spurhalflife(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
                   latlong = TRUE, q = 15, nrep = 10000)
summary(hl)
```

    ## 
    ## Spatial Half-Life Confidence Interval Summary
    ## =========================================
    ## 
    ## Confidence level :  95 %
    ## CI Lower bound   :  417743 
    ## CI Upper bound   : inf
    ## Max distance     :  4549672 
    ## Normalized units :  FALSE 
    ## 
    ## Configuration:
    ##   q (components)  :  15 
    ##   MC replications :  10000

The half-life confidence interval for `am` is 417.7 km — indicating very
high spatial persistence (the upper bound is unbounded, a sign of an
(near) unit root).

## Testing Multiple Variables

We can assess spatial persistence across many variables simultaneously.
The table below shows I(1) and I(0) test p-values for several covariates
from Chetty et al. (2014).

``` r
vars <- c("am", "fracblack", "racseg", "segpov25", "hipc", "gini",
           "hsdrop", "colgrad", "manshare", "fracfor")

results <- do.call(rbind, lapply(vars, function(v) {
  set.seed(42)
  t1 <- spurtest(as.formula(paste(v, "~ 1")), data = chetty,
                 coords = ~ s_1 + s_2, type = "i1",
                 latlong = TRUE, q = 15, nrep = 5000)
  t0 <- spurtest(as.formula(paste(v, "~ 1")), data = chetty,
                 coords = ~ s_1 + s_2, type = "i0",
                 latlong = TRUE, q = 15, nrep = 5000)
  hl <- spurhalflife(as.formula(paste(v, "~ 1")), data = chetty,
                     coords = ~ s_1 + s_2,
                     latlong = TRUE, normdist = TRUE,
                     q = 15, nrep = 5000)
  data.frame(
    Variable = v,
    I1_pvalue = round(t1$p_value, 3),
    I0_pvalue = round(t0$p_value, 3),
    HL_lower = round(hl$ci_lower, 3),
    HL_upper = if(is.infinite(hl$ci_upper)) "∞" else round(hl$ci_upper, 3),
    Verdict = if(t1$p_value > 0.05 && t0$p_value < 0.05) "I(1)" else
              if(t1$p_value < 0.05 && t0$p_value > 0.05) "I(0)" else "ambiguous"
  )
}))

knitr::kable(results, caption = "Spatial persistence diagnostics for Chetty et al. variables")
```

| Variable  | I1_pvalue | I0_pvalue | HL_lower | HL_upper | Verdict   |
|:----------|----------:|----------:|---------:|:---------|:----------|
| am        |     0.388 |     0.001 |    0.092 | ∞        | I(1)      |
| fracblack |     0.119 |     0.008 |    0.031 | ∞        | I(1)      |
| racseg    |     0.011 |     0.135 |    0.001 | 0.284    | I(0)      |
| segpov25  |     0.284 |     0.027 |    0.051 | ∞        | I(1)      |
| hipc      |     0.136 |     0.142 |    0.011 | ∞        | ambiguous |
| gini      |     0.780 |     0.000 |    0.263 | ∞        | I(1)      |
| hsdrop    |     0.095 |     0.018 |    0.031 | ∞        | I(1)      |
| colgrad   |     0.039 |     0.023 |    0.001 | ∞        | ambiguous |
| manshare  |     0.212 |     0.003 |    0.062 | ∞        | I(1)      |
| fracfor   |     0.551 |     0.034 |    0.173 | ∞        | I(1)      |

Spatial persistence diagnostics for Chetty et al. variables

Most variables fail to reject the I(1) null and reject the I(0) null,
indicating strong spatial persistence — consistent with the findings in
Becker, Boll & Voth (2025, Appendix A).

## Correcting Spurious Regression

A key application is diagnosing and correcting spurious spatial
regression. We illustrate with a regression of `am` on `fracblack` and
`racseg`.

### Step 1: Naïve OLS

``` r
m_ols <- feols(am ~ fracblack + racseg, data = chetty, se = "standard")
summary(m_ols)
```

    ## OLS estimation, Dep. Var.: am
    ## Observations: 693
    ## Standard-errors: IID 
    ##             Estimate Std. Error   t value   Pr(>|t|)    
    ## (Intercept)  47.8631   0.283542 168.80421  < 2.2e-16 ***
    ## fracblack   -24.4455   1.367710 -17.87331  < 2.2e-16 ***
    ## racseg      -13.4511   1.717636  -7.83115 1.8264e-14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## RMSE: 4.33291   Adj. R2: 0.411517

### Step 2: Test residuals for a spatial unit root

``` r
set.seed(42)
t_resid <- spurtest(m_ols, coords = ~ s_1 + s_2,
                    type = "i1resid", latlong = TRUE,
                    q = 15, nrep = 5000)
t_resid
```

    ## 
    ## Spatial I(1) Residual Test Results
    ## ---------------------------------------
    ## Test Statistic :  1111 
    ## P-value        :  0.192 
    ## ---------------------------------------

The residual test suggests that spatial persistence remains in the
residuals, which can lead to misleading inference.

### Step 3: Apply LBM-GLS transformation

``` r
transformed <- spurtransform(~ am + fracblack + racseg, data = chetty,
                              coords = ~ s_1 + s_2,
                              prefix = "h_", latlong = TRUE)
```

### Step 4: Re-estimate on transformed data

``` r
m_trans <- feols(h_am ~ h_fracblack + h_racseg, data = transformed,
                 se = "standard")
summary(m_trans)
```

    ## OLS estimation, Dep. Var.: h_am
    ## Observations: 693
    ## Standard-errors: IID 
    ##                  Estimate Std. Error       t value   Pr(>|t|)    
    ## (Intercept) -6.940000e-16   0.999411 -6.950000e-16 1.0000e+00    
    ## h_fracblack -1.297061e+01   2.095557 -6.189575e+00 1.0348e-09 ***
    ## h_racseg    -1.151500e+01   1.107396 -1.039827e+01  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## RMSE: 26.3   Adj. R2: 0.222151

The R² drops from 0.413 to 0.224 after removing spatial trends,
indicating that much of the apparent relationship was driven by shared
spatial patterns — a classic spurious regression result.

## Comparison with Stata SPUR Package

The R implementation produces results consistent with the Stata SPUR
package (Becker, Boll & Voth, 2025). Key numerical checks:

| Test            | Stata article (§4.2)   | R (`spatialunitroot`)  |
|-----------------|------------------------|------------------------|
| I(0) conclusion | Reject (p = 0.0004)    | Reject (p \< 0.01)     |
| I(1) conclusion | Not reject (p = 0.574) | Not reject (p \> 0.05) |
| Half-life       | Very long / unbounded  | Unbounded above        |

Minor numerical differences arise from Monte Carlo simulation noise (±
1/√nrep) and differences in R vs. Stata random number generators, but
the qualitative conclusions are identical.

## References

- Becker, S. O., Boll, P. D., & Voth, H.-J. (2025). “Spatial Unit Roots
  in Regressions: A Practitioner’s Guide and a Stata Package.” *Stata
  Journal*, forthcoming.
- Chetty, R., Hendren, N., Kline, P., & Saez, E. (2014). “Where is the
  Land of Opportunity?” *Quarterly Journal of Economics*, 129(4),
  1553–1623.
- Müller, U. K. & Watson, M. W. (2024). “Spatial Unit Roots and Spurious
  Regression.” *Econometrica*, 92, 1661–1695.
