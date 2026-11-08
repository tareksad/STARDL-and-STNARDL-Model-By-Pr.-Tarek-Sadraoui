# pystardl

**Smooth Transition ARDL (STARDL) and Smooth Transition Nonlinear ARDL (STNARDL) for Python**

![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg) ![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg) ![Version](https://img.shields.io/badge/version-0.1.0-green.svg)

Estimation, inference and diagnostics for ARDL error-correction models whose
coefficients change **smoothly** between regimes (STARDL), optionally combined with
**asymmetric** positive/negative partial-sum effects (STNARDL).

**Author:** Pr. Tarek Sadraoui  
**Email:** Tarek.sadraoui@gmail.com  
**GitHub:** https://github.com/tareksad

---

## Table of Contents
1. [Installation](#installation)
2. [Features](#features)
3. [The models](#the-models)
4. [Quick start](#quick-start)
5. [Estimators](#estimators)
6. [Tests and diagnostics](#tests-and-diagnostics)
7. [Long run, multipliers, plots](#long-run-multipliers-plots)
8. [API reference](#api-reference)
9. [Limitations](#limitations)
10. [References](#references)

## Installation
```
pip install pystardl            # once published
# or from source
git clone https://github.com/tareksad/pystardl.git
cd pystardl && pip install -e .
```
Dependencies: numpy, pandas, scipy, tabulate, matplotlib (pytest for the test-suite).

## Features
- **STARDL** and **STNARDL** with logistic (LSTR1), exponential (ESTR) and second-order logistic (LSTR2) transitions
- Transition variable: time trend, any column, or its first difference, with any delay
- Choose which coefficients move: `nonlinear='all' | 'lr' | 'sr' | [names]`
- Grid search + Nelder-Mead refinement; full NLS standard errors (nonrobust, HC1, HAC)
- Linearity test (Terasvirta LM-type, third-order Taylor) with LSTR1 / ESTR-LSTR2 suggestion
- Bounds test with simulated critical values for the exact specification
- Long-run coefficients per regime with delta-method SE; time-varying long-run path
- Long-run and short-run symmetry Wald tests per regime (STNARDL)
- Regime-dependent cumulative dynamic multipliers (m+, m-, asymmetry)
- Lag selection (AIC/BIC/HQIC) and comparison of transition shapes
- Diagnostics: Breusch-Godfrey, Breusch-Pagan, ARCH, Jarque-Bera, RESET, Durbin-Watson, CUSUM / CUSUMSQ
- Stata-style tables (`tabulate`), LaTeX / HTML export, figures

## The models
ECM form, with `G_t = G(s_t; gamma, c) in [0, 1]`:

```
dy_t = a + rho*y_{t-1} + theta'x_{t-1} + sum_i phi_i dy_{t-i} + sum_i pi_i' dx_{t-i}
       + G_t * [ a* + rho*·y_{t-1} + theta*'x_{t-1} + sum_i phi*_i dy_{t-i} + sum_i pi*_i' dx_{t-i} ] + e_t
```

Transition functions: `logistic 1/(1+exp(-gamma(s-c)))`, `exponential 1-exp(-gamma(s-c)^2)`,
`logistic2 1/(1+exp(-gamma(s-c1)(s-c2)))`.

**STNARDL** replaces `x` by `x+ = cumsum(max(dx,0))` and `x- = cumsum(min(dx,0))`
(Shin, Yu and Greenwood-Nimmo, 2014). The long run in regime `g` is
`LR_j(g) = -(theta_j + g*theta*_j) / (rho + g*rho*)`.

p and q are ARDL orders: the ECM has `p-1` lags of `dy` and `q` terms `dx_t ... dx_{t-q+1}` (q = 0 means none).

## Quick start
```python
from pystardl import STARDL, STNARDL, linearity_test, select_lags, simulate_stardl

df = simulate_stardl(300)                                  # or your DataFrame in time order
print(linearity_test(df, "y", "x", candidates=["time", "D.x"]))

p, q, table = select_lags(df, "y", "x", max_p=3, max_q=2, criterion="bic")
res = STARDL(df, "y", "x", p=p, q=q, transition_var="time", transition="logistic").fit()
print(res)                                   # coefficients, gamma, c, long run, speed of adjustment
print(res.bounds_test(n_sim=400))
print(res.diagnostics())
```

## Estimators
```python
from pystardl import STARDL, STNARDL

# STARDL: smooth regime change in every coefficient, transition driven by an external variable
STARDL(data, depvar="gdp", indepvars=["oil", "fdi"], p=2, q=2,
       transition_var="D.oil", delay=1, transition="exponential").fit()

# only the long-run terms move
STARDL(data, "gdp", ["oil"], nonlinear="lr").fit()

# only the coefficients of selected variables move (all their lags)
STARDL(data, "gdp", ["oil", "fdi"], nonlinear=["oil"]).fit()

# STNARDL: oil asymmetric, fdi symmetric
STNARDL(data, "gdp", ["oil", "fdi"], asym_vars=["oil"], p=1, q=1,
        transition_var="time", cov_type="HAC").fit()
```
`model.column_names` lists the regressors, i.e. the names accepted in `nonlinear=[...]`.
Pass `time_col="year"` if the frame is not already ordered.

## Tests and diagnostics
| Method | Purpose |
|---|---|
| `linearity_test(...)` | LM-type test of linearity against smooth transition; picks `s`, suggests LSTR1 vs ESTR/LSTR2 |
| `compare_transitions(...)` | fits all three shapes, ranks by BIC |
| `res.bounds_test(n_sim)` | F and t bounds test, simulated I(0)/I(1) critical values |
| `res.asymmetry_tests()` | long-run / additive short-run symmetry, per regime and jointly |
| `res.wald_test(R, r)`, `res.test_coefficients({...})` | arbitrary linear restrictions |
| `res.diagnostics()` | serial correlation, heteroskedasticity, ARCH, normality, RESET, DW |
| `res.cusum()` / `res.plot_cusum()` | parameter stability |

## Long run, multipliers, plots
```python
res.long_run(g=(0, 0.5, 1))                 # -theta/rho by regime, delta-method SE
res.long_run_path()                         # LR coefficients along G_t
res.dynamic_multipliers("x", horizon=20)    # plus / minus / asym per regime
res.plot_transition(); res.plot_regime_path(); res.plot_long_run()
res.plot_multipliers("x"); res.plot_cusum(); res.plot_fit()
res.to_dataframe(); res.to_latex(); res.to_html()
```
Multiplier convention: `plus` and `minus` are the responses of y to a unit change in `x+` and `x-`;
`asym = plus - minus` is zero under symmetry.

## API reference
**Constructor arguments** (STARDL / STNARDL): `data, depvar, indepvars, p=1, q=1, transition_var='time',
delay=0, transition='logistic', nonlinear='all', constant=True, trend=False, time_col=None,
cov_type='nonrobust', hac_lags=None, n_gamma=9, n_c=15, gamma_max=300` (+ `asym_vars` for STNARDL).

**Result attributes:** `coefficients, std_errors, transition_params, gamma, c, gamma_se, c_se, G, residuals,
fitted, r2, adj_r2, aic, bic, hqic, loglik, n, K, cov_full, index, names`.

## Limitations
- `gamma` and `c` are not identified under linearity, so Wald tests on the G* terms (`test_nonlinear_terms`) are descriptive only; use `linearity_test`.
- Bounds critical values are simulated conditional on the estimated transition path; they are specification-specific, not a replacement for published tables.
- Standard errors treat `gamma` and `c` as estimated by NLS; when `gamma` is very large (abrupt break) they are unreliable.
- Time-series models only; a panel extension is not included.

## References
- Terasvirta, T. (1994). Specification, estimation and evaluation of smooth transition autoregressive models. *JASA*, 89, 208-218.
- Luukkonen, R., Saikkonen, P., Terasvirta, T. (1988). Testing linearity against smooth transition autoregressive models. *Biometrika*, 75, 491-499.
- Pesaran, M.H., Shin, Y., Smith, R.J. (2001). Bounds testing approaches to the analysis of level relationships. *Journal of Applied Econometrics*, 16, 289-326.
- Shin, Y., Yu, B., Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework. In *Festschrift in Honor of Peter Schmidt*, Springer, 281-314.
- van Dijk, D., Terasvirta, T., Franses, P.H. (2002). Smooth transition autoregressive models - a survey of recent developments. *Econometric Reviews*, 21, 1-47.
- Davies, R.B. (1977). Hypothesis testing when a nuisance parameter is present only under the alternative. *Biometrika*, 64, 247-254.
- Brown, R.L., Durbin, J., Evans, J.M. (1975). Techniques for testing the constancy of regression relationships over time. *JRSS B*, 37, 149-192.

## Citation
```
@software{pystardl,
  author = {Sadraoui, Tarek},
  title  = {pystardl: Smooth Transition ARDL and NARDL for Python},
  year   = {2026},
  url    = {https://github.com/tareksad/pystardl}
}
```
## License
MIT
