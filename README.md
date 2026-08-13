# Variance Reduction for Monte-Carlo Barrier Option Pricing

Pricing a down-and-out barrier call by Monte Carlo simulation, with two variance
reduction techniques implemented from scratch in R.

**Control variates cut estimator variance by 97.02%. Antithetic variates cut it
by 73.06%.** Both estimators remain unbiased, with all three prices agreeing to
within 0.10 of each other.


## The Problem

Under the risk-neutral measure, an option's fair value is a discounted
expectation. Discretely-monitored barrier options have no closed form, so that
expectation has to be estimated by simulation and plain Monte Carlo converges
at O(1/√N), meaning a halving of the error costs four times the paths. Scaling N
is the expensive answer. Reducing the variance of the estimator is the cheap one.

## Approach

Asset paths follow geometric Brownian motion, discretised by Euler-Maruyama over
M = 252 steps to mirror daily barrier monitoring. Parameters: S₀ = K = 100,
B = 85, r = 0.05, σ = 0.20, T = 1, N = 100.000 paths, seed 888.

**Antithetic variates.** Each shock matrix Z is paired with its reflection −Z.
Since −W is also a standard Brownian motion, the mirrored paths are valid draws,
and payoffs monotone in W give the pair a negative covariance which is
precisely the condition under which the paired average has lower variance.

**Control variates.** A European call with the same strike and maturity is highly
correlated with the barrier payoff and has a known Black-Scholes price, giving an
exact expectation to correct against. The optimal coefficient
c\* = Cov(X, Y) / Var(Y) is estimated from the simulated sample, and minimises
variance to Var(X)(1 − ρ²).

All three estimators are driven by the same normal shock matrix, so the
measured variance differences are attributable to the techniques rather than to
sampling noise.

## Results

| Method | Price | Variance | Variance reduction | Gain at equal compute |
| Standard Monte Carlo | 9.9510 | 216.7461 | — | — |
| Antithetic variates | 10.0101 | 58.3868 | 73.06% | 46.1% |
| Control variates | 10.0457 | 6.4691 | 97.02% | 97.0% |


The last column matters. An antithetic *sample* consumes two simulated paths, so
its per-sample variance drop overstates the benefit at a fixed computational
budget — on equal compute the honest figure is about 46%, not 73%. Control
variates reuse the paths already generated, so their cost is essentially
unchanged and the 97% survives the correction. Chart 2 above plots the antithetic
estimator against 2n paths for this reason.

Control variates dominate here largely because the barrier at 85 sits well below
the strike at 100: the two payoffs coincide in almost every scenario that pays
out at all, so ρ is close to 1.

## Reproduce

```r
# from the repository root
source("R/01_simulation.R")   
```

Runs in base R.

## Limitations

- **Discrete monitoring bias.** Checking the barrier at 252 points rather than
  continuously misses intra-day breaches, biasing the price upward. A
  Brownian-bridge correction would address this and is not implemented.
- **c\* estimated in-sample.** The same paths estimate the control coefficient
  and the price, introducing a small bias that vanishes as N → ∞. A pilot run on
  a held-out subset would remove it.
- **Constant volatility.** No local or stochastic volatility; no jumps.
- **Correlation-dependent gains.** Move the barrier toward the strike and ρ falls,
  taking the control variate's advantage with it. The 97% is specific to this
  parameter set, not a general property of the method.
- **Antithetic gains require monotone payoffs**, which the barrier indicator
  weakens; the 73% would shrink further for payoffs less monotone in W.

## Repository

```
R/01_simulation.R   simulation, three estimators, results table
docs/report.pdf     full write-up with derivations
results/            summary CSV
```

Full theoretical treatment, including the variance proofs for both methods, is in
[docs/report.pdf](docs/report.pdf).
