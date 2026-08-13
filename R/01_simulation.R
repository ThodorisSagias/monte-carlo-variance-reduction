# We will start by defining the parameters that we will need in the procedure

S0 <- 100 # Initial Stock Price
K <- 100 # Strike Price
B <- 85 # Barrier Level
r <- 0.05 # Risk-free Rate
s <- 0.2 # Volatility
TT <- 1.0 # Time to Maturity in Years
M <- 252 # Discrete Time Steps (Trading Days in a year)
N <- 100000 # Number of Simulated Paths
dt <- TT / M

set.seed(888)  

# Standard Monte Carlo & Euler-Maruyama Engine
# Generate N x M matrix with standard normal variables

Z <- matrix(rnorm(N * M), nrow = N, ncol = M)

# Calculate discrete market shocks with using Euler-Maruyama

shocks <- exp((r - 0.5 * s^2) * dt + s * sqrt(dt) * Z)

# Building the asset paths

paths <- cbind(S0, S0 * t(apply(shocks, 1, cumprod)))

# Find minimum price of each path in order to check if we fell below the barrier or not

min_S <- apply(paths, 1, min)

# Standard Monte Carlo Payoff (It will take 0 as value if it is below the barrier)

payoffMC <- exp(-r * TT) * pmax(paths[, M+1] - K, 0) * (min_S > B)


###############################################################

# Antithetic Variates by using -Z
# Shocks using the mirrored matrix this time -Z

shocksMirrored <- exp((r - 0.5 * s^2) * dt + s * sqrt(dt) * (-Z))
pathsMirrored <- cbind(S0, S0 * t(apply(shocksMirrored, 1, cumprod)))
min_SMirrored <- apply(pathsMirrored, 1, min)
payoffMirrored <- exp(-r * TT) * pmax(pathsMirrored[, M+1] - K, 0) * (min_SMirrored > B)

# Antithetic Estimator which is the average of the standard and mirrored paths

payoffAV <- (payoffMC + payoffMirrored) / 2

###############################################################

# Control Variates 
# Calculate European Call payoff using the same standard paths that we created above

payoffEC <- exp(-r * TT) * pmax(paths[, M+1] - K, 0)

# Exact Analytical Price of European Call using Black-Scholes

d1 <- (log(S0 / K) + (r + 0.5 * s^2) * TT) / (s * sqrt(TT))
d2 <- d1 - s * sqrt(TT)

bs_price <- S0 * pnorm(d1) - K * exp(-r * TT) * pnorm(d2)

# Calculation of the optimal control coefficient

c_star <- cov(payoffMC, payoffEC) / var(payoffEC)

# Control Variate Estimator

payoffCV <- payoffMC - c_star * (payoffEC - bs_price)

# Results and Variance Calculation
# Average Payoff for each approach

mean(payoffMC)
mean(payoffAV)
mean(payoffCV)
var(payoffMC)

# Calculate percentage of variance decrease

ReductionAV <- (1 - (var(payoffAV) / var(payoffMC))) * 100
ReductionCV <- (1 - (var(payoffCV) / var(payoffMC))) * 100
print(c(mean(payoffMC), var(payoffMC)))
print(c(mean(payoffAV), var(payoffAV), ReductionAV))
print(c(mean(payoffCV), var(payoffCV), ReductionCV))