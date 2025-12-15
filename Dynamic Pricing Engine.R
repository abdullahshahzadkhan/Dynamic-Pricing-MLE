# Load Necessary Libraries
library("bayesm")
library("dplyr")
library(ggplot2)

# Load Data
data(orangeJuice)

str(orangeJuice)

oj_data <- orangeJuice$yx
head(oj_data)
str(oj_data)

# Filter for Brand: Tropicana
tropicana <- oj_data %>% 
  filter(brand == 1) %>% 
  mutate(log_price = log(price1))

# Visualize the Relation Between Price and Quantity Sold
plot(tropicana$log_price, tropicana$logmove)

ggplot(tropicana, aes(x = log_price, y = logmove)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(title = "Tropicana Demand Curve", x = "Log Price", y = "Log Quantity")

# Benchmark Linear Regression Model
benchmark_model <- lm(logmove ~ log_price, tropicana)
summary(benchmark_model)
coef(benchmark_model)
# For every 1% increase in price, Tropicana loses 2.71% of its volume

# Manual Maximum Likelihood Estimation (MLE)
# μi=β0+β1*xi
# Negative Log Likelihood Function. Minimizing the Negative Likelihood is the same as maximizing the Likelihood.
nll_gaussian <- function(theta, y, x) {
  # Predicted Value
  y_hat <- theta[1] + theta[2]*x
  
  #Probability Density of actual y
  sum_prob_y <- sum(dnorm(y, log = TRUE, mean = y_hat, sd = theta[3]))
  negative_sum_prob_y <- (-1*sum_prob_y)
  
  return(negative_sum_prob_y)
}

# Optimization
start_vals <- c(0, 0, 1) # Intial vals for intercept, slope and sigma

mle_results <- optim(par = start_vals, 
            method = "BFGS", 
            fn = nll_gaussian,
            y = tropicana$logmove,
            x = tropicana$log_price)
print(mle_results$par)

# Compare with the benchmark model
cat("Manual MLE Slope:", mle_results$par[2], "\n")
cat("Benchmark Slope:", benchmark_model$coefficients[2])
# Both coefficients are very close: Manual MLE Slope: -2.651496, Benchmark Slope: -2.711666

# Fisher Information (Hessian Matrix)
# SE(θ)=(diagonal(Inverse Hessian))^-2
mle_results <- optim(par = start_vals, 
                     method = "BFGS", 
                     fn = nll_gaussian,
                     y = tropicana$logmove,
                     x = tropicana$log_price,
                     hessian = TRUE)
mle_results$hessian

# Invert the hessian matrix to get the covariance matrix
cov_matrix <- solve(mle_results$hessian)
variance_vec <- diag(cov_matrix) # Variance
variance_vec
std_error_vec <- sqrt(variance_vec)
std_error_vec
# Problems: Negative variance, SE too high, Extremely small hessian values

# Fix: Use an advanced optimization function i.e. numDeriv
library("numDeriv")

# Redefining function with different variable names as Hessian function has an x parameter already
nll_gaussian <- function(theta, y, log_prices) {
  # Predicted Value
  y_hat <- theta[1] + theta[2]*log_prices
  
  #Probability Density of actual y
  sum_prob_y <- sum(dnorm(y, log = TRUE, mean = y_hat, sd = theta[3]))
  negative_sum_prob_y <- (-1*sum_prob_y)
  
  return(negative_sum_prob_y)
}

# Hessian Matrix
hessian_mat <- numDeriv::hessian(func = nll_gaussian,
                  x = mle_results$par,
                  y = tropicana$logmove,
                  log_prices = tropicana$log_price)

# Invert the matrix and get SE
fisher_info <- solve(hessian_mat)
std_errors <- sqrt(diag(fisher_info))
std_errors
# Result is still the same. Optim needs further fixing.

# Taking initial guesses from our benchmark model so that MLE doesnt get lost
theta_start <- c(coef(benchmark_model)[1],  # Intercept from lm
                 coef(benchmark_model)[2],  # Slope from lm
                 1)

# Optim with more restrictions
mle_results <- optim(par = theta_start, 
                     fn = nll_gaussian, 
                     y = tropicana$logmove, 
                     log_prices = tropicana$log_price,
                     method = "L-BFGS-B",
                     lower = c(-Inf, -Inf, 0.0001), # Sigma must be > 0
                     control = list(reltol = 1e-12, maxit = 10000))

hessian_mat <- hessian(func = nll_gaussian, 
                       x = mle_results$par, 
                       y = tropicana$logmove, 
                       log_prices = tropicana$log_price)

fisher_info <- solve(hessian_mat)
std_errors <- sqrt(diag(fisher_info))

# Results
results_table <- data.frame(
  Parameter = c("Intercept", "Price_Elasticity", "Sigma"),
  Estimate = mle_results$par,
  Std_Error = std_errors
)
print(results_table)
cat("\nBenchmark SE for Slope:", summary(benchmark_model)$coefficients[2,2])

# Revenue Maximization (R = P * Q)
# Demand Equation: ln(Q)=β0+β1ln(P) ---> Q = e^B0 * P^B1
# If we maximize R here by substituting the value for Q, we would get: R=A/P^1.7
# This says that the best move would be to reduce the price to zero
# Solution: Instead of maximizing revenue, we will maximize profit
# Profit(Π)=(Price−Cost)×Quantity
# P* = Cost x (B/1+B)
marg_cost <- 0.04 # Assumed marginal cost

beta <- mle_results$par[2]

opt_price <- marg_cost*(beta/(1+beta))
opt_price # What the price should be: 

mean(tropicana$price1) # What the price is

cat("Current Price:", mean(tropicana$price1),
    "\nOptimal Price:", opt_price,
    "\nDifference:", opt_price-mean(tropicana$price1))

# Conclusion
# Actual Price: 0.04485128 
# Optimal Price: 0.06336905 
# Tropicana is underpriced by nearly 40%
# If they raised price,  they would lose some volume but their total profit will increase due to higher per unit price

# Business Impact Analysis (How much are they losing at current price?)
b0 <- mle_results$par[1]
b1 <- mle_results$par[2]
cost <- 0.04

q_current <- exp(b0) * (mean(tropicana$price1))^b1
q_optimal <- exp(b0) * (opt_price)^b1

profit_current <- (mean(tropicana$price1) - cost) * q_current
profit_optimal <- (opt_price - cost) * q_optimal

profit_current
profit_optimal

# Change in Quantity Sold
cat("Current Quantity:", q_current,
    "\nOptimal Quantity:", q_optimal,
    "\n%-change:", ((q_optimal-q_current)/q_current) * 100)

# Change in Profit
cat("Current Profit:", profit_current,
    "\nOptimal Profit:", profit_optimal,
    "\n%-change:", ((profit_optimal-profit_current)/profit_current) * 100)

# But what happens if a competitor (minute maid) changes their price?
# Cross Price Elasticity
tropicana$log_price_mm <- log(tropicana$price5)

# Visualize the relationship between competitor price and quantitiy sold
plot(tropicana$log_price_mm, tropicana$logmove)

ggplot(tropicana, aes(x = log_price_mm, y = logmove)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Cross-Elasticity: Tropicana Sales vs Minute Maid Price",
       x = "Log Price (Minute Maid)",
       y = "Log Quantity (Tropicana)")
# We see a downward trend which is weird so we will run regression on both our price and competitor price and 
# analyse the coefficients when everything else is held constant

# Benchmark Model
benchmark_model_comp <- lm(logmove ~ log_price + log_price_mm, data = tropicana)
summary(benchmark_model_comp)
# The coeff of price of mm: 0.49429 suggests that for every 1 unit increase in price of mm, quantity sold of tropicana increases by 0.49%

# Multivariate MLE
# μi=β0+β1⋅ln(Pown)+β2⋅ln(Pcompetitor)
# Update our fucntion
nll_gaussian <- function(theta, y, log_prices, log_price_comp) {
  # Predicted Value
  y_hat <- theta[1] + theta[2]*log_prices + theta[3]*log_price_comp
  
  #Probability Density of actual y
  sum_prob_y <- sum(dnorm(y, log = TRUE, mean = y_hat, sd = theta[4]))
  negative_sum_prob_y <- (-1*sum_prob_y)
  
  return(negative_sum_prob_y)
}

# Update Optim
theta_start <- c(coef(benchmark_model_comp)[1],
                 coef(benchmark_model_comp)[2],
                 coef(benchmark_model_comp)[3],
                 1)

mult_mle_results <- optim(par = theta_start, 
                     fn = nll_gaussian, 
                     y = tropicana$logmove, 
                     log_prices = tropicana$log_price,
                     log_price_comp = tropicana$log_price_mm,
                     method = "L-BFGS-B",
                     lower = c(-Inf, -Inf, -Inf, 0.0001), # Sigma must be > 0
                     control = list(reltol = 1e-12, maxit = 10000))
mult_mle_results$par
# Coeffs match with those of benchmark model

# Key insight: Our own elasticity changed from -2.71 to -2.83 when we took competitor's price into account. 
# The univariate model was not considering price changes of minute maid. Assume tropicana raised their prices, minute maid did so too.
# Our model saw a drop in sales but not by much since competitor also raised their price so the actual impact reduced.

# Calculate the Optimal Price
beta <- mult_mle_results$par[2]

opt_price_mult <- marg_cost*(beta/(1+beta))
cat("Optimal Price Before:", opt_price,
    "Optimal Price Now:", opt_price_mult)
# Optimal price is now lower than before since customers are more price sensitive.