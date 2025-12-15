# Dynamic Pricing Engine: A 'From Scratch' Implementation 🍊
### [➡️ Read the Full Report](https://abdullahshahzadkhan.github.io/Dynamic-Pricing-MLE/)

## 📌 Project Overview
This project bridges the gap between theoretical economics and applied data science. Using scanner data from **Dominick’s Finer Foods**, I built a dynamic pricing algorithm to optimize the revenue of Tropicana Premium Orange Juice.

Unlike standard analyses that rely on "black-box" linear regression packages (`lm` or `glm`), this project builds the statistical engine from the ground up. I manually implemented **Maximum Likelihood Estimation (MLE)** to derive demand parameters and used calculus to solve for the mathematically optimal price point.

## 🎯 The Business Problem
* **Objective:** Determine the optimal price for Tropicana 64oz to maximize weekly store-level profit.
* **The "Trap":** A simple analysis suggested the product was performing well.
* **The Reality:** My model revealed the product was **underpriced by ~40%**, sacrificing significant margin for volume.
* **The Solution:** A corrected pricing strategy projected to **increase weekly profits by 81%**.

## 🚀 "Hard Mode" Methodology
This project demonstrates proficiency in statistical programming and mathematical derivation:

### 1. Manual Maximum Likelihood Estimation (MLE)
Instead of using pre-built libraries, I wrote a custom **Negative Log-Likelihood function** for the Gaussian demand curve and optimized it using the **L-BFGS-B** algorithm.
```r
# Snippet of the custom engine
nll_gaussian_multi <- function(theta, y, log_prices, log_price_comp) {
  y_hat <- theta[1] + theta[2]*log_prices + theta[3]*log_price_comp
  sum_prob_y <- sum(dnorm(y, log = TRUE, mean = y_hat, sd = theta[4]))
  return(-sum_prob_y)
}
