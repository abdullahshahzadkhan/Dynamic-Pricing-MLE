# Dynamic Pricing Engine: Algorithmic Revenue Optimization
### [➡️ Read the Full Report](https://abdullahshahzadkhan.github.io/Dynamic-Pricing-MLE/)

## 📌 Project Overview
This project implements a **dynamic pricing engine** to optimize retail revenue, using scanner data from the Dominick’s Finer Foods dataset.

The core objective was to move beyond standard regression packages (`lm`) to build a transparent, white-box statistical model. By manually implementing **Maximum Likelihood Estimation (MLE)** and **Numerical Optimization**, this project derives precise demand elasticities and calculates mathematically optimal price points to maximize store-level profitability.

## 💼 The Business Case
* **Context:** Retail CPG products often face "Omitted Variable Bias" where competitor pricing disguises true customer price sensitivity.
* **Problem:** Initial analysis suggested the product was performing optimally.
* **Diagnosis:** After controlling for competitive dynamics, the model revealed the product was **underpriced by ~40%**.
* **Outcome:** The algorithm proposes a corrected pricing strategy projected to **increase weekly profits by 81%** by capturing margin expansion.

## ⚙️ Technical Methodology

### 1. Custom Likelihood Estimation
To ensure full control over the statistical assumptions, a custom **Negative Log-Likelihood** function was developed for the demand curve. This allows for flexible extension to non-normal distributions (e.g., Poisson for count data) in future iterations.
* **Optimization:** Parameters were estimated using the **L-BFGS-B** algorithm via R's `optim()` function.

### 2. Statistical Inference
Standard errors were not taken from a package but derived from first principles. The **Fisher Information Matrix** (inverse Hessian) was calculated using the `numDeriv` library to ensure coefficient estimates were statistically significant.

### 3. Price Optimization (Calculus)
The optimal price point ($P^*$) was derived analytically by solving the first-order condition of the profit function:

$$P^* = \text{Cost} \times \left( \frac{\beta}{1 + \beta} \right)$$

## 📊 Impact Analysis

| Metric | Current Status | Optimized Model | Change |
| :--- | :--- | :--- | :--- |
| **Unit Price** | $0.045 / oz | $0.062 / oz | **+40%** |
| **Sales Volume** | 8,632 units | 3,476 units | **-60%** |
| **Weekly Profit** | $41.88 | $75.94 | **+81.3%** |

## 🛠️ Technologies
* **Language:** R
* **Core Libraries:** `bayesm`, `numDeriv`, `optim`, `dplyr`
* **Techniques:** Maximum Likelihood Estimation, Matrix Algebra, Econometric Modeling.

---
