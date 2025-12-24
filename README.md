# Generalized-Linear-Model-GLM
A statistical pipeline for Seed Length analysis using Generalized Linear Models (GLM) with a Gamma distribution. Features automated outlier detection, orthogonal contrasts, AIC/BIC model evaluation, and multi-factor visualization of QTL effects


# Seed Length Analysis — Generalized Linear Model (GLM) Pipeline

This repository contains a robust analytical pipeline in R to investigate the effects of QTL presence, planting conditions (Normal vs. Late), and biomass levels on **Seed Length** in wheat.

## 📊 Statistical Framework: Why GLM?
Standard linear regression assumes a normal distribution. However, biological traits like seed length are strictly positive and often exhibit a non-normal distribution. 

This pipeline uses a **Gamma GLM with a Log Link** to model the data:
$$E(Y) = \mu$$
$$\log(\mu) = \beta_0 + \beta_1 X_1 + \dots + \beta_k X_k$$

* **Gamma Distribution:** Chosen because seed length data is continuous and strictly positive.
* **Log Link:** Ensures that predicted values remain positive and allows for multiplicative effects between factors.

## 🚀 Key Features
* **Automated Outlier Detection:** Implements Tukey’s Rule (1.5 * IQR) to detect and report outliers across every Year × Planting × QTL combination.
* **Orthogonal Contrasts:** Custom contrast matrices are defined for categorical variables (Year, Planting, QTL Status) to allow for specific, independent hypothesis testing.
* **Model Evaluation:** Provides **AIC** (Akaike Information Criterion) and **BIC** values to assess model fit, alongside standard diagnostic plots.
* **Automated Result Export:** Filters for significant coefficients (p < 0.05) and automatically saves them to a structured CSV file.

## 📈 Visualization Suite
The pipeline generates high-quality `ggplot2` visualizations, including:
* **EMM Plots:** Visualizes Estimated Marginal Means (EMMs) for complex interactions.
* **Distribution Analysis:** Combined Violin and Boxplots to show raw data density and quartiles faceted by Year and Planting.
* **Significance Comparison:** Color-coded plots highlighting the difference between "WITH QTL" and "WITHOUT QTL" groups.

## 🛠 Tech Stack
* **Language:** R
* **Modeling:** `glm()` (Gamma family), `emmeans`
* **Data Processing:** `dplyr`, `janitor`
* **Visualization:** `ggplot2`

## 📂 Project Structure
* `Seed_Length_GLM.R` - The primary analysis script.
* `/Results/` - Directory for significant coefficient tables and diagnostic PNGs.
