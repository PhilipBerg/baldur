<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

![R-CMD-check](https://github.com/PhilipBerg/baldur/actions/workflows/check-standard.yaml/badge.svg)
![Rhub](https://github.com/PhilipBerg/baldur/actions/workflows/rhub.yaml/badge.svg)
![DOI Badge](https://img.shields.io/badge/10.1016%2Fj.mcpro.2023.100658-green?style=palstic&label=DOI%3A&labelColor=grey&link=https%3A%2F%2Fdoi.org%2F10.1016%2Fj.mcpro.2023.100658)
[![CRAN status](https://www.r-pkg.org/badges/version/baldur)](https://CRAN.R-project.org/package=baldur)
<!-- badges: end -->

# Baldur

Baldur is a hierarchical Bayesian model for the analysis of proteomics data. By leveraging empirical Bayes methods, Baldur estimates hyperparameters for variance and measurement-specific uncertainty. It then computes the posterior difference in means between conditions for each peptide, protein, or PTM, and integrates the posterior to estimate error probabilities.

## Features

- **Hierarchical Bayesian modeling** of proteomics data
- **Empirical Bayes estimation** of variance and uncertainty
- **Posterior probability calculations** for differential analysis
- Supports **peptide, protein, and PTM** data

## Installation

Install the stable release from CRAN:

```r
install.packages('baldur')
```

Or, install the development version from GitHub (after installing `rstan`):

Follow the instructions for installing `rstan`:
<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

Then:

```r
devtools::install_github('PhilipBerg/baldur', build_vignettes = TRUE)
```

> **Note:**  
> - On Ubuntu, [pandoc](https://pandoc.org/) may be needed to build vignettes.  
> - On Windows, sometimes the development version of `rstan` is required.

## Usage

For detailed examples, see the package vignettes:

```r
vignette('baldur_yeast_tutorial')
vignette('baldur_ups_tutorial')
```

## Main Modeling Work and Equations

**Baldur** implements a hierarchical Bayesian framework for label-free proteomics quantification, designed to robustly estimate differential abundance while accounting for the mean-variance relationship in mass spectrometry data.
For exact details please see the original paper.

### 1. Observation Model

For each feature (peptide, protein, or PTM) $i$ in sample $j$, the observed intensity $y_{ij}$ is modeled as:

$$y_{ij}\sim\text{Normal}(\mu_{j},\sigma u_{ij})$$
$$\mu_{j}\sim\text{Normal}(\mu_{0j}+\eta_j\sigma,\sigma)$$

### 2. Mean-Variance Modeling

#### Gamma Regression

The measurement standard deviation $s_{j}$ is not constant, but depends on the mean intensity. This relationship is modeled with gamma regression:
$$s_{j} \sim \Gamma(\alpha, \frac{\alpha}{\beta(\bar{y}_j)})$$
where:
- $\alpha$: shape parameter (estimated empirically)
- $\beta(\bar{y}_j)$: rate parameter as a function of peptide/protein mean intensity

#### Latent Gamma Mixture Regression

##### Regression Function

For each observation, the expected mean-variance relationship is modeled as:

$$\beta_i=\kappa\cdot\exp(\theta_i\cdot(I_L-S_Lx_i))+\exp(I-S\bar{y}_i)$$

where:
- $S, S_L$: slope parameters (common and latent)
- $I, I_L$: intercepts (common and latent)
- $\bar{y}_i$: mean
- $\theta_i$: feature-specific mixture parameter

##### Likelihood

Given the expected mean-variance, the observed standard deviation $\sigma_i$ is modeled as:
$$\sigma_i\sim\Gamma(\alpha,\frac{\alpha}{\beta_i})$$
where:
- $\alpha$: gamma shape parameter
- $\beta_i$: expected mean-variance for observation $i$

##### Priors

- $\alpha\sim\text{Cauchy}(0,25)$
- $\eta\sim\text{Normal}(0,1)$
- $I_L\sim\text{SkewNormal}(2,15,35)$
- $\theta_i\sim\text{Uniform}(0,1)$

##### NRMSE (Model Fit Metric)

The normalized root-mean-square error (NRMSE) is calculated for model diagnostics.

### 3. Hierarchical Modeling of Condition Means

#### **Empirical Bayes Prior**

- **What is it?**  
  A prior that is *estimated from your actual data*, rather than set by hand.
- **How does it work?**  
  - Baldur looks at the spread and center of observed means across features.
  - It sets the mean hyper-prior for each group to match the average observed mean, and its uncertainty to match the variability in your data.
- **Why use it?**  
  - *Strength:* Provides "shrinkage" toward realistic values, reducing noise and false positives, especially with small sample sizes.
- **Mathematical form:**  
  $\mu_{0j}\sim\text{Normal}(\bar{y}_j,\sigma n_R)$
  - Here, $\bar{y}_j$ is the estimated mean for group $j$, and $\sigma n_R$ is the estimated standard deviation for the prior.

#### **Weakly Informative Prior**

- **What is it?**  
  A broad, generic prior that doesn't make strong assumptions—it's like saying "I have no idea what the mean should be, but it's probably not infinite."
- **How does it work?**  
  - The prior mean is set to zero for all groups.
  - The uncertainty (standard deviation) is set to a large value (e.g., $10$), meaning the model expects almost any value is possible.
- **Why use it?**  
  - *Strength:* Maximizes flexibility and lets the data speak for itself, at the cost of potentially adding more noise or less stability if data is limited.
- **Mathematical form:**  
  $$\mu_{0j}\sim\text{Normal}(0,10)$$
  - For group $j$, the mean is $0$ and the standard deviation is $10$.

### 4. Differential Abundance

For differential analysis, Baldur estimates the posterior distribution of the difference in means between conditions:
$$\boldsymbol{D}\sim\mathcal{N}(\boldsymbol{\mu}^\text{T}\boldsymbol{K},\sigma\boldsymbol{\xi}),\quad \xi_{m}=\sqrt{\sum_{i=1}^{C}\frac{|k_{im}|}{n_i}}$$

where:
- $\boldsymbol{K}$: contrast matrix
- $k_{im}$: contrast coefficient for condition $i$ in contrast $m$
- $n_i$: number of samples in condition $i$
- $\boldsymbol{\xi}$: scaling factor for each contrast

The probability of error for contrast $c$ is then:

$$P(\mathrm{error}) = 2\Phi(-|\mu_{D_c} - \mu_{h_0}| \odot \tau_{D_c})$$

where:
- $\Phi$: cumulative distribution function (CDF) of the standard normal
- $\mu_{h_0}$: null hypothesis mean (often zero)
- $\boldsymbol{\tau}_{\boldsymbol{D}}$: precision (inverse standard deviation) for each contrast
- $\odot$: element-wise multiplication

---

**Summary:**  
Baldur combines hierarchical modeling, mean-variance trend estimation via gamma regression, and empirical Bayes to robustly quantify differential abundance and propagate uncertainty from individual measurements to protein/PTM level, outputting interpretable error probabilities for each feature.

For full details, see the reference publication.

## Reference

Berg, Philip, and George Popescu.  
“Baldur: Bayesian Hierarchical Modeling for Label-Free Proteomics with Gamma Regressing Mean-Variance Trends.”  
_Molecular & Cellular Proteomics_ (2023): 2023-12.  
[https://doi.org/10.1016/j.mcpro.2023.100658](https://doi.org/10.1016/j.mcpro.2023.100658)
