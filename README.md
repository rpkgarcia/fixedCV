# fixedCV

Welcome! `fixedCV` is an R package that provides robust statistical inference for time series and other dependent data.

## What does it do?

When working with correlated data (like time series), standard statistical methods can give misleading results. This package provides tools for:

- **Robust hypothesis testing** that accounts for unknown correlation structures
- **Reliable confidence intervals** for regression coefficients in dependent data
- **Long-run variance estimation** using state-of-the-art kernel methods

The package implements fixed-b critical values and multiple kernel-based estimators to ensure your statistical inferences are valid even when observations are correlated.

## Installation

You can install `fixedCV` directly from GitHub using the `devtools` package:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install fixedCV from GitHub
devtools::install_github("rpkgarcia/fixedCV")
```

Alternatively, you can use the `remotes` package:

```r
# Install remotes if you haven't already
install.packages("remotes")

# Install fixedCV from GitHub
remotes::install_github("rpkgarcia/fixedCV")
```

## Quick Start

Once installed, load the package and you're ready to go:

```r
library(fixedCV)

# Fit a linear model
model <- lm(y ~ x, data = your_data)

# Get robust inference that accounts for autocorrelation
robust_results <- robust_lm(model)
```

The `robust_lm()` function automatically selects appropriate bandwidth parameters and provides robust standard errors, t-statistics, and p-values that are valid under general dependence structures.

## Getting Help

For detailed documentation and examples:

```r
# View package documentation
help(package = "fixedCV")

# See examples for the main function
?robust_lm

# Access the package vignette
vignette("fixedCV-vignette")
```

## Questions or Issues?

If you encounter any problems or have questions, please open an issue on the [GitHub repository](https://github.com/rpkgarcia/fixedCV/issues).
