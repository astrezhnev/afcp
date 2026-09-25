# afcp

Implementation of the Average Feature Choice Probability (AFCP) estimator for conjoint experiments. 

For more details on the estimator see Abramson, Scott F., Korhan Kocak, Asya Magazinnik, and Anton Strezhnev. "Aggregation, Interpretation, and Estimation of Preferences in Conjoint Experiments." (2026). (https://osf.io/preprints/socarxiv/xjre9_v3)

# Installation

Install with `remotes::install_github`

First, install the `remotes` package if you do not have it already.

```{r}
install.packages("remotes")
```

Then install the development version directly from github using

```{r}
remotes::install_github("astrezhnev/afcp", build_vignettes = TRUE)
```

# Usage

See the built-in vignette `vignette("afcp")` for a guide on how to use the package.