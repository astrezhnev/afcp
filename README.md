# afcp

Implementation of the Average Feature Choice Probability (AFCP) estimator for conjoint experiments. 

For more details on the estimator see Abramson, Scott F., Korhan Kocak, Asya Magazinnik, and Anton Strezhnev. "Detecting Preference Cycles in Forced-Choice Conjoint Experiments." (2023). (https://osf.io/preprints/socarxiv/xjre9/)

# Installation

Install with `remotes::install_github`

First, install the `remotes` package if you do not have it already.

```{r}
install.packages("remotes")
```

Then install the development version directly from github using

```{r}
remotes::install_github("astrezhnev/afcp")
```

# Usage

See the built-in vignette `vignette("afcp")` for a guide on how to use the package.