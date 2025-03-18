gamlss2.ISD: Modeling of size spectra under the GAMLSS framework
================
Alexis Lazaris

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

The goal of gamlss2.ISD is to estimate the exponent of the Individual
Size Distribution for size spectra as described with a bounded Power-Law
distribution. The exponent can be estimated from raw body size data
within the GAMLSS framework as applied by the `gamlss2` package. All the
functionalities of the `gamlss2` package can be used.

## Installation

The `gamlss2.dist` package requires installation of `gamlss2`:

``` r
install.packages("gamlss2",
  repos = c("https://gamlss-dev.R-universe.dev",
            "https://cloud.R-project.org"))
```

It can then be installed using `devtools`:

``` r
devtools::install_github("alexislaz/gamlss2.ISD")
```

## Examples

First load relevant packages:

``` r
library(dplyr)
library(ggplot2)
library(gamlss2)
library(gamlss2.ISD)
```

## Estimate the ISD exponent from a single sample/population

We can simulate a dataset of body sizes from the bounded Power-Law
distribution using the `sizeSpectra` package:

``` r
# remotes::install_github("andrew-edwards/sizeSpectra")
d = data.frame(size = sizeSpectra::rPLB(1000, -2, 1, 1000))
```

And then estimate the exponent using GAMLSS:

``` r
m = gamlss2(size ~ 1, data = d, family = bPL(1, 1, 1000))

coef(m)
```

## Estimate the ISD exponent from multiple populations

As a more complicated example, we can simulate body sizes for different
populations with different characteristics:

``` r
pop1_size_range = c(10, 800)
pop2_size_range = c(20, 1000)
pop3_size_range = c(1, 300)

pop1_b = -1.2
pop2_b = -1.5
pop3_b = -1.8

d = rbind(
  data.frame(pop = "pop1", 
             size = round(sizeSpectra::rPLB(1e4, pop1_b, pop1_size_range[1], pop1_size_range[2])),
             mn = pop1_size_range[1], mx = pop1_size_range[2]),
  data.frame(pop = "pop2", 
             size = round(sizeSpectra::rPLB(1e4, pop2_b, pop2_size_range[1], pop2_size_range[2])),
             mn = pop2_size_range[1], mx = pop2_size_range[2]),
  data.frame(pop = "pop3", 
             size = round(sizeSpectra::rPLB(1e4, pop3_b, pop3_size_range[1], pop3_size_range[2])),
             mn = pop3_size_range[1], mx = pop3_size_range[2])
)
d$pop = factor(d$pop)
```

And the estimate of the exponents (b) of each population:

``` r
m = gamlss2(size ~ pop, data = d, family = bPL(1, d$mn, d$mx))

summary(m)$elapsed

predict(m, newdata = data.frame(pop = c("pop1", "pop2", "pop3")), se.fit = TRUE)
```

Alternatively, for faster fitting, we can aggregate equal sizes with
their respective abundance:

``` r
d_counts = d |> 
  group_by(pop, size, mn, mx) |>
  summarize(n = n()) |>
  ungroup()

m = gamlss2(size ~ pop, data = d_counts, family = bPL(d_counts$n, d_counts$mn, d_counts$mx))

summary(m)$elapsed

predict(m, newdata = data.frame(pop = c("pop1", "pop2", "pop3")), se.fit = TRUE)
```

We can, also, utilize the random effects structure that `gamlss2`
offers:

``` r
m = gamlss2(size ~ s(pop, bs = "re"), data = d_counts, family = bPL(d_counts$n, d_counts$mn, d_counts$mx))

predict(m, newdata = data.frame(pop = factor(levels(m$model$pop), levels(m$model$pop))), se.fit = TRUE)

plot(m)
```

## Estimate the effect of a predictor on the ISD exponent

Finally, we can add a predictor that affects the exponent of the ISD
while using random intercepts:

``` r
# setup global interncept and slope and a hardcoded 'random' intercept
set.seed(1)

intercept = -1.5   # global
slope = 0.1        # global

x.lvls = seq(-1, 1, length.out = 20)  # 'x' covariate
g.lvls = factor(c("pop1", "pop2"))    # different samples

# setup grid
x = rep(x.lvls, each = length(g.lvls))
g = rep(g.lvls, length(x.lvls))

b0 = c(0.15, -0.2) # 'random' intercept deviations from global
# the random intercepts are:
intercept + b0

X = model.matrix( ~ x)
Z0 = model.matrix( ~ -1 + g)

y = (X %*% c(intercept, slope)) + (Z0 %*% b0) + rnorm(length(x), 0, 0.01)

# collect in a data.frame
dat = data.frame(y = y, x = x, g = g)

# visualize the dependence of the exponent on 'x'
ggplot(dat) + 
  theme_classic(base_size = 15) + 
  geom_point(aes(x = x, y = y, colour = g)) + 
  geom_smooth(aes(x = x, y = y, colour = g), method = "lm", se = FALSE) +
  labs(x = "x", y = "ISD exponent", colour = "population")


# simulate sizes from the grid of exponents
dat_sizes = do.call(rbind,
                    Map(function(beta, x, g) 
                      data.frame(g = g, 
                                 x = x, 
                                 size = sizeSpectra::rPLB(500, beta, 1, 1000)), 
                      dat$y, dat$x, dat$g))

# visualize our raw data of individual sizes
# pop2 has a steeper size spectra "slope" (b) and, hence, more small individuals
ggplot(dat_sizes) + 
  theme_classic(base_size = 15) + 
  facet_wrap(vars(g)) + 
  geom_histogram(aes(x = size))
```

And fit the mixed model:

``` r
# run model
m = gamlss2(size ~ x + re(fixed = ~x, random = ~1|g), data = dat_sizes, family = bPL(1, 1, 1000))

summary(m)

plot(m)

coef(specials(m)$model)

# the estimated intercepts are:
coef(specials(m)$model)[, 1] + coef(m)[1]
# the actual intercepts are:
intercept + b0

# the estimated slopes are:
coef(specials(m)$model)[, 2] + coef(m)[2]
# the actual slope is:
slope
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-26dbcc1d-72d2-3940-bff0-891fe87c4e21" class="csl-entry">

Edwards, Andrew M., James P. W. Robinson, Julia L. Blanchard, Julia K.
Baum, and Michael J. Plank. 2020. “Accounting for the Bin Structure of
Data Removes Bias When Fitting Size Spectra.” *Marine Ecology Progress
Series* 636: pp. 19–33. <https://www.jstor.org/stable/26920653>.

</div>

<div id="ref-https://doi.org/10.1111/2041-210X.12641" class="csl-entry">

Edwards, Andrew M., James P. W. Robinson, Michael J. Plank, Julia K.
Baum, and Julia L. Blanchard. 2017. “Testing and Recommending Methods
for Fitting Size Spectra to Data.” *Methods in Ecology and Evolution* 8
(1): 57–67. https://doi.org/<https://doi.org/10.1111/2041-210X.12641>.

</div>

<div id="ref-10.1111/j.1467-9876.2005.00510.x" class="csl-entry">

Rigby, R. A., and D. M. Stasinopoulos. 2005. “Generalized Additive
Models for Location, Scale and Shape.” *Journal of the Royal Statistical
Society Series C: Applied Statistics* 54 (3): 507–54.
<https://doi.org/10.1111/j.1467-9876.2005.00510.x>.

</div>

<div id="ref-https://doi.org/10.1111/2041-210X.14312" class="csl-entry">

Wesner, Jeff S., Justin P. F. Pomeranz, James R. Junker, and Vojsava
Gjoni. 2024. “Bayesian Hierarchical Modelling of Size Spectra.” *Methods
in Ecology and Evolution* 15 (5): 856–67.
https://doi.org/<https://doi.org/10.1111/2041-210X.14312>.

</div>

</div>
