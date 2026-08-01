# Simulated truth ratings of true and fake news headlines

A simulated trial-level dataset in a classic signal-detection design,
inspired by studies of psychological inoculation against misinformation:
participants in two conditions (control / inoculation) rated true and
fake news headlines before and after the intervention on a 6-point
perceived-truth scale.

## Usage

``` r
sdt_ratings
```

## Format

A data frame with 4,000 rows and 5 variables:

- id:

  Participant identifier (factor, 50 levels).

- condition:

  Between-participant condition: `control` or `inoculation`.

- time:

  Measurement occasion: `pre` or `post` intervention.

- target:

  Type of headline: `fake` or `true`.

- value:

  Perceived-truth rating on an ordered 1–6 scale (1 = "definitely fake",
  6 = "definitely true").

## Details

The data were generated from an unequal-variance probit SDT model (fake
headlines centered at 0 with unit SD; true headlines shifted by d' with
SD 1.2), with a participant random intercept. In the simulation,
inoculation increases post-intervention discrimination (d') and makes
the response criteria stricter. See `data-raw/sdt_ratings.R` in the
source repository for the exact generating code.

## Examples

``` r
fit <- ordinal::clm(value ~ target * time, data = sdt_ratings, link = "probit")
summary(fit)
#> formula: value ~ target * time
#> data:    sdt_ratings
#> 
#>  link   threshold nobs logLik   AIC      niter max.grad cond.H 
#>  probit flexible  4000 -6415.74 12847.48 5(0)  8.51e-11 7.0e+01
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> targettrue           1.19291    0.04852  24.584   <2e-16 ***
#> timepost            -0.09002    0.04631  -1.944   0.0519 .  
#> targettrue:timepost  0.11911    0.06640   1.794   0.0729 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Threshold coefficients:
#>     Estimate Std. Error z value
#> 1|2 -0.99508    0.03926 -25.343
#> 2|3 -0.26616    0.03532  -7.535
#> 3|4  0.39475    0.03531  11.180
#> 4|5  1.04480    0.03751  27.857
#> 5|6  1.66810    0.04126  40.431
```
