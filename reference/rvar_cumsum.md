# Cumulative sum for posterior::rvar objects

Computes the running sum of a
[posterior::rvar](https://mc-stan.org/posterior/reference/rvar.html)
vector, preserving the full posterior distribution of each partial sum.

## Usage

``` r
rvar_cumsum(rvars)
```

## Arguments

- rvars:

  A [posterior::rvar](https://mc-stan.org/posterior/reference/rvar.html)
  vector.

## Value

An rvar vector of the same length holding the cumulative sums.
