# ebayes_shrink

`ebayes_shrink` is a Stata command for applying Empirical Bayes shrinkage to one or more already-estimated fixed effects using per-observation variance-covariance estimates.

The command estimates the signal covariance matrix by method of moments and computes posterior mean shrunken effects. It can shrink fixed effects jointly, using the full covariance structure, or separately, using only diagonal variance entries.

## Installation

Install the package directly from GitHub with:

```stata
net install ebayes_shrink, from("https://raw.githubusercontent.com/martin-andresen/ebayes_shrink/main") replace
```

## Basic usage

Suppose `alpha` and `theta` are estimated fixed effects, and `S11`, `S21`, and `S22` contain the per-observation lower-triangle variance-covariance entries for `(alpha, theta)`. Then run:

```stata
ebayes_shrink alpha theta, wcov(S)
```

To shrink each effect separately using only the diagonal variance entries, run:

```stata
ebayes_shrink alpha theta, wcov(S) method(separate)
```

To choose custom prefixes for the generated shrunken effects and shrinkage weights:

```stata
ebayes_shrink alpha theta, wcov(S) shprefix(eb) wprefix(ww)
```

See the Stata help file for full syntax and options:

```stata
help ebayes_shrink
```

Citation

If you use `ebayes_shrink`, please cite:

Martin Eckhoff Andresen, Manudeep Bhuller, and Alfred Lovgren. *Pay Beliefs and the Amenity-Pay Tradeoff*.

License

This project is licensed under the MIT License. See `LICENSE` for details.
