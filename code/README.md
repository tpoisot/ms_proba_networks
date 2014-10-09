# About the code

The code is in `julia`

I will package it as a module later

## A quick overview of the files

- `proba_utils` has functions to calculate the mean and variance, as well
as the variance of cumulative additive or multiplicative events -- at the
moment it is necessary to `include` it in each file, but that will not be
the case once packaged in a module.

- `matrix_utils` has (at the moment) the code to make bipartite networks
unipartites

- `degree` has all sorts of code related to the degree distribution.

