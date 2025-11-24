# Define the stationary SPDE cgeneric model for INLA.

Define the stationary SPDE cgeneric model for INLA.

## Usage

``` r
cgeneric_sspde(mesh, alpha, control.priors, constr = FALSE, ...)
```

## Arguments

- mesh:

  triangulation mesh to discretize the model.

- alpha:

  integer used to compute the smoothness parameter.

- control.priors:

  named list with parameter priors. This shall contain `prange` and
  `psigma` each one as a length two vector with (U, a) to define the
  PC-prior parameters such that P(range\<U)=a and P(sigma\>U)=a,
  respectively. See Fuglstad et. al. (2019) \<DOI:
  10.1080/01621459.2017.1415907\>. If a=0 or a=NA, then U is taken to be
  the fixed value of the parameter.

- constr:

  logical, default is FALSE, to indicate if the integral of the field
  over the domain is to be constrained to zero.

- ...:

  additional arguments that will be passed on to
  [`INLAtools::cgenericBuilder()`](https://rdrr.io/pkg/INLAtools/man/cgeneric-class.html),
  such as: `debug` : logical/integer, default is FALSE/0.
  `useINLAprecomp` logical, default is TRUE, indicating if it is to be
  used the shared object pre-compiled by INLA.

## Value

objects to be used in the f() formula term in INLA.

## Note

This is the stationary case of
[`INLA::inla.spde2.pcmatern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.pcmatern.html)
with slight change on the marginal variance when the domain is the
sphere, following Eq. (23) in Lindgren et. al. (2024).

## References

Geir-Arne Fuglstad, Daniel Simpson, Finn Lindgren & Håvard Rue (2019).
Constructing Priors that Penalize the Complexity of Gaussian Random
Fields. Journal of the American Statistical Association, V. 114, Issue
525.

Finn Lindgren, Haakon Bakka, David Bolin, Elias Krainski and Håvard Rue
(2024). A diffusion-based spatio-temporal extension of Gaussian Matérn
fields. [SORT vol. 48, no. 1, pp.
3-66](https://raco.cat/index.php/SORT/article/view/428665) \<doi:
10.57645/20.8080.02.13\>
