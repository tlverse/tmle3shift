
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tmle3shiftshift`

[![Travis-CI Build
Status](https://travis-ci.org/tlverse/tmle3shift.svg?branch=master)](https://travis-ci.org/tlverse/tmle3shift)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/tmle3shift?branch=master&svg=true)](https://ci.appveyor.com/project/tlverse/tmle3shift)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/tmle3shift/master.svg)](https://codecov.io/github/tlverse/tmle3shift?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

> TML Estimation of Shift Intervention Effects with the `tmle3`
> framework

-----

## What’s `tmle3shift`?

`tmle3shift` is an adapter / plug-in R package for use with `tmle3`, a
general framework that supports the implementation of a range of TMLE
parameters with a unified interface. `tmle3shift` extends support to the
estimation of parameters defined as shift interventions of
continuous-valued treatments by way of stochastic treatment regimes. For
a detailed description of the parameter, TML estimator, and algorithm
implemented in `tmle3shift`, the interested reader is invited to consult
Muñoz and van der Laan (2012) and Díaz and van der Laan (2018).

For a general discussion of the framework of targeted minimum loss-based
estimation and the role this methodology plays in statistical causal
inference, the canonical references are van der Laan and Rose (2011) and
van der Laan and Rose (2018).

-----

## Installation

You can install the development version of `tmle3shift` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with

``` r
devtools::install_github("tlverse/tmle3shift")
```

-----

## Getting Started

The best place to get started is the “Framework Overview” document,
which describes the individual components of the `tmle3shift` framework.
It may be found here:
<https://tmle3shift.tlverse.org/articles/framework.html>.

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/tmle3shift/issues).

-----

## License

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-diaz2018stochastic">

Díaz, Iván, and Mark J van der Laan. 2018. “Stochastic Treatment
Regimes.” In *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*, 167–80. Springer Science & Business
Media.

</div>

<div id="ref-munoz2012population">

Muñoz, Iván Díaz, and Mark J van der Laan. 2012. “Population
Intervention Causal Effects Based on Stochastic Interventions.”
*Biometrics* 68 (2). Wiley Online Library: 541–49.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vdl2018targeted">

———. 2018. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.

</div>

</div>
