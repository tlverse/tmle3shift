
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tmle3shift`

[![Travis-CI Build
Status](https://travis-ci.org/tlverse/tmle3shift.svg?branch=master)](https://travis-ci.org/tlverse/tmle3shift)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/tmle3shift?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/tmle3shift)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/tmle3shift/master.svg)](https://codecov.io/github/tlverse/tmle3shift?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

> Targeted Learning and Variable Importance with Stochastic
> Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Jeremy
Coyle](https://github.com/jeremyrcoyle), and [Mark van der
Laan](https://vanderlaan-lab.org)

-----

## What’s `tmle3shift`?

`tmle3shift` is an adapter/extension R package in the `tlverse`
ecosystem that exposes support for the estimation of target parameters
defined as shifts of interventions of continuous-valued treatments by
way of stochastic treatment regimes. As an adapter package, `tmle3shift`
builds upon the core `tlverse` grammar introduced by `tmle3`, a general
framework that supports the implementation of a range of TMLE parameters
through a unified interface. For a detailed description of the target
parameter, TML estimator, and algorithm implemented in `tmle3shift`, the
interested reader is invited to consult
(<span class="citeproc-not-found" data-reference-id="munoz2012population">**???**</span>)
and Díaz and van der Laan (2018).

For a general discussion of the framework of targeted minimum loss-based
estimation and the role this methodology plays in statistical and causal
inference, the canonical references are van der Laan and Rose (2011) and
van der Laan and Rose (2018).

`vimshift` is an R package in the `tlverse` ecosystem that implements a
strategy for assessing variable importance measures (VIMs) based on a
parameter defined as the mean counterfactual outcome under a specified
set of shifts of a given set of continuous-valued variables of interest.
The parameter at the heart of the proposed VIM is defined through
stochastic treatment regimes. For a detailed description of the target
parameter, as well as the corresponding TML estimator, the interested
reader is invited to consult Díaz and van der Laan (2012) and Díaz and
van der Laan (2018) (n.b., the `tlverse` implementation of the relevant
TML estimator is available throught the [`tmle3shift` R
package](https://tmle3shift.tlverse.org)). The `vimshift` package builds
upon the core `tlverse` grammar introduced by `tmle3`, a general
framework that supports the implementation of a range of TMLE parameters
through a unified interface, building upon this framework to introduce
an end-to-end methodology for variable importance.

For a general discussion of the framework of targeted minimum loss-based
estimation and the role this methodology plays in statistical and causal
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

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/tmle3shift/issues).

-----

## Related

  - [R/`txshift`](https://github.com/nhejazi/txshift) - An R package
    providing an independent implementation of the TML estimation
    procedure and statistical methodology as is made available here,
    without reliance on the `tlverse` grammar provided by `tmle3`.

-----

## Funding

The development of this software was supported in part through a grant
from the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

-----

## License

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-diaz2012population">

Díaz, Iván, and Mark J van der Laan. 2012. “Population Intervention
Causal Effects Based on Stochastic Interventions.” *Biometrics* 68 (2).
Wiley Online Library: 541–49.

</div>

<div id="ref-diaz2018stochastic">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

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
