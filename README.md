
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
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4603372.svg)](https://doi.org/10.5281/zenodo.4603372)

> Targeted Learning of the Causal Effects of Stochastic Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Jeremy
Coyle](https://github.com/jeremyrcoyle), and [Mark van der
Laan](https://vanderlaan-lab.org)

-----

## What’s `tmle3shift`?

`tmle3shift` is an adapter/extension R package in the `tlverse`
ecosystem that exposes support for the estimation of a target parameter
defined as the mean counterfactual outcome under a posited shift of the
natural value of a continuous-valued intervention, using the formalism
of stochastic treatment regimes. As an adapter package, `tmle3shift`
builds upon the core `tlverse` grammar introduced by `tmle3`, a general
framework that supports the implementation of a range of TMLE parameters
through a unified interface. For a detailed description of the target
parameter, TML estimator, and algorithm implemented in `tmle3shift`, the
interested reader is invited to consult Dı́az and van der Laan (2012)
and Dı́az and van der Laan (2018). For a general discussion of the
framework of targeted minimum loss-based estimation and the role this
methodology plays in statistical and causal inference, the canonical
references are van der Laan and Rose (2011) and van der Laan and Rose
(2018).

Building on the original work surrounding the TML estimator for the
aforementioned target parameter, `tmle3shift` additionally implements a
set of techniques for variable importance analysis, allowing for a
sequence of mean counterfactual outcomes, estimated under a sequence of
posited shifts, to be summarized via a working marginal structural model
(MSM). The goal of this work is to build upon the `tlverse` framework
and the estimation methodology implemented for a single mean
counterfactual outcome in order to introduce an end-to-end methodology
for variable importance analyses.

-----

## Installation

You can install the development version of `tmle3shift` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes) with

``` r
remotes::install_github("tlverse/tmle3shift")
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/tmle3shift/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/tlverse/tmle3shift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `tmle3shift` R package, please cite the following:

``` 
    @software{hejazi2021tmle3shift-rpkg,
      author = {Hejazi, Nima S and Coyle, Jeremy R and {van der Laan}, Mark
        J},
      title = {{tmle3shift}: {Targeted Learning} of the Causal Effects of
        Stochastic Interventions},
      year = {2021},
      howpublished = {\url{https://github.com/tlverse/tmle3shift}},
      note = {{R} package version 0.2.0},
      url = {https://doi.org/10.5281/zenodo.4603372},
      doi = {10.5281/zenodo.4603372}
    }
```

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

Dı́az, Iván, and Mark J van der Laan. 2012. “Population Intervention
Causal Effects Based on Stochastic Interventions.” *Biometrics* 68 (2):
541–49.

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
