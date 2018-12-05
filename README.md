[![Build Status](https://travis-ci.org/RGLab/flowCore.svg?branch=trunk)](https://travis-ci.org/RGLab/flowCore)

NEWS
=====
`spillover_ng` has been moved to the `flowStats` package to eliminate a dependency
due needed for pregating. `spillover_match` has been added in `flowCore` to still provide the simplified
matching functionality.

We have a new spillover method `spillover_ng` that provides a simplified API. Matching
between single stained controls and channels is provided via a csv file that maps filenames to channels.


flowCore
========

Core flow infrastructure

Install the devtools package and then do
`devtools::install_github("RGLab/flowCore",ref="trunk")`

You may need to install other dependencies as well.
