# UnfoldRIDE

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldRIDE.jl/dev)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/workflows/Test/badge.svg)](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions)
[![Test workflow status](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/unfoldtoolbox/UnfoldRIDE.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldRIDE.jl)
<!--[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)-->
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/unfoldtoolbox/UnfoldRIDE.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

A re-implementation of the [RIDE](https://cns.hkbu.edu.hk/RIDE.htm) algorithm in Julia with an extension to replace the RIDEs iterative decomposition with an [Unfold](https://github.com/unfoldtoolbox/Unfold.jl) deconvolution.

## Install

### Installing Julia

<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).
It allows you to, e.g., easily update Julia at a later point, but also test out alpha/beta versions etc.

TL:DR; If you dont want to read the explicit instructions, just copy the following command

#### Windows

AppStore -> JuliaUp,  or `winget install julia -s msstore` in CMD

#### Mac & Linux

`curl -fsSL https://install.julialang.org | sh` in any shell
</details>

### Installing UnfoldRIDE

```julia
using Pkg
Pkg.add("UnfoldRIDE")
```

## Quickstart

```Julia
#config for ride algorithm
cfg = RideConfig(
    #sfreq is the sampling frequency of the data
    sfreq = 100,
    #ranges for the individual components have to be determined through manual inspection of the data
    s_range = [-0.1, 0.3],
    r_range = [0, 0.4],
    c_range = [-0.4, 0.4],
    #the range in which the initial peak estimation for the C component is performed
    c_estimation_range = [0, 0.9],
    #the range for one epoch
    epoch_range = [-0.1, 1]
)
#run the ride algorithm
resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
resultsUnfold = ride_algorithm(UnfoldMode, data_noisy, evts_without_c, cfg)
```

## How to Cite

If you use UnfoldRIDE.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/blob/main/CITATION.cff).

## Contributing
Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

### How-to Contribute

You are very welcome to raise issues and start pull requests!

### Adding Documentation

1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
