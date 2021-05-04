# PSSFSS - analysis of polarization and frequency selective surfaces in Julia

**Author: Peter S. Simon (@simonp0420)**

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simonp0420.github.io/PSSFSS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simonp0420.github.io/PSSFSS.jl/dev)
[![Build Status](https://travis-ci.com/simonp0420/PSSFSS.jl.svg?branch=master)](https://travis-ci.com/simonp0420/PSSFSS.jl)
[![Coverage](https://codecov.io/gh/simonp0420/PSSFSS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/simonp0420/PSSFSS.jl)



`PSSFSS` is a Julia package for analyzing 
[polarization selective surfaces](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=polarization+selective+surface&btnG=) (PSSs), [frequency selective surfaces](https://en.wikipedia.org/wiki/Frequency_selective_surface) (FSSs), 
[reflectarray](https://en.wikipedia.org/wiki/Reflectarray_antennahttps://en.wikipedia.org/wiki/Reflectarray_antenna) elements, 
[radomes](https://en.wikipedia.org/wiki/Radome), and similar structures.  It is intended to be useful to antenna design engineers and others who work in applied electromagnetic engineering.

The user specifies the geometry to be analyzed as a `Vector` containing two or more dielectric [`Layer`](@ref)s 
and zero or more [`Sheet`](@ref) objects that define the PSS/FSS surfaces.  Due to the included plot recipes, the surfaces 
and their associated triangulations can be conveniently visualized using Julia's standard 
[`Plots`](https://github.com/JuliaPlots/Plots.jl) package. After also specifying the scan angles or
unit cell incremental phasings, frequencies to be analyzed, and optionally selecting performance parameters to be written
to [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) file(s), 
the user then invokes the [`analyze`](@refs) function to perform the analysis.  Post-processing and plotting of results can be
performed in the same analysis script using the immensely powerful Julia programming language.


## Features

* Designed to be useful and accessible to working engineers.
* Accommodates planar FSS/PSS surfaces with no limits to number of dielectric layers or FSS/PSS sheets.
* Automatically chooses number of modes needed for cascading multiple FSS/PSS sheets using
  generalized scattering matrices (GSMs).
* Supports (approximate) cascading multiple sheets of different periodicities, as in a multilayer
  meanderline polarizer.
* Simple specification of geometry to be analyzed.
* Solution of mixed-potential integral equation using Rao-Wilton-Glisson triangle subdomain basis functions 
  and multi-threaded method of moments.
* Fast analysis for frequency sweeps without approximations or interpolation using a wide-band expansion of the 
  potential Green's functions for a stratified medium with quasi-periodic excitation.
* Automatic triangulation of sheet geometries to user-specified number of triangles.
* Exploits redundancies inherent in structured meshes for greater numerical efficiency.
* Easy extraction of useful engineering performance parameters, including 
    * Reflection and transmission coefficient magnitudes and/or phases or complex coefficients for the field components of 
        * TE/TM 
        * Vertical/horizontal (Ludwig 3)
        * LHCP/RHCP (circular polarization)
    * Delta insertion phase delay (ΔIPD)
    * Delta insertion loss (ΔIL)
    * Axial ratio 

## Limitations

* Only zero-thickness FSS/PSS sheets are currently supported.
* Frequency sweeps are fast for normal incidence or for the case where unit cell incremental phase shifts ψ₁ and ψ₂ are constant
  with frequency (as in a waveguide).  Frequency sweeps where the angle of incidence is held constant are typically much slower (except for normal incidence).


## Installation
PSSFSS is not yet registered, so it must currently be installed by the Julia package manager with an explicit URL:

```Julia
julia> ]
(v1.6) pkg> add add https://github.com/simonp0420/PSSFSS.jl
```

## Documentation
- The theory documentation is [here](https://github.com/simonp0420/PSSFSS.jl/blob/main/docs/TheoryDocs/theorydoc.pdf)
- The user manual is [here](https://simonp0420.github.io/PSSFSS.jl/stable)

## Community
Help from the community is actively sought and greatly appreciated!  There are several open issues which you might
want to tackle, and the documentation could always be improved. Pull requests are welcome.  Feel free to open more issues, whether for 
basic capability, performance, examples, documentation, etc.
