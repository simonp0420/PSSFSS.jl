```@meta
CurrentModule = PSSFSS
```

# [PSSFSS](https://github.com/simonp0420/PSSFSS) - analysis of polarization and frequency selective surfaces in Julia


[PSSFSS](https://github.com/simonp0420/PSSFSS) is a Julia package for analyzing 
[polarization selective surfaces](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=polarization+selective+surface&btnG=) (PSSs), [frequency selective surfaces](https://en.wikipedia.org/wiki/Frequency_selective_surface) (FSSs), 
[reflectarray](https://en.wikipedia.org/wiki/Reflectarray_antennahttps://en.wikipedia.org/wiki/Reflectarray_antenna) elements, 
[radomes](https://en.wikipedia.org/wiki/Radome), and similar structures.  It is intended to be useful to antenna design engineers and others who work in applied electromagnetic engineering.

The user specifies the geometry to be analyzed as a `Vector` containing two or more dielectric `Layer`s 
and zero or more `RWGSheet` objects that define the PSS/FSS surfaces.  Due to the included plot recipes, the surfaces 
and their associated triangulations can be conveniently visualized using Julia's standard 
[`Plots`](https://github.com/JuliaPlots/Plots.jl) package. After also specifying the scan angles or
unit cell incremental phasings, frequencies to be analyzed, and optionally selecting performance parameters to be written
to [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) file(s), 
the user then invokes the `analyze` function to perform the analysis.  Post-processing and plotting of results can be
performed in the same analysis script using the immensely powerful Julia programming language.


## Features

* Designed to be useful and accessible to working engineers.
* Accommodates planar FSS/PSS surfaces with no limits to number of dielectric layers or FSS/PSS sheets.
* Simple specification of geometry to be analyzed.
* Automatically chooses number of modes needed for rigorously cascading multiple FSS/PSS sheets of identical
  periodicities using generalized scattering matrices (GSMs).
* Supports (approximate) cascading multiple sheets of different periodicities, as in a multilayer
  meanderline polarizer.
* Solution of mixed-potential integral equation using Rao-Wilton-Glisson triangle subdomain basis
  functions and multi-threaded method of moments.
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

* Only zero-thickness, planar FSS/PSS sheets are currently supported.
* Frequency sweeps are fast for normal incidence or for the case where unit cell 
  incremental phase shifts *ψ₁* and *ψ₂* are  constant with frequency (as in a waveguide).
  Frequency sweeps where the angle of incidence is held constant with frequency are generally 
  much slower (except for normal incidence).

## Installation
You can obtain PSSFSS using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):

```julia
(v1.6) pkg> add PSSFSS
```

or with `using Pkg; Pkg.add("PSSFSS")`.


## Documentation
- The theory documentation is [here](https://github.com/simonp0420/PSSFSS.jl/blob/main/docs/TheoryDocs/theorydoc.pdf)
- The user manual is [here](https://simonp0420.github.io/PSSFSS.jl/stable)
- If you prefer to explore the documentation more interactively, Jupyter notebooks containing documentation can be found
  [here](https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/notebooks) 


## Community
Help from the community is actively sought and greatly appreciated!  There are several open issues which you might
want to tackle, and the documentation could always be improved. Pull requests are welcome.  Feel free to open more issues, whether for 
basic capability, performance, examples, documentation, etc.
