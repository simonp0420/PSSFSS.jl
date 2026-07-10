# PSSFSS - Analysis of polarization and frequency selective surfaces in Julia


| **Documentation**   |  **Tests**     | **CodeCov**  |
|:--------:|:---------------:|:-------:|
|[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simonp0420.github.io/PSSFSS.jl/stable/manual)  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simonp0420.github.io/PSSFSS.jl/dev/manual)| [![Build Status](https://github.com/simonp0420/PSSFSS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simonp0420/PSSFSS.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![codecov.io](https://codecov.io/github/simonp0420/PSSFSS.jl/coverage.svg?branch=main)](https://codecov.io/github/simonp0420/PSSFSS.jl?branch=main) |

`PSSFSS` is a Julia package (Python users should check out [pyPSSFSS](https://github.com/simonp0420/pypssfss)) for analyzing planar 
[polarization selective surfaces](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=polarization+selective+surface&btnG=) (PSSs), [frequency selective surfaces](https://en.wikipedia.org/wiki/Frequency_selective_surface) (FSSs), 
[reflectarray](https://en.wikipedia.org/wiki/Reflectarray_antennahttps://en.wikipedia.org/wiki/Reflectarray_antenna) elements, 
[radomes](https://en.wikipedia.org/wiki/Radome), and similar structures.  It is intended to be useful to antenna design engineers and others who work in applied electromagnetic engineering.

The user specifies the geometry to be analyzed as a `Vector` containing two or more dielectric `Layer`s 
and zero or more `Sheet` objects that define the PSS/FSS surfaces.  Due to the included plot recipes, the surfaces 
and their associated triangulations can be conveniently visualized using Julia's standard 
[`Plots`](https://github.com/JuliaPlots/Plots.jl) package. After also specifying the scan angles or
unit cell incremental phasings, frequencies to be analyzed, and optionally selecting performance parameters to be written
to [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) file(s), 
the user then invokes the `analyze` function to perform the analysis.  Post-processing and plotting of results can be performed in the same analysis script using the immensely powerful 
[Julia programming language](https://julialang.org/).


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
* Fast analysis for frequency sweeps using an extremely robust rational function interpolation algorithm.
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
* Analysis results can be exported to TICRA-compatible TEP (tabulated electrical properties) files.

## Limitations

* Only zero-thickness FSS/PSS sheets are currently supported.
* Frequency sweeps are fastest for normal incidence or for the case where unit cell incremental phase shifts ψ₁ and ψ₂ are
  constant with frequency (as in a waveguide).  This is due to the use of a wide-band expansion of the 
  potential Green's functions for a stratified medium with quasi-periodic excitation. Frequency sweeps for non-normal
  angle of incidence are typically slower. However, as of PSSFSS version 1.1, all frequency sweeps are now much faster, 
  often by more than an order of magnitude, compared to previous versions.  The speedup is due to the use of a fast interpolated sweep by default.

## News
- Version 1.1: A highly reliable fast sweep is now the default, resulting in order-of-magnitude speedups.
- Version 1.2: Sheet resistance has been replaced by a possibly complex sheet impedance.  Also, sheet
  conductivity and surface roughness may now be specified.
- Version 1.3: `loadedcross`, `jerusalemcross`, and 4-sided `polyring` elements are now triangulated 
  using a structured mesh, by default, resulting in reduced execution times.  To obtain the old, 
  unstructured mesh on these elements, specify `structuredtri = false` in the constructor argument list.
- Version 1.4: New, chiral `manji` element added.
- Version 1.5: `sinuous` element added.
- Version 1.6: Added `export_sheet` for exporting FSS/PSS sheets to STL-format CAD files.
- Version 1.7: Function `extract_result_file` is deprecated in favor of a new method for `extract_result`.
- Version 1.8: Improved multi-threading in matrix fill functions.
- Version 1.9: Added `orient` keyword to `rectstrip`, allowing geometry rotation within the stationary unit cell.
- Version 1.10: New `res2tep` function allows saving results in the form of a TICRA-compatible TEP 
  (tabulated electrical properties) file.
- Version 1.11: New `res2fresnel` function allows saving results in the form of an HFSS SBR+ Fresnel table.
- Version 1.12: Added `pixels` and `sympixels` elements.
- Version 1.13: Added `write_sheet_data` function.
  
## Installation
You can obtain and install PSSFSS using Julia's Pkg REPL-mode (hitting `]` as the first character at the command prompt):

```julia
(@v1.8) pkg> add PSSFSS
```

(and then hitting `<Backspace>` to return to the REPL) or with 
```julia
using Pkg: Pkg
Pkg.add("PSSFSS")
```


## Documentation
- The theory documentation is [here](https://github.com/simonp0420/PSSFSS.jl/blob/main/docs/TheoryDocs/theorydoc.pdf)
- The user manual is [here](https://simonp0420.github.io/PSSFSS.jl/stable/manual)

## Citing PSSFSS
If you use PSSFSS for a scientific publication, please cite the 2024
[ACES Journal paper](https://journals.riverpublishers.com/index.php/ACES/article/view/22443) in the following way:

> P. S. Simon, “PSSFSS—An Open-source Code for Analysis of Polarization and Frequency Selective Surfaces”, ACES Journal, vol. 39, no. 02, pp. 139–148, Feb. 2024.

or use this BibTeX entry:

```bib
@article{simo:24,
  title={{PSSFSS}---An Open-source Code for Analysis of Polarization and Frequency Selective Surfaces},
  volume={39},
  url={https://journals.riverpublishers.com/index.php/ACES/article/view/22443},
  DOI={https://doi.org/10.13052/2024.ACES.J.390207},
  number={02},
  journal={The Applied Computational Electromagnetics Society Journal (ACES)},
  author={Simon, Peter S.},
  year={2024},
  month={feb},
  pages={139--148}
}
```

## Community
Help from the community is actively sought and greatly appreciated!  There are several open issues which you might
want to tackle, and the documentation could always be improved. Pull requests are welcome.  Feel free to open more issues, whether for 
basic capability, performance, examples, documentation, etc.
