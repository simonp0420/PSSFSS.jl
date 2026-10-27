# AMR implementation note

Initial design snapshot: 2026-10-02. This is a development companion to
[the proposed theory appendix](TheoryDocs/appendix_amr.tex). Mathematical
choices belong in that appendix; code integration and validation belong here.

## Status and scope

- **Implemented:** `tricharge!`, `tricharge2!`, and centroid `tricurrent!` in
  `src/AutoMeshRefine.jl`; output initialization uses `fill!` for immutable SVectors.
- **Validated:** the focused file last passed 35 assertions, including 34 new
  analytic assertions. This is historical targeted validation, not a full-suite
  pass or an implemented AMR feature.
- **Agreed:** centroid current dot products with testing vectors carrying length;
  electric- and magnetic-current sheets are both in initial scope; use the
  Green-function interface factors described below for both source classes.
- **Proposed:** all-retained-excitation maxima, area-weighted marking, neighbouring
  footprint, and a common mesh over selected adaptation conditions.
- **Open:** normalization, testing-vector identification, degenerate material
  handling, physical-edge handling, exact refinement algorithm, sampling, and stopping.

Do not silently promote proposals to requirements. Update this record after
mathematical agreement or implementation; record validation evidence separately.

## Proposed data flow

1. Reconstruct charge and centroid current for one solved excitation.
2. Build edge scores using neighbours in a common periodic image frame.
3. Average edge scores into cell indicators and apply a normalization yet to be agreed upon.
4. Accumulate a cellwise maximum across retained TE/TM modes on both Gblock sides
   and, when selected, across adaptation frequencies/scan conditions.
5. Apply the selected marking controller (area weighting is optional).
6. Expand the footprint, enforce periodic partner refinement, and form a conforming mesh.
7. Rebuild derived data, solve again, and compare scattering outputs.
8. Stop with an explicit reason: converged, iteration limit, unknown limit,
   minimum size, or inability to refine further.

Keep these stages separately callable and testable. Proposed APIs should not be
introduced until their inputs, units, and boundary behavior are settled.

## Integration points

- `src/PSSFSS.jl`: `calculate_jtype_gsm` / `calculate_mtype_gsm` solve loops.
  Immediately after `ldiv!`, `rwgdat.rhs` contains the current coefficients
  (amperes for electric, volts for magnetic sources). Compute indicators before
  the next right-hand side overwrites this buffer. Preserve GSM extraction.
- `src/RWG.jl`: use `bff`, `bfe`, `ebf`, and `eci` to build ordinary and periodic
  neighbour relations. Multiple contributions may connect the same face pair;
  retain edge identity and translation, not just neighbour face indices.
- `src/Sheets.jl`: preserve geometry/material metadata and the opposite-vertex
  relationship between `fv` and `fe`. Mesh changes invalidate cached integrals.
- `src/Meshsub.jl`, `src/structuredtri.jl`: candidate triangulation entry points.
  Preserve holes, physical boundaries, and paired periodic segmentation. Existing
  routines are not yet an adaptive remesher.
- `src/FillZY.jl`, `src/Zint.jl`: review all caches and reuse conditions after
  connectivity changes. Rebuild RWG data, matrices, Fourier storage, unique-face
  pair data, and mesh-dependent integral arrays as needed.
- `src/FastSweep.jl`: design adaptation and final interpolation stages explicitly.
  Reuse of GSMs and duplicate sheets/Gblocks must remain consistent after remeshing.  Since the same sheet object can occur in different Gblocks, a shared adaptive mesh needs indicators aggregated across all relevant occurrences before refinement, including occurrences with different surrounding layers. Duplicate-Gblock reuse must also propagate or reuse indicator data, not just the GSM.

## Electric and magnetic sources

Use the same reconstruction and geometry machinery for both source classes.
Document magnetic output units (V/m and Wb/m^2) alongside electric units.
The agreed interface weights use the immediately adjacent layers:

```text
mu_tilde    = 2 * mu_left * mu_right / (mu_left + mu_right)
epsilon_bar = (epsilon_left + epsilon_right) / 2

electric_score = abs(omega * mu_tilde * B_e) + abs(delta_q_e / epsilon_bar)
magnetic_score = abs(omega * epsilon_bar * B_m) + abs(delta_q_m / mu_tilde)
```

`B_e` / `B_m` are the centroid current/testing-vector contractions defined in
the appendix. Use absolute SI material parameters, including their complex
values for lossy media. Average the complex parameters first, then take the
magnitudes of the weighted contributions. Do not compute `k` or `eta` just to
form their products; the direct expressions avoid square-root branch choices.

Electric factors follow `eq:finalEGF` in Chapter 5. Magnetic factors follow the
sum of one-sided Green functions in Chapter 7 and their leading coefficients
in Chapter 5: arithmetic permittivity and arithmetic **reciprocal** permeability.
Thus magnetic charge uses `1/mu_tilde`, not `1/mu_bar`. A common omitted factor
of two for magnetic side sums must be treated consistently when fixing normalization.

`FillZY.jl` computes relative averages; restore `mu_0` and `epsilon_0` for SI
weights. Equal-material cases reduce to homogeneous weights. Zero/singular
effective parameters need a specified policy. Remote layers enter through the
solved response, not these local prefactors. The policy is agreed, but its
accuracy over material contrasts and losses still needs numerical validation.
Magnetic duality establishes the algebraic construction, not numerical
performance on all aperture geometries.

Existing tests independently validate the reconstruction algebra; they do not
yet exercise a magnetic GSM-to-indicator integration path. Add that validation
as solver integration develops.

## Verification plan

### Reconstruction and indicator tests

- Retain analytic unequal-area, SI conversion, dilation, inverse-frequency,
  zero-output, periodic half-basis phase, and uniform-current/zero-charge tests.
- Check proposed testing-vector identities on non-equilateral triangles.
- Check complex dot products, ordinary-edge orientation reversal, electric and
  magnetic dimensions, and a hand-computed nonzero CWCD example.
- Test interface averages with equal materials, equal permeabilities, unequal
  permeabilities, and complex lossy parameters. Check reciprocal averaging,
  SI factors, nonnegative scores, and the eventual degenerate-material policy.
- Test periodic image translation in both lattice directions without double
  applying phases; include a nonrectangular lattice.
- Specify and test zero/negligible normalization and physical boundary policies.
- Check streaming excitation aggregation against stored indicator maxima.

### Mesh and integration tests

- Conformity, positive areas, no lost/filled holes, preserved sheet area, edge/vertex
  connectivity, and matched periodic boundary segments after refinement.
- Cache invalidation even when scan phase and frequency are unchanged.
- Both J- and M-class solve paths, unequal source-side mode counts, evanescent
  modes, and repeated-sheet/Gblock reuse.
- A simple isolated sheet, then a closely spaced multilayer example; compare
  complex GSM and final S-parameters against converged uniform-mesh references.
- Check frequency/scan coverage and separately control Floquet mode truncation.
- Measure total time and unknown count, including repeated solves and loss of
  structured-mesh efficiencies. Fewer triangles alone does not prove a speedup.

Focused command, from the repository root:

```sh
julia --project=. -e 'include("test/AutoMeshRefine_test.jl")'
```

Full suite:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Use a runtime compatible with the current environment. Do not commit sandbox
loader hacks or report tests of an in-memory corrected implementation as tests
of the repository. Fresh-process validation complements persistent REPL testing.

## Decision record

| Date | Status | Decision or question |
| --- | --- | --- |
| 2026-10-01 | Implemented | Added analytic reconstruction tests and corrected SVector output initialization with `fill!`. |
| 2026-10-02 | Agreed | Start with centroid-based current testing; no silent change to exact path integration. |
| 2026-10-02 | Agreed | Electric and magnetic sheets both belong to initial AMR scope. |
| 2026-10-02 | Proposed | Per-excitation normalization and cellwise maximum over all retained incident modes. |
| 2026-10-02 | Proposed | Optional area-weighted marking; literature value 0.7 is not an accuracy tolerance. |
| 2026-10-02 | Open | Fix the exact normalization and negligible-response policy. |
| 2026-10-02 | Open | Confirm length-carrying testing vectors against the chosen centroid formulation. |
| 2026-10-02 | Agreed | Electric weights: `omega*mu_tilde` and `1/epsilon_bar`; magnetic weights: `omega*epsilon_bar` and `1/mu_tilde`, from adjacent-layer Green-function factors. |
| 2026-10-02 | Open | Resolve degenerate material handling and physical boundary contributions. |
| 2026-10-02 | Open | Choose refinement, sampling, output convergence norms, and resource limits. |

## Document integration

`docs/TheoryDocs/appendix_amr.tex` is self-contained for standalone review.
`theorydoc.tex` defines `AMRAppendixEmbedded` before inputting it, suppressing
its standalone preamble. It is inserted after existing appendices and before
back matter. It uses its own symbol definitions and source-reference list so
initial review does not require changes to the shared macro or bibliography files.
