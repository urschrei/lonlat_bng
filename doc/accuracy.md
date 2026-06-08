# Accuracy and precision of `lonlat_bng`

This note describes the numeric precision and the measured accuracy of the
coordinate transformations implemented by `lonlat_bng`. It is intended for
surveyors, geodesists, and engineers who need to state the accuracy of the
library in a technical report.

## Summary

`lonlat_bng` converts between ETRS89/WGS84 geographic coordinates and OSGB36
British National Grid eastings and northings using the Ordnance Survey OSTN15
transformation. The transformation has two parts: a Transverse Mercator
projection on the GRS80 ellipsoid, and the OSTN15 grid-shift correction. Two
implementations of the projection step are provided:

- **Default** — the truncated series published by Ordnance Survey ("Redfearn"
  series). It reproduces OS-published results and is the method used by other OS
  software such as Grid InQuest II.
- **`karney_tm` (optional Cargo feature)** — Karney's Krüger *n*-series, accurate
  to a few nanometres over the whole grid. It removes the series-truncation error
  of the default method at the extreme west of the grid.

Measured against the 40-point OSTN15/OSGM15 developer-pack reference:

| Direction | Default (Redfearn) | `karney_tm` |
|:----------|:-------------------|:------------|
| Forward (lon/lat → OSGB36), max | 1.4 mm (≤ 1 mm per axis) | 1.0 mm (≤ 1 mm per axis) |
| Round trip (lon/lat → OSGB36 → lon/lat), mean | 0.42 mm | 0.29 mm |
| Round trip, max | **4.97 mm** (St Kilda, TP31) | **0.61 mm** (TP31) |

Both methods are sub-millimetre across mainland Great Britain. They differ only
at the far-western isles, where the default method's truncation error is largest.

## 1. The transformation pipeline

Forward (geographic to grid):

```
ETRS89 lon/lat
  → Transverse Mercator (GRS80 ellipsoid, National Grid parameters)
  → ETRS89 grid easting/northing
  → OSTN15 bilinear grid shift
  → OSGB36 easting/northing
```

Reverse (grid to geographic):

```
OSGB36 easting/northing
  → iterative reverse OSTN15 grid shift
  → ETRS89 grid easting/northing
  → inverse Transverse Mercator
  → ETRS89 lon/lat
```

National Grid projection parameters: central meridian 2°W, latitude of true
origin 49°N, scale factor on the central meridian F0 = 0.9996012717, false
easting 400 000 m, false northing −100 000 m, GRS80 ellipsoid
(a = 6 378 137.0 m, b = 6 356 752.3141 m).

**OSTN15 is the definitive realisation of OSGB36.** The shift grid is applied
exactly by bilinear interpolation, so the datum transformation itself is not an
approximation. The only approximations the library introduces are:

1. the Transverse Mercator projection series (the subject of this note);
2. IEEE-754 double-precision (binary64) floating-point arithmetic; and
3. rounding of grid coordinates to the millimetre on output.

## 2. Numeric precision

All internal computation uses `f64` (IEEE-754 binary64, approximately 15–16
significant decimal digits). Intermediate coordinates are carried at full
precision; in particular, the ETRS89 grid coordinate passed between the
projection and the OSTN15 step is not rounded.

- Grid eastings and northings are rounded to the nearest millimetre (half away
  from zero) on output.
- Longitude and latitude are returned at full `f64` precision (not rounded).

Iteration tolerances:

- The inverse projection's footpoint latitude is solved by iteration to a
  tolerance of 1 × 10⁻⁵ m.
- The reverse OSTN15 step iterates until successive grid shifts agree to within
  1 × 10⁻⁴ m (0.1 mm), to a maximum of ten iterations; coordinates that do not
  converge return an error.

## 3. The two projection implementations

### 3.1 Default: truncated OS (Redfearn) series

This is the series given in the Ordnance Survey *A guide to coordinate systems in
Great Britain* (Annex C) and in the *OSGM15 transformation and user guide*
(Annex B). The forward projection is a power series in the longitude difference
from the central meridian, truncated after the sixth power; the inverse is a
power series in the easting difference from the false origin, truncated after the
seventh power.

This is the method Ordnance Survey specifies, and it reproduces the output of OS
software. Its truncation error grows steeply with distance from the central
meridian (2°W). Because the points furthest from 2°W within Great Britain are the
north-western Scottish isles, the error appears to be latitude-dependent, but the
controlling variable is the east–west distance from the central meridian.

### 3.2 Optional: Karney Krüger *n*-series (`karney_tm`)

Enabled with the `karney_tm` Cargo feature. This implements the Transverse
Mercator projection as a series in the third flattening *n*, carried to order
*n*⁶, following Karney (2011). It is accurate to a few nanometres across the
whole British National Grid extent, with no truncation growth toward the grid
edges. The implementation is self-contained (no external dependency) and only
replaces the GRS80 projection legs; the OSTN15 step is unchanged.

Build and test with the feature:

```
cargo build --features karney_tm
cargo test  --features karney_tm
```

## 4. Measured accuracy

**Reference data.** Ordnance Survey publishes a 40-point test dataset with the
OSTN15/OSGM15 developer pack, spanning the whole grid from the Isles of Scilly to
Shetland. The published coordinates are produced by the OS authoritative
computation and are used here as ground truth. Errors below are horizontal ground
distances in millimetres.

### 4.1 Forward: lon/lat → OSGB36

| | Default (Redfearn) | `karney_tm` |
|:--|:--|:--|
| mean | 0.14 mm | 0.18 mm |
| max  | 1.41 mm | 1.00 mm |
| max per axis | 1 mm | 1 mm |

Both methods agree with the published OSGB36 coordinates to the last published
digit. The residual is the millimetre rounding of the published values, not a
projection error.

### 4.2 Round trip: lon/lat → OSGB36 → lon/lat

| | Default (Redfearn) | `karney_tm` |
|:--|:--|:--|
| mean | 0.42 mm | 0.29 mm |
| max  | 4.97 mm (TP31) | 0.61 mm (TP31) |

TP31 is on St Kilda (8.58°W), the test point furthest from the central meridian.
This round trip is the clearest measure of an implementation's internal
consistency, because it does not depend on the millimetre rounding of any
published intermediate.

The default method's maximum round-trip error as a function of longitude offset
from 2°W:

| Point | Longitude | Offset from 2°W | Round-trip error |
|:--|--:|--:|--:|
| TP31 (St Kilda) | 8.58°W | 6.58° | 4.97 mm |
| TP32 | 7.59°W | 5.59° | 1.31 mm |
| TP33 | 6.26°W | 4.26° | 0.66 mm |
| mainland points | — | < 4° | < 0.6 mm |

With `karney_tm`, every point is below 0.61 mm.

### 4.3 A note on the reverse direction

The developer pack's reverse test vectors (`OSGB36_to_ETRS.csv` and the OSGB→ETRS
input file) are constructed so that the OS *truncated* inverse returns the
original coordinates. They therefore embed the truncated method's own
inconsistency: at TP31 the published forward output easting (9587.909 m) and the
published reverse input easting (9587.906 m) differ by ~3 mm. The default
(Redfearn) build reproduces these vectors, by design. The `karney_tm` build does
not reproduce them at the far west, because it does not reproduce the truncation
error they encode — this is expected, not a regression.

A fairer measure of inverse accuracy is to invert the *exact* grid coordinate (the
published forward output) and compare with the original geographic coordinate. At
TP31 this gives 4.67 mm for the default method and 1.50 mm for `karney_tm` (the
latter limited by the millimetre rounding of the published grid coordinate).

### 4.4 Relationship to Grid InQuest II

Grid InQuest II uses the same truncated series as the default build and the same
`f64` (double-precision) arithmetic. A like-for-like comparison of the default
build against the full-precision Grid InQuest routines agrees to floating-point
level (forward: sub-nanometre; inverse: bit-identical). The default build is
therefore not less accurate than Grid InQuest; the far-western round-trip error
described above is shared by both, because it is inherent to the OS-specified
truncated projection. The `karney_tm` build is more accurate than both, at the
cost of diverging from the OS truncated convention by the error it removes.

## 5. Choosing an implementation

- **Use the default** when results must match OS-published values, OS software,
  or Grid InQuest II exactly. Accuracy is sub-millimetre across mainland Great
  Britain and up to ~5 mm in the round trip at the extreme western isles.
- **Use `karney_tm`** when the best absolute projection accuracy is required —
  for example, to keep the round trip sub-millimetre at the western isles, or
  where the reference is the true geodetic projection rather than the OS truncated
  series. Note that results will then differ from the OS truncated convention (and
  Grid InQuest) by a few millimetres at the far west.

## 6. Performance

Measured with `cargo bench --bench projection` over 10 000 representative points
on the development machine (Apple Silicon). Treat the figures as relative; the
ratio between the two methods is the portable quantity.

| Operation | Default | `karney_tm` | Overhead |
|:--|--:|--:|--:|
| Forward, projection only (`convert_etrs89`) | 67 ns/pt | 176 ns/pt | +163% |
| Inverse, projection only (`convert_etrs89_to_ll`) | 186 ns/pt | 310 ns/pt | +67% |
| Forward, full pipeline (`convert_osgb36`) | 141 ns/pt | 251 ns/pt | +77% |
| Inverse, full pipeline (`convert_osgb36_to_ll`) | 409 ns/pt | 510 ns/pt | +25% |

The Karney forward projection is the most affected, because its conformal-latitude
step evaluates several transcendental functions where the truncated series
evaluates a polynomial. End to end, the OSTN15 grid shift and the reverse
iteration dilute the overhead. All conversions remain well under a microsecond
per point.

## 7. Reproducing these figures

- Test suites (both projections):
  `cargo test` and `cargo test --features karney_tm`.
  Relevant tests: `dev_pack_forward_etrs_to_osgb`,
  `test_osgb36_to_etrs89_iterations_detailed`,
  `dev_pack_redfearn_roundtrip_bound`, `dev_pack_karney_roundtrip_submm`.
- Benchmark (A/B against a saved baseline):
  ```
  cargo bench --bench projection -- --save-baseline redfearn
  cargo bench --bench projection --features karney_tm -- --baseline redfearn
  ```
- Reference data: `test_inputs/` (the OSTN15/OSGM15 developer-pack CSVs).

## 8. References

- Ordnance Survey, *A guide to coordinate systems in Great Britain*, Annex C
  (Transverse Mercator projection).
- Ordnance Survey, *OSTN15 transformation and OSGM15 user guide*.
- Ordnance Survey, *OSTN15/OSGM15 developer pack* (test dataset).
- C. F. F. Karney (2011), "Transverse Mercator with an accuracy of a few
  nanometers", *Journal of Geodesy* 85(8): 475–485.
  doi:[10.1007/s00190-011-0445-3](https://doi.org/10.1007/s00190-011-0445-3).
- OSTN15 data: © Crown copyright, Ordnance Survey and the Ministry of Defence
  (MOD) 2016.
