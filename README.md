# nano-h3

A single-header C++17 converter from lat/lng to an [H3](https://github.com/uber/h3)
cell index at one compile-time resolution. `nanoh3::Grid<11>::cell()` returns
exactly what H3 v4.1.0's `latLngToCell()` returns, 2.1x faster on scattered
points and 2.3x faster on a GPS trace.

Use it if you convert millions of points at a resolution you know when you
compile, and especially if consecutive points are near each other. That is the
map-matching workload it was pulled out of. If you need polygons, directed
edges, compaction, or a resolution chosen at runtime, use H3 instead. This
library has three functions.

```cpp
#include "nanoh3/nanoh3.hpp"

auto cell = nanoh3::Grid<11>::cell_deg(45.52, -122.68);
auto same = nanoh3::Grid<11>::cell(0.7945, -2.1412);           // radians

double lat, lng;
nanoh3::Grid<11>::center(cell, lat, lng);             // == cellToLatLng

std::uint64_t ring[7];
nanoh3::Grid<11>::ring1(cell, ring);                  // == gridDisk(cell, 1)
```

No build step and no dependencies. Copy `include/nanoh3/` into your project, or
see [Installing](#installing). The header includes `<cmath>` and `<cstdint>`
and nothing else. H3 itself is needed only to run the tests.

## Speed, and where it comes from

AMD Ryzen Threadripper 7960X, gcc 13.3.0, `-O3`, governor `performance`, median
of 3 runs. One op is one conversion over 65,536 precomputed points.

| | ns/op | vs H3 |
|---|---:|---:|
| `latLngToCell`, scattered points | 808.6 | 1.00x |
| `Grid<11>::cell`, scattered globally | 387.5 | **2.09x** |
| `Grid<11>::cell`, scattered on one face | 323.1 | |
| `latLngToCell`, GPS trace | 619.0 | 1.00x |
| `Grid<11>::cell`, trace | 273.9 | **2.26x** |
| `Grid<11>::cell_fast`, trace | 150.1 | **4.12x** |

**Fixed resolution** is the mechanism that always pays. `Grid<Res>` is a
template, so the resolution-dependent branches and the scale loop fold at
compile time, and the integer digit walk replaces H3's `long double` rounding
with exact integer arithmetic. That is the whole of the scattered-points gain.

**Spatial locality** is worth another 113 ns on top, and the reason is worth
stating precisely because the obvious explanation is wrong. Confining scattered
points to a single icosahedron face saves 60 ns; making them sequential saves 53
ns more. Branch mispredictions fall from 6.86% to 1.19%. The gain is real and it
is delivered by the CPU's branch predictor and warm lookup tables, for free.

There used to be an explicit single-face cache here, on the theory that a track
stays on one face so the 20-face search could be skipped. The theory was correct
and the optimisation was worthless. The cache hit 99.9999% of the time and was
still worth between −1.1% and +2.8% depending on function and optimisation level,
because twenty unrollable distance computations vectorise into almost nothing.
The hardware was already exploiting the locality the cache was built to exploit.

It was deleted in 0.2.0. Removing it removed the library's only mutable state,
so every entry point is now a pure function of its arguments, and the entire
question of thread safety disappeared with it.

Run-to-run spread on this machine is about 10% even for identical code. Treat
smaller differences as noise.

## Scope and limits

Read this before adopting. It converts a point to a cell at one resolution,
returns a cell's centre, and returns a cell's six neighbours. Everything else
H3 does is out of scope.

**Resolutions.** All 16, from 0 to 15. For resolutions 0 through 5 the
coverage is exhaustive rather than sampled: every one of the 2,352,972 cells
that exist at those resolutions is checked for `center`, `ring1` and the
centre-to-cell round trip. Above res 5 the space grows sevenfold per level, so
coverage there is sampled and adversarial. Resolutions 10 and 11 additionally
carry the deepest cases, because those are what the original workload used.

**Platforms.** Linux and macOS x86-64, gcc and clang, tested in CI on Release
and Debug. Windows is not supported: `M_PI` needs `_USE_MATH_DEFINES` under
MSVC in three places and nothing has been tested there. arm64 macOS is
untested.

**Threads.** Nothing to say, which is the point. Every entry point is a pure
function of its arguments. There is no shared mutable state, no per-thread setup,
nothing to own. The one internal static is a read-only per-face table
initialised on first use, and 16 threads released from a spin barrier to race
that initialisation run clean under ThreadSanitizer.

**Input validation: there is none.** H3 rejects non-finite or out-of-range
coordinates with `E_LATLNG_DOMAIN`. `cell()` returns a `uint64_t` and has no
error channel, so for NaN, infinity, or a latitude outside +/-pi/2 it returns a
cell index that is **structurally valid and completely meaningless**. There is
no sentinel to check. Measured over NaN, +/-inf, 1e308 and 300,000 random bit
patterns: no crash and no undefined behaviour, under ASan and under UBSan with
`float-cast-overflow` enabled, but every out-of-domain input produced a
plausible-looking wrong answer. Validate your own coordinates.

**`ring1` output layout.** Do not assume a layout. `ring1` reproduces
`gridDisk`'s output slot for slot rather than tidying it up, and `gridDisk` has
two paths. On the fast path you get seven cells with the origin at slot 0. Near
a pentagon it falls back to a hash set keyed by `origin % 7`, and then:

- the origin can be at **any** slot, not slot 0;
- a pentagon has five neighbours, so one slot is 0, and that empty slot can also
  be at any index.

Measured at resolution 2: 61 cells place the origin somewhere other than slot 0,
and only 10 of those are pentagons themselves. The other 51 are ordinary
hexagons in a pentagon's distortion area. So code that reads `out[0]` as the
origin, or skips index 0 to iterate "just the neighbours", is wrong for those
cells. Iterate all seven slots, skip zeros, and compare against your origin
explicitly. H3 behaves identically; this is faithfulness, not a quirk.

**Floating-point flags.** Bit-identity means reproducing H3's arithmetic
exactly, so flags that licence the compiler to rewrite that arithmetic break
it. Measured on the full suite:

| Flags | Result |
|---|---|
| `-O3` | all 24 cases pass |
| `-O3 -march=native` | 7 cases fail |
| `-O3 -march=native -ffp-contract=off` | all 24 cases pass |
| `-O3 -ffast-math` | 7 cases fail |
| `-O3 -mlong-double-64` | 7 cases fail |

If you use `-march=native`, add `-ffp-contract=off`. The only culprit there is
FMA contraction, and turning it off costs nothing measurable.

Do not use `-ffast-math` or `-mlong-double-64` if you need the contract. They
break `center()` on ordinary inputs and `cell()` on cell vertices. `cell()`
still agrees with H3 on randomly sampled points under those flags, which is
precisely why that is not a safe thing to rely on.

## What "bit-identical" means here

`cell()`, `center()` and `ring1()` return exactly what H3 v4.1.0's
`latLngToCell()`, `cellToLatLng()` and `gridDisk(k=1)` return, `ring1`'s output
ordering included. The two halves of that rest on different evidence.

The integer back end is proven. The aperture-7 digit walk replaces H3's
`lroundl((3i - j) / 7.0L)` with exact integer round-half-away division. For
integer `n`, `n/7` is never a half-integer, its distance from one is at least
1/14, and `long double` rounds `n/7` to within 2^-63. Both round to the same
integer, for every reachable input.

The floating-point front end is tested. It replicates H3's operation sequence
and its mixed double/long-double precision statement by statement, and a
differential suite enforces equality against real H3 over 88,926,394 assertions.
Two kinds of input go into that number, and the second kind matters more than
the first.

Sampled, to cover volume: uniform-sphere points, a regional distribution, 200
random walks at 15 m steps, dense sampling around all 12 pentagons, and every
base cell on a 0.5 degree global grid.

Constructed, to attack specific mechanisms, at every resolution:

- **cell boundaries found by bisection.** Halve toward the boundary between two
  cells until the endpoints are adjacent doubles, then probe both sides and the
  surrounding ulp grid. Unlike a fixed-step march this finds the exact decision
  boundary anywhere, not only at geometrically special points.
- **cell vertices**, equidistant from three cells: every vertex of 24,000 cell
  boundaries.
- **cell edge midpoints**, equidistant from two.
- **icosahedron face ties**, equidistant from two face centres, where the
  20-face search has no unique winner and the answer depends on iteration order.
- **poles, the antimeridian, signed zero and subnormals**, with their one-ulp
  neighbours.
- **exhaustive enumeration** of all 2,352,972 cells at resolutions 0 to 5, which
  is a complete proof for `center` and `ring1` at those resolutions rather than
  evidence about them.

No proof covers every representable double for the coordinate front end, so that
is where the claim stops.

That last input class earns its place, and it is the one most libraries skip.
Cell vertices are equidistant from three cells, which makes them the sharpest
probe available, and they are exactly what you get if you feed `cellToBoundary`
output back in. During development they caught a real bit-identity defect that
39 million randomly sampled assertions had not: the front end was truncating
H3's long-double constants to double, which is invisible everywhere except
within a few nanometres of a cell corner, where it was wrong 22% of the time.

`cell_fast()` is the deliberate exception. It trades exactness for a vector
gnomonic projection: three dot products and a divide against precomputed
per-face axes, at 4.11x instead of 2.25x. A last-ulp difference can flip a
point sitting within nanometres of a boundary, and the answer is then an
immediate neighbour of the true cell. The suite measures that rather than
assuming it: zero divergences in 2,000,000 uniform points at resolution 11,
with any divergence required to be an immediate neighbour. Reach for it when a
rare neighbour-cell answer is acceptable, and use `cell()` when it is not.

## Installing

```cmake
include(FetchContent)
FetchContent_Declare(nano-h3
  GIT_REPOSITORY https://github.com/JeremyBYU/nano-h3.git
  GIT_TAG        <a full 40-character commit SHA>
  GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(nano-h3)

target_link_libraries(your-target PRIVATE nano-h3::nano-h3)
```

Pin a full commit SHA rather than a tag. FetchContent stops contacting the
remote once `_deps/nano-h3-src` exists, so a moved tag silently gives a fresh
checkout and a working build directory different code.

Tests, benchmarks and install rules are gated on `PROJECT_IS_TOP_LEVEL`, so
consuming nano-h3 pulls in no H3, no doctest and no nanobench. `find_package`
works too, after `cmake --install`. CMake 3.21 or newer is required, for
`PROJECT_IS_TOP_LEVEL` itself.

To make a mismatched header a compile error rather than a wrong answer:

```cpp
static_assert(NANOH3_VERSION_MAJOR == 0 && NANOH3_VERSION_MINOR == 2 &&
                  NANOH3_VERSION_PATCH == 0,
              "nanoh3 version mismatch: the compiled header is not the pinned one");
```

Building the tests:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/nanoh3-tests          # ~24 s
./build/benchmarks/bench-nanoh3
```

## Status

Extracted from a private project and published because the differential suite
seemed worth sharing. It is used in production for one workload at resolution
11. Issues and pull requests are welcome, but treat response time as
best-effort. Any change that touches the front-end arithmetic needs the full
suite green, including the vertex case.

## Licence

Apache-2.0. nano-h3 is a derivative work of H3 v4.1.0: roughly 70% of the
header is a transcription of H3's C sources, and the data tables come from H3
as well. [NOTICE](NOTICE) records which parts derive from what, table by table.

H3 is a project of Uber Technologies, Inc. nano-h3 is an independent derivative
work and is not affiliated with, endorsed by, or supported by Uber
Technologies, Inc. This project is also unrelated to the NanoPi H3
single-board computer.
