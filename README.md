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

nanoh3::Cache cache;                                  // optional, holds the last face
auto cell = nanoh3::Grid<11>::cell_deg(45.52, -122.68, &cache);
auto same = nanoh3::Grid<11>::cell(0.7945, -2.1412, &cache);   // radians

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
| `latLngToCell`, scattered points | 804.2 | 1.00x |
| `Grid<11>::cell`, scattered, cache | 386.1 | **2.08x** |
| `latLngToCell`, GPS trace | 551.3 | 1.00x |
| `Grid<11>::cell`, trace, cache | 244.8 | **2.25x** |
| `Grid<11>::cell_fast`, trace, cache | 134.2 | **4.11x** |

Two mechanisms produce those two rows.

Scattered points win on **fixed resolution** alone. `Grid<Res>` is a template,
so the resolution-dependent branches and the scale loop fold at compile time.
The face cache almost never helps here, and the integer digit walk is the whole
gain.

The trace row adds **spatial locality**. Consecutive points in a track share an
icosahedron face essentially always, so `Cache` holds the last face and a point
inside that face's inscribed spherical cap skips the 20-face search. The bound
is half the minimum chord distance between face centres, squared, so a Voronoi
argument guarantees the cached face really is the nearest one. The fallback is
the full search, so the fast path cannot change the answer, and a test
recomputes the true minimum from the face table to confirm the bound.

Run-to-run spread on this machine is about 10% even for identical code. Treat
smaller differences as noise.

## Scope and limits

Read this before adopting. It converts a point to a cell at one resolution,
returns a cell's centre, and returns a cell's six neighbours. Everything else
H3 does is out of scope.

**Resolutions.** All 16, from 0 to 15, are checked against H3 for `cell`,
`center` and `ring1`. Resolutions 10 and 11 additionally carry the deep
adversarial cases, because those are what the original workload used.

**Platforms.** Linux and macOS x86-64, gcc and clang, tested in CI on Release
and Debug. Windows is not supported: `M_PI` needs `_USE_MATH_DEFINES` under
MSVC in three places and nothing has been tested there. arm64 macOS is
untested.

**Threads.** `Cache` is one `int` and carries no synchronization. Give each
thread its own, or pass `nullptr`, which is always correct and always does the
full face search.

**Pentagons.** Twelve cells per resolution have five neighbours rather than
six, so `ring1` fills six of its seven slots and zeroes the rest. The empty
slot is not at a fixed index. This matches H3 exactly, hole position included,
which is the point: `ring1` reproduces `gridDisk`'s output slot for slot rather
than tidying it up.

**Floating-point flags.** Bit-identity means reproducing H3's arithmetic
exactly, so flags that licence the compiler to rewrite that arithmetic break
it. Measured on the full suite:

| Flags | Result |
|---|---|
| `-O3` | all 15 cases pass |
| `-O3 -march=native` | 3 cases fail |
| `-O3 -march=native -ffp-contract=off` | all 15 cases pass |
| `-O3 -ffast-math` | 3 cases fail |
| `-O3 -mlong-double-64` | 3 cases fail |

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
differential suite enforces equality against real H3 over 41,893,661
assertions: uniform-sphere points, a regional distribution, 200 random walks at
15 m steps, an adversarial 1 cm march across cell boundaries, dense sampling
around all 12 pentagons, every base cell on a 0.5 degree global grid, the poles
and antimeridian, and every vertex of 24,000 cell boundaries across all 16
resolutions. No proof covers every representable double, so that is where the
claim stops.

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
static_assert(NANOH3_VERSION_MAJOR == 0 && NANOH3_VERSION_MINOR == 1 &&
                  NANOH3_VERSION_PATCH == 0,
              "nanoh3 version mismatch: the compiled header is not the pinned one");
```

Building the tests:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/nanoh3-tests          # ~20 s
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
