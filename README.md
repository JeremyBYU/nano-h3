# nano-h3

A single-header C++17 converter from lat/lng to an [H3](https://github.com/uber/h3)
cell index at one compile-time resolution. `nanoh3::Grid<11>::cell()` returns
exactly what H3 v4.1.0's `latLngToCell()` returns, about twice as fast.

Use it if you convert millions of points at a resolution you know when you
compile. That is the map-matching workload it was pulled out of. If you need
polygons, directed edges, compaction, or a resolution chosen at runtime, use
H3 instead. This library has three functions.

```cpp
#include "nanoh3/nanoh3.hpp"

auto cell = nanoh3::Grid<11>::cell_deg(45.52, -122.68);
auto same = nanoh3::Grid<11>::cell(0.7945, -2.1412);   // radians

double lat, lng;
nanoh3::Grid<11>::center(cell, lat, lng);             // == cellToLatLng

std::uint64_t ring[7];
nanoh3::Grid<11>::ring1(cell, ring);                  // == gridDisk(cell, 1)
```

No build step and no dependencies. Copy `include/nanoh3/` into your project, or
see [Installing](#installing). The header includes `<cmath>` and `<cstdint>`
and nothing else. H3 itself is needed only to run the tests. Every entry point
is a pure function of its arguments: no state, nothing to own, nothing to say
about threads.

## Speed

AMD Ryzen Threadripper 7960X, gcc 13.3.0, `-O3`, governor `performance`,
CPU-pinned, median of 3 runs of `bench-nanoh3`. One op is one conversion
over 65,536 precomputed points.

| | ns/op | vs H3 |
|---|---:|---:|
| `latLngToCell`, scattered points | 811.5 | 1.00x |
| `Grid<11>::cell`, scattered | 298.5 | **2.72x** |
| `latLngToCell`, GPS trace | 613.9 | 1.00x |
| `Grid<11>::cell`, trace | 260.7 | **2.35x** |
| `Grid<11>::cell_fast`, trace | 134.5 | **4.56x** |
| `Grid<11>::cell_fast`, scattered | 151.5 | **5.36x** |

The gain comes from the fixed resolution: `Grid<Res>` is a template, so the
resolution-dependent branches fold and the scale loop can unroll at compile
time. An integer digit walk replaces H3's `long double` rounding with exact
integer arithmetic, carries two coordinates instead of repeatedly normalizing
three, and encodes directions directly. A trace is faster than scattered
points because the CPU's branch predictor and warm tables reward locality on
their own; the library does nothing for it. Timings vary with compiler, CPU
flags, and machine load.

`cell_fast()` trades exactness for a vector gnomonic projection. A point
within nanometres of a boundary can land in the immediate neighbour of the
true cell, and the suite measures that: zero divergences in 2,000,000 uniform
points at resolution 11, any divergence required to be a neighbour. Use it
when a rare neighbour-cell answer is acceptable and `cell()` when it is not.

## What "bit-identical" means

`cell()`, `center()` and `ring1()` return exactly what H3 v4.1.0's
`latLngToCell()`, `cellToLatLng()` and `gridDisk(k=1)` return, `ring1`'s
output ordering included.

The integer back end is proven: the aperture-7 digit walk replaces H3's
`lroundl((3i - j) / 7.0L)` with exact integer round-half-away division, and
for integer `n` the two always round the same way.

The floating-point front end is tested. It replicates H3's operation sequence
and its mixed double/long-double precision statement by statement, and a
differential suite enforces equality against real H3 over 88 million
assertions: sampled inputs for volume (uniform sphere, random walks, dense
sampling around all 12 pentagons, every base cell), and constructed inputs
that attack specific mechanisms at every resolution: cell boundaries found by
bisection to adjacent doubles, cell vertices and edge midpoints, icosahedron
face ties, poles, the antimeridian, signed zero and subnormals, and exhaustive
enumeration of all 2,352,972 cells at resolutions 0 to 5. The vertex case is
the one that matters: it caught a defect (long-double constants truncated to
double) that 39 million random assertions had not. No proof covers every
representable double for the front end, so that is where the claim stops.

## Limits

- **Input validation: none.** `cell()` has no error channel. For NaN,
  infinity or a latitude outside ±π/2 it returns a structurally valid,
  meaningless index, with no sentinel to check. No crash and no undefined
  behaviour under ASan and UBSan, but validate your own coordinates.
- **`ring1` output layout.** It reproduces `gridDisk` slot for slot. On the
  fast path the origin is at slot 0. Near a pentagon the origin can be at any
  slot and a pentagon's empty slot can be anywhere. Iterate all seven slots,
  skip zeros, compare against your origin. H3 behaves identically.
- **Floating-point flags.** Bit-identity reproduces H3's arithmetic, so flags
  that rewrite it break it. `-O3` passes the suite; `-march=native` needs
  `-ffp-contract=off` (FMA contraction is the only culprit, and turning it
  off costs nothing measurable); `-ffast-math` and `-mlong-double-64` break
  `center()` on ordinary inputs and `cell()` on cell vertices.
- **Platforms.** Linux and macOS x86-64, gcc and clang, Release and Debug in
  CI. Windows is not supported (`M_PI` under MSVC); arm64 macOS is untested.

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

Pin a full commit SHA rather than a tag: FetchContent stops contacting the
remote once `_deps/nano-h3-src` exists, so a moved tag silently gives a fresh
checkout and a working build directory different code. Tests, benchmarks and
install rules are gated on `PROJECT_IS_TOP_LEVEL`, so consuming nano-h3 pulls
in no H3, no doctest and no nanobench; `find_package` works after
`cmake --install`. CMake 3.21 or newer.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/nanoh3-tests          # ~24 s
./build/benchmarks/bench-nanoh3
```

## Status

Extracted from a private project. Used in production for one workload at
resolution 11. Issues and pull requests are welcome; response time is
best-effort. Any change to the front-end arithmetic needs the full suite
green, including the vertex case.

## Licence

Apache-2.0. nano-h3 is a derivative work of H3 v4.1.0: roughly 70% of the
header is a transcription of H3's C sources, and the data tables come from H3
as well. [NOTICE](NOTICE) records which parts derive from what, table by table.

H3 is a project of Uber Technologies, Inc. nano-h3 is an independent derivative
work and is not affiliated with, endorsed by, or supported by Uber
Technologies, Inc. This project is also unrelated to the NanoPi H3
single-board computer.
