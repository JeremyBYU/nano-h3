// nanoh3 differential suite: the library's ONLY correctness authority is
// equality with the real H3 latLngToCell, enforced here over uniform-sphere,
// regional, boundary-adversarial, and cached-vs-uncached inputs, at an odd
// (Class III) and an even (Class II) resolution. Determinism: fixed seeds.
#include <doctest/doctest.h>
#include <h3api.h>

#include <cmath>
#include <cstdint>
#include <random>

#include "nanoh3/nanoh3.hpp"

namespace {

std::uint64_t h3_oracle(double lat, double lng, int res) {
  LatLng g{lat, lng};
  H3Index out = 0;
  REQUIRE(latLngToCell(&g, res, &out) == E_SUCCESS);
  return out;
}

template <int Res>
void check_pt(double lat, double lng, nanoh3::Cache* c = nullptr) {
  const auto want = h3_oracle(lat, lng, Res);
  const auto got = nanoh3::Grid<Res>::cell(lat, lng, c);
  if (got != want) {
    CAPTURE(lat);
    CAPTURE(lng);
    CAPTURE(Res);
  }
  REQUIRE(got == want);
}

}  // namespace

TEST_CASE("nanoh3: face-cache bound is provably safe") {
  // Voronoi ball argument needs bound <= (min pairwise chord / 2)^2.
  double min_sqd = 5.0;
  for (int a = 0; a < 20; ++a) {
    for (int b = a + 1; b < 20; ++b) {
      const auto& u = nanoh3::kFaceCenterPoint[a];
      const auto& v = nanoh3::kFaceCenterPoint[b];
      const double d = (u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) +
                       (u.z - v.z) * (u.z - v.z);
      if (d < min_sqd) min_sqd = d;
    }
  }
  REQUIRE(0.126 < min_sqd / 4.0);  // kSameFaceBound < (min chord / 2)^2
}

TEST_CASE("nanoh3: uniform sphere, res 11 and res 10") {
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int n = 0; n < 2'000'000; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);  // uniform on sphere
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    check_pt<11>(lat, lng);
    check_pt<10>(lat, lng);
  }
}

TEST_CASE("nanoh3: US region with warm cache (the matcher's regime)") {
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> lat_d(24.0, 49.0), lng_d(-125.0, -66.0);
  nanoh3::Cache c11, c10;
  constexpr double kD2R = M_PI / 180.0;
  for (int n = 0; n < 2'000'000; ++n) {
    const double lat = lat_d(rng) * kD2R, lng = lng_d(rng) * kD2R;
    check_pt<11>(lat, lng, &c11);
    check_pt<10>(lat, lng, &c10);
  }
}

TEST_CASE("nanoh3: trace-like random walks (locality + boundary crossings)") {
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> lat0(-65.0, 70.0), lng0(-180.0, 180.0);
  std::normal_distribution<double> step(0.0, 15.0 / 111000.0);  // ~15 m in deg
  constexpr double kD2R = M_PI / 180.0;
  for (int w = 0; w < 200; ++w) {
    double lat = lat0(rng), lng = lng0(rng);
    nanoh3::Cache c;
    for (int n = 0; n < 5000; ++n) {
      lat += step(rng);
      lng += step(rng);
      check_pt<11>(lat * kD2R, lng * kD2R, &c);
    }
  }
}

TEST_CASE("nanoh3: adversarial near-cell-boundary points") {
  // March across cell boundaries in tiny steps: take random points, step
  // toward a random bearing in 1 cm increments for 60 m — every res-11
  // boundary in the path is crossed within 1 cm of an evaluation.
  std::mt19937_64 rng(99);
  std::uniform_real_distribution<double> lat0(-65.0, 70.0), lng0(-180.0, 180.0),
      br(0.0, 2.0 * M_PI);
  constexpr double kD2R = M_PI / 180.0;
  for (int w = 0; w < 300; ++w) {
    double lat = lat0(rng), lng = lng0(rng);
    const double b = br(rng);
    const double dlat = std::sin(b) * 0.01 / 111000.0;
    const double dlng = std::cos(b) * 0.01 / 111000.0;
    nanoh3::Cache c;
    for (int n = 0; n < 6000; ++n) {
      lat += dlat;
      lng += dlng;
      check_pt<11>(lat * kD2R, lng * kD2R, &c);
    }
  }
}

TEST_CASE("nanoh3: pentagon neighborhoods") {
  // All 12 pentagon base cells: sample densely around each pentagon center
  // (the rotation/k-axis code paths live here).
  H3Index pentas[12] = {};
  REQUIRE(getPentagons(11, pentas) == E_SUCCESS);
  std::mt19937_64 rng(5);
  std::normal_distribution<double> jitter(0.0, 0.02);  // ~2 km spread in deg
  constexpr double kD2R = M_PI / 180.0;
  for (const auto p : pentas) {
    LatLng c;
    REQUIRE(cellToLatLng(p, &c) == E_SUCCESS);
    const double clat = c.lat / kD2R, clng = c.lng / kD2R;
    for (int n = 0; n < 20000; ++n) {
      check_pt<11>((clat + jitter(rng)) * kD2R, (clng + jitter(rng)) * kD2R);
    }
  }
}

TEST_CASE("nanoh3: cell_deg replicates degsToRads bit-for-bit") {
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> lat_d(24.0, 49.0), lng_d(-125.0, -66.0);
  nanoh3::Cache c;
  for (int n = 0; n < 500'000; ++n) {
    const double lat = lat_d(rng), lng = lng_d(rng);
    LatLng g{degsToRads(lat), degsToRads(lng)};
    H3Index want = 0;
    REQUIRE(latLngToCell(&g, 11, &want) == E_SUCCESS);
    REQUIRE(nanoh3::Grid<11>::cell_deg(lat, lng, &c) == want);
  }
}

namespace {

template <int Res>
void check_ring(std::uint64_t cell) {
  H3Index want[7] = {};
  const H3Error err = gridDisk(cell, 1, want);
  std::uint64_t got[7] = {};
  nanoh3::Grid<Res>::ring1(cell, got);
  if (err == E_SUCCESS) {
    for (int i = 0; i < 7; ++i) {
      if (got[i] != want[i]) CAPTURE(cell);
      REQUIRE(got[i] == want[i]);
    }
  }
}

}  // namespace

TEST_CASE("nanoh3: ring1 vs gridDisk, random cells") {
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int n = 0; n < 500'000; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    check_ring<11>(h3_oracle(lat, lng, 11));
    check_ring<10>(h3_oracle(lat, lng, 10));
  }
}

TEST_CASE("nanoh3: ring1 in pentagon distortion areas") {
  // Every cell within k=3 of every res-11 pentagon: exercises the deleted-k
  // subsequence, the safe fallback's hash-set order, and rotation
  // adjustments.
  H3Index pentas[12] = {};
  REQUIRE(getPentagons(11, pentas) == E_SUCCESS);
  for (const auto p : pentas) {
    H3Index disk[37] = {};
    int dist[37] = {};
    REQUIRE(gridDiskDistancesSafe(p, 3, disk, dist) == E_SUCCESS);
    for (const auto c : disk) {
      if (c != 0) check_ring<11>(c);
    }
  }
}

TEST_CASE("nanoh3: ring1 across base-cell boundaries") {
  // Res-11 cells on a 0.5 deg global grid hit every base cell and hundreds
  // of base-cell edges.
  constexpr double kD2R = M_PI / 180.0;
  for (double lat = -84.0; lat <= 84.0; lat += 0.5) {
    for (double lng = -180.0; lng < 180.0; lng += 0.5) {
      check_ring<11>(h3_oracle(lat * kD2R, lng * kD2R, 11));
    }
  }
}

TEST_CASE("nanoh3: center vs cellToLatLng (bit-exact)") {
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  auto check_center = [&](std::uint64_t cell, int res) {
    LatLng want;
    REQUIRE(cellToLatLng(cell, &want) == E_SUCCESS);
    double lat = 0, lng = 0;
    if (res == 11) nanoh3::Grid<11>::center(cell, lat, lng);
    else nanoh3::Grid<10>::center(cell, lat, lng);
    if (lat != want.lat || lng != want.lng) CAPTURE(cell);
    REQUIRE(lat == want.lat);
    REQUIRE(lng == want.lng);
  };
  for (int n = 0; n < 500'000; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    check_center(h3_oracle(lat, lng, 11), 11);
    check_center(h3_oracle(lat, lng, 10), 10);
  }
  // Pentagon neighborhoods: leading-5/leading-4 digit paths and overage.
  H3Index pentas[12] = {};
  REQUIRE(getPentagons(11, pentas) == E_SUCCESS);
  for (const auto p : pentas) {
    H3Index disk[37] = {};
    int dist[37] = {};
    REQUIRE(gridDiskDistancesSafe(p, 3, disk, dist) == E_SUCCESS);
    for (const auto c : disk) {
      if (c != 0) check_center(c, 11);
    }
  }
}

TEST_CASE("nanoh3: cell_fast divergence rate (opt-in path)") {
  // cell_fast is NOT bit-identical by design; this measures and bounds its
  // divergence. Correct axes diverge only within float-noise of boundaries
  // (empirically ~1e-7 of uniform points); a convention slip fails at ~100%.
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  long diverged = 0;
  const long N = 2'000'000;
  nanoh3::Cache c;
  for (long n = 0; n < N; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    const auto want = h3_oracle(lat, lng, 11);
    const auto got = nanoh3::Grid<11>::cell_fast(lat, lng, &c);
    if (got != want) {
      ++diverged;
      // Any divergence must still be a NEIGHBOR of the true cell.
      H3Index ring[7] = {};
      REQUIRE(gridDisk(want, 1, ring) == E_SUCCESS);
      bool adjacent = false;
      for (const auto r : ring) adjacent = adjacent || r == got;
      REQUIRE(adjacent);
    }
  }
  MESSAGE("cell_fast divergence: ", diverged, " / ", N);
  REQUIRE(diverged < N / 100'000);  // < 1e-5 rate, orders below noise floors
}

TEST_CASE("nanoh3: poles and antimeridian") {
  constexpr double kD2R = M_PI / 180.0;
  for (double lat : {89.9999, 89.0, -89.0, -89.9999, 0.0}) {
    for (double lng : {179.9999, -179.9999, 180.0, -180.0, 0.0, 90.0}) {
      check_pt<11>(lat * kD2R, lng * kD2R);
      check_pt<10>(lat * kD2R, lng * kD2R);
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Breadth across every resolution Grid<Res> advertises.
//
// Everything above is DEPTH at res 11 (Class III) and res 10 (Class II), which
// is what the map matcher this came from actually used. The library claims
// 0..15, so the other fourteen need evidence too. This is deliberately
// shallower per resolution: the expensive adversarial cases stay at 10/11, and
// this adds ~1 s rather than ~2 min.
// ─────────────────────────────────────────────────────────────────────────────
namespace {

template <int Res>
void sweep_one_res(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  nanoh3::Cache cache;

  for (int n = 0; n < 12'000; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    const auto want = h3_oracle(lat, lng, Res);

    const auto got_cold = nanoh3::Grid<Res>::cell(lat, lng);
    const auto got_warm = nanoh3::Grid<Res>::cell(lat, lng, &cache);
    if (got_cold != want || got_warm != want) {
      CAPTURE(lat);
      CAPTURE(lng);
      CAPTURE(Res);
    }
    REQUIRE(got_cold == want);
    REQUIRE(got_warm == want);  // the face cache can never change the answer

    LatLng want_c;
    REQUIRE(cellToLatLng(want, &want_c) == E_SUCCESS);
    double clat = 0, clng = 0;
    nanoh3::Grid<Res>::center(want, clat, clng);
    if (clat != want_c.lat || clng != want_c.lng) CAPTURE(Res);
    REQUIRE(clat == want_c.lat);
    REQUIRE(clng == want_c.lng);

    std::uint64_t got_ring[7] = {};
    H3Index want_ring[7] = {};
    nanoh3::Grid<Res>::ring1(want, got_ring);
    if (gridDisk(want, 1, want_ring) == E_SUCCESS) {
      for (int i = 0; i < 7; ++i) REQUIRE(got_ring[i] == (std::uint64_t)want_ring[i]);
    }
  }

  // Pentagons AT THIS RESOLUTION, plus their k=2 disks: the rotation and
  // deleted-k paths differ per resolution class.
  H3Index pentas[12] = {};
  REQUIRE(getPentagons(Res, pentas) == E_SUCCESS);
  for (const auto p : pentas) {
    H3Index disk[19] = {};
    int dist[19] = {};
    REQUIRE(gridDiskDistancesSafe(p, 2, disk, dist) == E_SUCCESS);
    for (const auto c : disk) {
      if (c == 0) continue;
      LatLng cc;
      REQUIRE(cellToLatLng(c, &cc) == E_SUCCESS);
      REQUIRE(nanoh3::Grid<Res>::cell(cc.lat, cc.lng) == c);
      double clat = 0, clng = 0;
      nanoh3::Grid<Res>::center(c, clat, clng);
      REQUIRE(clat == cc.lat);
      REQUIRE(clng == cc.lng);
      std::uint64_t got_ring[7] = {};
      H3Index want_ring[7] = {};
      nanoh3::Grid<Res>::ring1(c, got_ring);
      if (gridDisk(c, 1, want_ring) == E_SUCCESS) {
        for (int i = 0; i < 7; ++i) REQUIRE(got_ring[i] == (std::uint64_t)want_ring[i]);
      }
    }
  }

  // Poles and antimeridian, where lng wrapping and the polar pentagons live.
  constexpr double kD2R = M_PI / 180.0;
  for (double lat : {89.9999, 89.0, 0.0, -89.0, -89.9999}) {
    for (double lng : {179.9999, -179.9999, 180.0, -180.0, 0.0}) {
      REQUIRE(nanoh3::Grid<Res>::cell(lat * kD2R, lng * kD2R) ==
              h3_oracle(lat * kD2R, lng * kD2R, Res));
      REQUIRE(nanoh3::Grid<Res>::cell_deg(lat, lng) ==
              h3_oracle(degsToRads(lat), degsToRads(lng), Res));
    }
  }
}

template <std::size_t... Rs>
void sweep_all(std::mt19937_64& rng, std::index_sequence<Rs...>) {
  (sweep_one_res<static_cast<int>(Rs)>(rng), ...);
}

}  // namespace

TEST_CASE("nanoh3: every resolution 0..15, cell/center/ring1 vs H3") {
  std::mt19937_64 rng(20260816);
  sweep_all(rng, std::make_index_sequence<16>{});
}

// ─────────────────────────────────────────────────────────────────────────────
// Cell vertices: the input class that every other case in this file misses.
//
// This exists because it FOUND A BUG. cell() truncated H3's long-double
// constants to double, and a comment claimed the suite had measured that
// bit-identical. It had not: the divergence band around a vertex is only a few
// nanometres wide, and random sampling plus a 1 cm boundary march never lands
// in it. Feeding cellToBoundary vertices straight back into cell() showed
// 42,253 of 192,727 probes wrong (21.9%), at every resolution except 0.
//
// Vertices are the natural adversarial input because they are equidistant from
// three cells, so they maximise the chance that a last-ulp difference changes
// the answer. They are also a realistic input: anyone tiling or filling
// polygons feeds boundary coordinates back in.
// ─────────────────────────────────────────────────────────────────────────────
namespace {

template <int Res>
void vertex_probe_one_res(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int n = 0; n < 1'500; ++n) {
    const double lat = std::asin(2.0 * u(rng) - 1.0);
    const double lng = (2.0 * u(rng) - 1.0) * M_PI;
    CellBoundary cb;
    REQUIRE(cellToBoundary(h3_oracle(lat, lng, Res), &cb) == E_SUCCESS);
    for (int v = 0; v < cb.numVerts; ++v) {
      const auto want = h3_oracle(cb.verts[v].lat, cb.verts[v].lng, Res);
      const auto got = nanoh3::Grid<Res>::cell(cb.verts[v].lat, cb.verts[v].lng);
      if (got != want) {
        CAPTURE(Res);
        CAPTURE(cb.verts[v].lat);
        CAPTURE(cb.verts[v].lng);
      }
      REQUIRE(got == want);
    }
  }
}

template <std::size_t... Rs>
void vertex_probe_all(std::mt19937_64& rng, std::index_sequence<Rs...>) {
  (vertex_probe_one_res<static_cast<int>(Rs)>(rng), ...);
}

}  // namespace

TEST_CASE("nanoh3: cell() on exact cell vertices, every resolution") {
  std::mt19937_64 rng(4242);
  vertex_probe_all(rng, std::make_index_sequence<16>{});
}
