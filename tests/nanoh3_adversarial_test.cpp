// Constructed adversarial inputs, at every resolution.
//
// nanoh3_test.cpp samples: uniform points, random walks, a 1 cm boundary march.
// That is how a real bit-identity defect survived 39 million assertions. The
// divergence band around a cell vertex is a few nanometres wide, and no amount
// of random sampling reliably lands in a band that thin. More points do not
// fix this; different points do.
//
// So the cases here do not sample. Each one CONSTRUCTS the input that is
// hardest for a specific mechanism in the library:
//
//   bisection      the exact decision boundary between two cells, anywhere,
//                  found by halving until two doubles are adjacent
//   edge midpoint  equidistant from exactly two cells
//   face tie       equidistant from two icosahedron face centres, where the
//                  20-face argmin has no unique winner
//   degenerate     poles, antimeridian, signed zero, subnormals, and their
//                  one-ulp neighbours
//
// Every entry point is a pure function of its arguments, so there is no state
// to corrupt and nothing to check about call ORDER. That was not true while the
// face cache existed, and two test cases here existed only to police it.
//
// The vertex case (equidistant from three cells) lives in nanoh3_test.cpp
// because it is the one that found the bug and belongs next to the contract.
#include <doctest/doctest.h>
#include <h3api.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

#include "nanoh3/nanoh3.hpp"

namespace {

template <int Res>
void require_matches(double lat, double lng) {
  LatLng g{lat, lng};
  H3Index want = 0;
  if (latLngToCell(&g, Res, &want) != E_SUCCESS) return;  // out of domain
  const auto got = nanoh3::Grid<Res>::cell(lat, lng);
  if (got != static_cast<std::uint64_t>(want)) {
    CAPTURE(Res);
    CAPTURE(lat);
    CAPTURE(lng);
  }
  REQUIRE(got == static_cast<std::uint64_t>(want));
}

// Halve toward the true boundary until the two endpoints are adjacent doubles,
// then check both sides and the surrounding ulp grid. This is the sharpest
// probe available and, unlike the vertex case, it works ANYWHERE rather than
// only at the handful of geometrically special points.
template <int Res>
void bisect_to_boundary(std::mt19937_64& rng, int trials) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int t = 0; t < trials; ++t) {
    double alat = std::asin(2.0 * u(rng) - 1.0), alng = (2.0 * u(rng) - 1.0) * M_PI;
    double blat = alat + 1e-4, blng = alng + 1e-4;
    LatLng ga{alat, alng}, gb{blat, blng};
    H3Index ca = 0, cb = 0;
    if (latLngToCell(&ga, Res, &ca) != E_SUCCESS) continue;
    if (latLngToCell(&gb, Res, &cb) != E_SUCCESS) continue;
    if (ca == cb) continue;  // no boundary in this interval at this resolution

    for (int i = 0; i < 200; ++i) {
      const double mlat = 0.5 * (alat + blat), mlng = 0.5 * (alng + blng);
      if ((mlat == alat && mlng == alng) || (mlat == blat && mlng == blng)) break;
      LatLng gm{mlat, mlng};
      H3Index cm = 0;
      if (latLngToCell(&gm, Res, &cm) != E_SUCCESS) break;
      if (cm == ca) { alat = mlat; alng = mlng; } else { blat = mlat; blng = mlng; }
    }

    const double inf_d = std::numeric_limits<double>::infinity();
    require_matches<Res>(alat, alng);
    require_matches<Res>(blat, blng);
    for (double dl : {std::nextafter(alat, -inf_d), alat, std::nextafter(alat, inf_d)}) {
      for (double dg : {std::nextafter(alng, -inf_d), alng, std::nextafter(alng, inf_d)}) {
        require_matches<Res>(dl, dg);
      }
    }
  }
}

template <int Res>
void edge_midpoints(std::mt19937_64& rng, int cells) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int n = 0; n < cells; ++n) {
    LatLng g{std::asin(2.0 * u(rng) - 1.0), (2.0 * u(rng) - 1.0) * M_PI};
    H3Index c = 0;
    if (latLngToCell(&g, Res, &c) != E_SUCCESS) continue;
    CellBoundary cb;
    if (cellToBoundary(c, &cb) != E_SUCCESS) continue;
    for (int v = 0; v < cb.numVerts; ++v) {
      const LatLng& a = cb.verts[v];
      const LatLng& b = cb.verts[(v + 1) % cb.numVerts];
      const double mlat = 0.5 * (a.lat + b.lat), mlng = 0.5 * (a.lng + b.lng);
      require_matches<Res>(mlat, mlng);
      const double inf_d = std::numeric_limits<double>::infinity();
      require_matches<Res>(std::nextafter(mlat, -inf_d), mlng);
      require_matches<Res>(std::nextafter(mlat, inf_d), mlng);
    }
  }
}

// Points exactly between two adjacent face centres, where the 20-face argmin
// has no unique winner and the result depends on iteration order and on which
// comparison the strict < resolves first. H3 makes the same choice; this pins
// that it keeps doing so.
template <int Res>
void face_ties() {
  for (int a = 0; a < 20; ++a) {
    for (int b = a + 1; b < 20; ++b) {
      const auto& u = nanoh3::kFaceCenterPoint[a];
      const auto& v = nanoh3::kFaceCenterPoint[b];
      const double d = (u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) +
                       (u.z - v.z) * (u.z - v.z);
      if (d > 0.52) continue;  // adjacent faces only
      double x = u.x + v.x, y = u.y + v.y, z = u.z + v.z;
      const double n = std::sqrt(x * x + y * y + z * z);
      x /= n; y /= n; z /= n;
      const double lat = std::asin(z), lng = std::atan2(y, x);

      require_matches<Res>(lat, lng);
    }
  }
}

template <int Res>
void degenerate_inputs() {
  const double tiny = std::numeric_limits<double>::denorm_min();
  const double smallest = std::numeric_limits<double>::min();
  for (double lat : {M_PI_2, -M_PI_2, std::nextafter(M_PI_2, 0.0),
                     std::nextafter(-M_PI_2, 0.0), 0.0, -0.0, smallest, -smallest,
                     tiny, std::nextafter(0.0, 1.0), std::nextafter(0.0, -1.0)}) {
    for (double lng : {M_PI, -M_PI, std::nextafter(M_PI, 0.0),
                       std::nextafter(-M_PI, 0.0), 0.0, -0.0, M_PI_2, -M_PI_2, tiny}) {
      require_matches<Res>(lat, lng);
    }
  }
}

template <std::size_t... Rs>
void all_res_bisect(std::mt19937_64& rng, std::index_sequence<Rs...>) {
  (bisect_to_boundary<static_cast<int>(Rs)>(rng, 250), ...);
}
template <std::size_t... Rs>
void all_res_edges(std::mt19937_64& rng, std::index_sequence<Rs...>) {
  (edge_midpoints<static_cast<int>(Rs)>(rng, 200), ...);
}
template <std::size_t... Rs>
void all_res_ties(std::index_sequence<Rs...>) {
  (face_ties<static_cast<int>(Rs)>(), ...);
}
template <std::size_t... Rs>
void all_res_degenerate(std::index_sequence<Rs...>) {
  (degenerate_inputs<static_cast<int>(Rs)>(), ...);
}

}  // namespace

TEST_CASE("nanoh3: bisection to the exact cell boundary, every resolution") {
  std::mt19937_64 rng(20260817);
  all_res_bisect(rng, std::make_index_sequence<16>{});
}

TEST_CASE("nanoh3: cell edge midpoints, every resolution") {
  std::mt19937_64 rng(20260817);
  all_res_edges(rng, std::make_index_sequence<16>{});
}

TEST_CASE("nanoh3: icosahedron face ties, every resolution") {
  all_res_ties(std::make_index_sequence<16>{});
}

TEST_CASE("nanoh3: poles, antimeridian, signed zero, subnormals") {
  all_res_degenerate(std::make_index_sequence<16>{});
}

// Documents behaviour rather than asserting a contract, because there is no
// contract to assert: cell() has no error channel. H3 rejects non-finite and
// out-of-range input with E_LATLNG_DOMAIN; nanoh3 cannot, so it returns
// SOMETHING. Measured over NaN, +/-inf, 1e308 and 300k random bit patterns:
// no crash and no undefined behaviour under UBSan including
// float-cast-overflow, but the returned index is meaningless, non-zero, and
// indistinguishable from a real cell. Callers must validate their own input.
TEST_CASE("nanoh3: out-of-domain input does not crash, and is not detectable") {
  const double nan = std::nan("");
  const double inf = std::numeric_limits<double>::infinity();
  long garbage = 0;
  for (double lat : {nan, inf, -inf, 1e308, -1e308, 6.0, -6.0}) {
    for (double lng : {nan, inf, -inf, 1e308, 4.0 * M_PI, 0.0}) {
      LatLng g{lat, lng};
      H3Index want = 0;
      const H3Error err = latLngToCell(&g, 11, &want);
      const std::uint64_t got = nanoh3::Grid<11>::cell(lat, lng);
      if (err != E_SUCCESS) {
        ++garbage;
        // The point of the test: no sentinel is returned. If this ever starts
        // failing because a validity check was added, that is an improvement,
        // and the README paragraph about it must change with it.
        REQUIRE(isValidCell(got) == 1);
      } else {
        REQUIRE(got == static_cast<std::uint64_t>(want));
      }
    }
  }
  MESSAGE("out-of-domain inputs returning a valid-looking but meaningless cell: ",
          garbage);
  REQUIRE(garbage > 0);
}
