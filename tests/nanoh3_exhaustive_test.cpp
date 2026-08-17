// Exhaustive coverage, res 0 to 5. Not sampling: EVERY cell.
//
// Sampling answers "is it probably right". Enumeration answers "is it right",
// and for the low resolutions the cell space is small enough to just ask:
//
//   res 0        122 cells        res 3     41,162
//   res 1        842              res 4    288,122
//   res 2      5,882              res 5  2,016,842
//
// That is 2.35 million cells and about 3 s. For center() and ring1(), whose
// input is a cell rather than a coordinate, this is a complete proof over the
// whole domain at those resolutions, not evidence about it. Above res 5 the
// space grows 7x per level and enumeration stops being possible, which is what
// the sampled and adversarial cases are for.
//
// Three things are checked per cell, plus structural invariants that hold
// independently of H3, so a shared misunderstanding of H3 cannot make both
// sides agree on something wrong.
#include <doctest/doctest.h>
#include <h3api.h>

#include <cstdint>
#include <set>
#include <vector>

#include "nanoh3/nanoh3.hpp"

namespace {

template <int Res>
long exhaustive_one_res() {
  std::vector<H3Index> res0(res0CellCount());
  REQUIRE(getRes0Cells(res0.data()) == E_SUCCESS);

  long cells = 0;
  for (const H3Index base : res0) {
    std::int64_t n = 0;
    REQUIRE(cellToChildrenSize(base, Res, &n) == E_SUCCESS);
    std::vector<H3Index> kids(static_cast<std::size_t>(n));
    REQUIRE(cellToChildren(base, Res, kids.data()) == E_SUCCESS);

    for (const H3Index c : kids) {
      if (c == 0) continue;  // pentagon children are sparse
      ++cells;

      // 1. center() is bit-exact against cellToLatLng.
      LatLng want_c;
      REQUIRE(cellToLatLng(c, &want_c) == E_SUCCESS);
      double lat = 0, lng = 0;
      nanoh3::Grid<Res>::center(c, lat, lng);
      if (lat != want_c.lat || lng != want_c.lng) CAPTURE(Res);
      REQUIRE(lat == want_c.lat);
      REQUIRE(lng == want_c.lng);

      // 2. Round trip. A cell's own centre must land back in that cell. This
      //    is the one invariant that ties cell() and center() together, and it
      //    needs no oracle at all.
      REQUIRE(nanoh3::Grid<Res>::cell(want_c.lat, want_c.lng) ==
              static_cast<std::uint64_t>(c));

      // 3. ring1() slot for slot against gridDisk.
      H3Index want_r[7] = {};
      std::uint64_t got_r[7] = {};
      nanoh3::Grid<Res>::ring1(c, got_r);
      if (gridDisk(c, 1, want_r) != E_SUCCESS) continue;
      for (int i = 0; i < 7; ++i) {
        if (got_r[i] != static_cast<std::uint64_t>(want_r[i])) {
          CAPTURE(Res);
          CAPTURE(i);
        }
        REQUIRE(got_r[i] == static_cast<std::uint64_t>(want_r[i]));
      }

      // Structural invariants, independent of H3.
      //
      // NOTE the origin is NOT required at slot 0. On the fast spiral path it
      // is, but pentagons and their distortion areas take the safe hash-set
      // fallback, which places the origin at origin % 7. Measured at res 2: 61
      // cells put it elsewhere, and only 10 of those are pentagons themselves.
      // H3 does exactly the same, so requiring slot 0 here would be asserting
      // something false about H3 rather than something true about nanoh3.
      int nonzero = 0, origin_count = 0;
      std::set<std::uint64_t> distinct;
      for (int i = 0; i < 7; ++i) {
        if (got_r[i] == 0) continue;
        ++nonzero;
        distinct.insert(got_r[i]);
        if (got_r[i] == static_cast<std::uint64_t>(c)) ++origin_count;
      }
      REQUIRE(origin_count == 1);                        // present exactly once
      REQUIRE(static_cast<int>(distinct.size()) == nonzero);  // no duplicates
      REQUIRE((nonzero == 7 || nonzero == 6));           // 6 only at a pentagon

      // Adjacency is symmetric: if b is in our ring, we are in b's ring.
      for (int i = 0; i < 7; ++i) {
        if (got_r[i] == 0 || got_r[i] == static_cast<std::uint64_t>(c)) continue;
        std::uint64_t back[7] = {};
        nanoh3::Grid<Res>::ring1(got_r[i], back);
        bool found = false;
        for (int j = 0; j < 7; ++j) found = found || back[j] == static_cast<std::uint64_t>(c);
        if (!found) {
          CAPTURE(Res);
          CAPTURE(got_r[i]);
        }
        REQUIRE(found);
      }
    }
  }
  return cells;
}

}  // namespace

TEST_CASE("nanoh3: EVERY cell at res 0-5, center/roundtrip/ring1 and invariants") {
  long total = 0;
  total += exhaustive_one_res<0>();
  total += exhaustive_one_res<1>();
  total += exhaustive_one_res<2>();
  total += exhaustive_one_res<3>();
  total += exhaustive_one_res<4>();
  total += exhaustive_one_res<5>();
  MESSAGE("cells enumerated exhaustively: ", total);
  REQUIRE(total == 122 + 842 + 5882 + 41162 + 288122 + 2016842);
}
