#pragma once
/*
 * Portions of this file are transcribed from Uber H3 v4.1.0, whose upstream
 * copyright notices are retained here as required by Apache-2.0 Section 4(c):
 *
 *   src/h3lib/lib/faceijk.c   Copyright 2016-2021 Uber Technologies, Inc.
 *   src/h3lib/lib/coordijk.c  Copyright 2016-2018, 2020-2022 Uber Technologies, Inc.
 *   src/h3lib/lib/h3Index.c   Copyright 2016-2021 Uber Technologies, Inc.
 *   src/h3lib/lib/latLng.c    Copyright 2016-2021 Uber Technologies, Inc.
 *   src/h3lib/lib/algos.c     Copyright 2016-2021 Uber Technologies, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
///////////////////////////////////////////////////////////////////////////////
// MODIFICATIONS made by Jeremy Castagno, Copyright 2026
// Licensed under the Apache License, Version 2.0 (see LICENSE).
// This file has been changed from the original H3 v4.1.0 sources.
//
// Summary of Modifications:
//   Transcribed the latLngToCell / cellToLatLng / gridDisk(k=1) call chains
//   from H3's faceijk.c, coordijk.c, h3Index.c, latLng.c and algos.c into one
//   dependency-free C++17 header, then:
//     - templated the resolution so res-dependent branches fold at compile time
//     - replaced _upAp7/_upAp7r's lroundl((3i-j)/7.0L) with exact integer
//       round-half-away division
//     - added cell_fast, a vector gnomonic projection with no H3 counterpart
//       (NOT bit-identical to H3)
//
// H3 is a project of Uber Technologies, Inc. nano-h3 is an independent
// derivative work and is not affiliated with, endorsed by, or supported by
// Uber Technologies, Inc.
///////////////////////////////////////////////////////////////////////////////
//
// nanoh3 — a specialized, self-contained lat/lng -> H3 cell converter for
// ONE compile-time resolution, exploiting the two things a map matcher's
// query stream guarantees that the general H3 API cannot assume:
//
//   1. FIXED RESOLUTION. The res-dependent branches and scale loop fold at
//      compile time (Grid<11> for this project; any 0..15 works).
//   2. SPATIAL LOCALITY, which the hardware exploits for free. Consecutive
//      points share an icosahedron face essentially always, so the 20-face
//      argmin, the hex2d rounding branch tree and the digit walk all take the
//      same paths repeatedly: branch mispredictions fall from 6.9% to 1.2% and
//      the trace regime runs 113 ns/cell faster than scattered global points.
//      An explicit single-face cache USED to sit here. It was removed: it hit
//      99.9999% of the time and was still worth between -1.1% and +2.8%,
//      because twenty unrollable distance computations vectorise into almost
//      nothing. Removing it also removed the library's only mutable state, so
//      every entry point below is now a pure function.
//
// BIT-IDENTITY CONTRACT: Grid<R>::cell(), ::center() and ::ring1() return
// exactly what H3 v4.1.0's latLngToCell(), cellToLatLng() and gridDisk(k=1)
// return. The two halves of that claim rest on different evidence, and the
// difference is worth stating rather than blurring:
//
//   The integer back end (the aperture-7 digit walk) is PROVEN. It replaces
//   lroundl((3i-j)/7.0L) with exact integer round-half-away division: for
//   integer n, n/7 is never a half-integer, its distance from one is >= 1/14,
//   and long double rounds n/7 to within 2^-63, so both round to the same
//   integer, always.
//
//   The floating-point front end (face search, gnomonic projection, hex2d
//   cube-round) is TESTED, not proven. It replicates H3's operation sequence
//   and precision statement by statement, and the differential suite in
//   tests/ enforces equality over 41.8 million assertions spanning uniform,
//   regional, boundary-adversarial, pentagon, digit-walk and cell-vertex
//   inputs, at every resolution from 0 to 15. That is strong evidence, not a
//   proof over every representable double.
//
// cell_fast() is explicitly OUTSIDE this contract and is NOT bit-identical to
// H3; see its own comment.
//
// Derived from Uber H3 v4.1.0 (Apache License 2.0). See NOTICE for which
// parts derive from what. This header depends on nothing but <cmath> and
// <cstdint>; H3 itself is linked only by the tests.

// Bump on any behavioural change. Consumers that pin a version can turn a
// silently-swapped header into a compile error with a static_assert on these.
#define NANOH3_VERSION_MAJOR 0
#define NANOH3_VERSION_MINOR 2
#define NANOH3_VERSION_PATCH 0

#include <cmath>
#include <cstdint>

namespace nanoh3 {

struct LatLngRad {
  double lat, lng;
};
struct Vec3 {
  double x, y, z;
};
struct BaseCellRot {
  int baseCell, ccwRot60;
};
struct CoordIJK {
  int i, j, k;
};
struct FaceOrient {
  int face;
  CoordIJK translate;
  int ccwRot60;
};

// The .inl resolves through the include search path, not relative to this
// file, so a stale copy of nanoh3_tables.inl elsewhere on the path can shadow
// just the tables while this header comes from the intended checkout. That
// split compiles silently and the header's own version macros do not catch it,
// because they live here. Hence a separate version stamped inside the .inl.
#include "nanoh3/nanoh3_tables.inl"
#ifndef NANOH3_TABLES_VERSION
#error "nanoh3_tables.inl did not define NANOH3_TABLES_VERSION: a stale or foreign copy is shadowing it on the include path"
#endif
static_assert(NANOH3_TABLES_VERSION == 1,
              "nanoh3_tables.inl version does not match this nanoh3.hpp");

// H3 constants, verbatim (long-double literals matter: the front end must
// reproduce H3's exact mixed-precision arithmetic).
inline constexpr long double kEpsilon = 0.0000000000000001L;
inline constexpr long double kM2Pi = 6.28318530717958647692528676655900576839433L;
inline constexpr long double kSin60 = 0.8660254037844386467637231707529361834714L;
inline constexpr long double kAp7RotRads = 0.333473172251832115336090755351601070065900389L;
inline constexpr long double kRes0UGnomonic = 0.38196601125010500003L;
inline constexpr long double kSqrt7 = 2.6457513110645905905016157536392604257102L;

// H3 index bit layout (h3Index.h).
inline constexpr std::uint64_t kH3Init = 35184372088831ull;  // res 15, all digit-7
inline constexpr int kModeOffset = 59, kResOffset = 52, kBcOffset = 45, kPerDigit = 3;
inline constexpr int kCellMode = 1, kMaxRes = 15;
inline constexpr int kKAxesDigit = 1, kInvalidDigit = 7;

namespace detail {

inline int digit_at(std::uint64_t h, int r) {
  return static_cast<int>((h >> ((kMaxRes - r) * kPerDigit)) & 7u);
}
inline std::uint64_t set_digit(std::uint64_t h, int r, int d) {
  const int off = (kMaxRes - r) * kPerDigit;
  return (h & ~(7ull << off)) | (static_cast<std::uint64_t>(d) << off);
}

// Exact integer round-half-away division by 7 (== lroundl(n / 7.0L) for all
// int n reachable here; see the header contract note).
inline int divround7(int n) {
  return n >= 0 ? (n + 3) / 7 : -((-n + 3) / 7);
}

inline void normalize(CoordIJK& c);

inline void up_ap7(CoordIJK& c) {  // integer-exact _upAp7 (incl. its normalize)
  const int i = c.i - c.k, j = c.j - c.k;
  c.i = divround7(3 * i - j);
  c.j = divround7(i + 2 * j);
  c.k = 0;
  normalize(c);
}
inline void up_ap7r(CoordIJK& c) {  // integer-exact _upAp7r (incl. its normalize)
  const int i = c.i - c.k, j = c.j - c.k;
  c.i = divround7(2 * i + j);
  c.j = divround7(3 * j - i);
  c.k = 0;
  normalize(c);
}

inline void normalize(CoordIJK& c) {
  if (c.i < 0) { c.j -= c.i; c.k -= c.i; c.i = 0; }
  if (c.j < 0) { c.i -= c.j; c.k -= c.j; c.j = 0; }
  if (c.k < 0) { c.i -= c.k; c.j -= c.k; c.k = 0; }
  int m = c.i;
  if (c.j < m) m = c.j;
  if (c.k < m) m = c.k;
  if (m > 0) { c.i -= m; c.j -= m; c.k -= m; }
}

inline void down_ap7(CoordIJK& c) {
  const CoordIJK r{3 * c.i + 1 * c.j + 0 * c.k,
                   0 * c.i + 3 * c.j + 1 * c.k,
                   1 * c.i + 0 * c.j + 3 * c.k};
  c = r;
  normalize(c);
}
inline void down_ap7r(CoordIJK& c) {
  const CoordIJK r{3 * c.i + 0 * c.j + 1 * c.k,
                   1 * c.i + 3 * c.j + 0 * c.k,
                   0 * c.i + 1 * c.j + 3 * c.k};
  c = r;
  normalize(c);
}

inline int unit_ijk_to_digit(const CoordIJK& in) {
  constexpr CoordIJK kUnit[7] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                                 {1, 0, 0}, {1, 0, 1}, {1, 1, 0}};
  CoordIJK c = in;
  normalize(c);
  for (int d = 0; d < 7; ++d) {
    if (c.i == kUnit[d].i && c.j == kUnit[d].j && c.k == kUnit[d].k) return d;
  }
  return kInvalidDigit;
}

// Digit rotations (coordijk.c _rotate60ccw/_rotate60cw):
// ccw K->IK->I->IJ->J->JK->K, cw is the inverse cycle.
inline int rotate60ccw(int d) {
  constexpr int t[8] = {0, 5, 3, 1, 6, 4, 2, 7};
  return t[d];
}
inline int rotate60cw(int d) {
  constexpr int t[8] = {0, 3, 6, 2, 5, 1, 4, 7};
  return t[d];
}

inline int leading_nonzero_digit(std::uint64_t h, int res) {
  for (int r = 1; r <= res; ++r) {
    if (digit_at(h, r) != 0) return digit_at(h, r);
  }
  return 0;
}

inline std::uint64_t rotate60ccw_index(std::uint64_t h, int res) {
  for (int r = 1; r <= res; ++r) h = set_digit(h, r, rotate60ccw(digit_at(h, r)));
  return h;
}
inline std::uint64_t rotate60cw_index(std::uint64_t h, int res) {
  for (int r = 1; r <= res; ++r) h = set_digit(h, r, rotate60cw(digit_at(h, r)));
  return h;
}
inline std::uint64_t rotate_pent60ccw_index(std::uint64_t h, int res) {
  bool found = false;
  for (int r = 1; r <= res; ++r) {
    h = set_digit(h, r, rotate60ccw(digit_at(h, r)));
    if (!found && digit_at(h, r) != 0) {
      found = true;
      if (leading_nonzero_digit(h, res) == kKAxesDigit) h = rotate60ccw_index(h, res);
    }
  }
  return h;
}

// _hex2dToCoordIJK, verbatim (mixed double/long-double exactly as H3).
inline void hex2d_to_ijk(double vx, double vy, CoordIJK& h) {
  double a1, a2, x1, x2, r1, r2;
  int m1, m2;
  h.k = 0;
  a1 = fabsl(vx);
  a2 = fabsl(vy);
  x2 = a2 / kSin60;
  x1 = a1 + x2 / 2.0L;
  m1 = static_cast<int>(x1);
  m2 = static_cast<int>(x2);
  r1 = x1 - m1;
  r2 = x2 - m2;
  if (r1 < 0.5L) {
    if (r1 < 1.0L / 3.0L) {
      if (r2 < (1.0L + r1) / 2.0L) { h.i = m1; h.j = m2; }
      else { h.i = m1; h.j = m2 + 1; }
    } else {
      if (r2 < (1.0L - r1)) h.j = m2;
      else h.j = m2 + 1;
      if ((1.0L - r1) <= r2 && r2 < (2.0 * r1)) h.i = m1 + 1;
      else h.i = m1;
    }
  } else {
    if (r1 < 2.0L / 3.0L) {
      if (r2 < (1.0L - r1)) h.j = m2;
      else h.j = m2 + 1;
      if ((2.0L * r1 - 1.0L) < r2 && r2 < (1.0L - r1)) h.i = m1;
      else h.i = m1 + 1;
    } else {
      if (r2 < (r1 / 2.0L)) { h.i = m1 + 1; h.j = m2; }
      else { h.i = m1 + 1; h.j = m2 + 1; }
    }
  }
  if (vx < 0.0L) {
    if ((h.j % 2) == 0) {  // even
      const long long axisi = h.j / 2;
      const long long diff = h.i - axisi;
      h.i = static_cast<int>(h.i - 2.0 * diff);
    } else {
      const long long axisi = (h.j + 1) / 2;
      const long long diff = h.i - axisi;
      h.i = static_cast<int>(h.i - (2.0 * diff + 1));
    }
  }
  if (vy < 0.0L) {
    h.i = h.i - (2 * h.j + 1) / 2;
    h.j = -1 * h.j;
  }
  normalize(h);
}

inline double pos_angle(double rads) {
  double tmp = ((rads < 0.0L) ? rads + kM2Pi : rads);
  if (rads >= kM2Pi) tmp -= kM2Pi;
  return tmp;
}

inline double geo_azimuth(const LatLngRad& p1, double lat, double lng) {
  return atan2(cos(lat) * sin(lng - p1.lng),
               cos(p1.lat) * sin(lat) - sin(p1.lat) * cos(lat) * cos(lng - p1.lng));
}

}  // namespace detail

namespace detail {

inline bool is_pentagon_index(std::uint64_t h, int res) {
  const int bc = static_cast<int>((h >> kBcOffset) & 127u);
  return kIsPentagon[bc] && leading_nonzero_digit(h, res) == 0;
}

// h3NeighborRotations (algos.c), verbatim logic on the ported tables.
// Returns 0 on E_PENTAGON (deleted k direction), which the callers treat
// exactly as H3 does.
inline std::uint64_t neighbor_rotations(std::uint64_t origin, int dir,
                                        int& rotations, int res) {
  std::uint64_t current = origin;
  rotations = rotations % 6;
  for (int i = 0; i < rotations; ++i) dir = rotate60ccw(dir);

  int new_rotations = 0;
  const int old_bc = static_cast<int>((current >> kBcOffset) & 127u);
  const int old_leading = leading_nonzero_digit(current, res);

  int r = res - 1;
  for (;;) {
    if (r == -1) {
      int nbc = kBaseCellNeighbors[old_bc][dir];
      new_rotations = kBaseCellNeighbor60CCWRots[old_bc][dir];
      if (nbc == 127) {  // INVALID_BASE_CELL: deleted k vertex at base level
        nbc = kBaseCellNeighbors[old_bc][5 /*IK*/];
        new_rotations = kBaseCellNeighbor60CCWRots[old_bc][5];
        current = rotate60ccw_index(current, res);
        rotations = rotations + 1;
      }
      current = (current & ~(127ull << kBcOffset)) |
                (static_cast<std::uint64_t>(nbc) << kBcOffset);
      break;
    }
    const int old_digit = digit_at(current, r + 1);
    int next_dir;
    if (old_digit == kInvalidDigit) return 0;
    if ((r + 1) % 2 == 1) {  // Class III res -> "II" tables per H3's naming
      current = set_digit(current, r + 1, kNEW_DIGIT_II[old_digit][dir]);
      next_dir = kNEW_ADJUSTMENT_II[old_digit][dir];
    } else {
      current = set_digit(current, r + 1, kNEW_DIGIT_III[old_digit][dir]);
      next_dir = kNEW_ADJUSTMENT_III[old_digit][dir];
    }
    if (next_dir != 0) {
      dir = next_dir;
      --r;
    } else {
      break;
    }
  }

  const int new_bc = static_cast<int>((current >> kBcOffset) & 127u);
  if (kIsPentagon[new_bc]) {
    bool adjusted_k = false;
    if (leading_nonzero_digit(current, res) == kKAxesDigit) {
      if (old_bc != new_bc) {
        const bool cw = kCwOffsetPent[new_bc][0] == kHomeFace[old_bc] ||
                        kCwOffsetPent[new_bc][1] == kHomeFace[old_bc];
        current = cw ? rotate60cw_index(current, res)
                     : rotate60ccw_index(current, res);
        adjusted_k = true;
      } else {
        if (old_leading == 0) return 0;  // E_PENTAGON: k deleted from here
        if (old_leading == 3 /*JK*/) {
          current = rotate60ccw_index(current, res);
          rotations = rotations + 1;
        } else if (old_leading == 5 /*IK*/) {
          current = rotate60cw_index(current, res);
          rotations = rotations + 5;
        } else {
          return 0;
        }
      }
    }
    for (int i = 0; i < new_rotations; ++i)
      current = rotate_pent60ccw_index(current, res);
    if (old_bc != new_bc) {
      if (new_bc == 4 || new_bc == 117) {  // polar pentagons
        if (old_bc != 118 && old_bc != 8 &&
            leading_nonzero_digit(current, res) != 3 /*JK*/) {
          rotations = rotations + 1;
        }
      } else if (leading_nonzero_digit(current, res) == 5 /*IK*/ &&
                 !adjusted_k) {
        rotations = rotations + 1;
      }
    }
  } else {
    for (int i = 0; i < new_rotations; ++i)
      current = rotate60ccw_index(current, res);
  }
  rotations = (rotations + new_rotations) % 6;
  return current;
}

inline void ijk_rotate60ccw(CoordIJK& c) {
  const CoordIJK r{c.i * 1 + c.j * 0 + c.k * 1, c.i * 1 + c.j * 1 + c.k * 0,
                   c.i * 0 + c.j * 1 + c.k * 1};
  c = r;
  normalize(c);
}
inline void ijk_rotate60cw(CoordIJK& c) {
  const CoordIJK r{c.i * 1 + c.j * 1 + c.k * 0, c.i * 0 + c.j * 1 + c.k * 1,
                   c.i * 1 + c.j * 0 + c.k * 1};
  c = r;
  normalize(c);
}

inline void neighbor_digit(CoordIJK& c, int digit) {
  constexpr CoordIJK kUnit[7] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                                 {1, 0, 0}, {1, 0, 1}, {1, 1, 0}};
  if (digit > 0 && digit < 7) {
    c.i += kUnit[digit].i;
    c.j += kUnit[digit].j;
    c.k += kUnit[digit].k;
    normalize(c);
  }
}

// _ijkToHex2d, verbatim mixed precision (0.5L, M_SQRT3_2 long double).
inline void ijk_to_hex2d(const CoordIJK& h, double& x, double& y) {
  const int i = h.i - h.k, j = h.j - h.k;
  x = static_cast<double>(i - 0.5L * j);
  y = static_cast<double>(j * kSin60);
}

inline double constrain_lng(double lng) {
  while (lng > M_PI) lng = lng - (2 * M_PI);
  while (lng < -M_PI) lng = lng + (2 * M_PI);
  return lng;
}

// _geoAzDistanceRads, verbatim.
inline void geo_az_distance(const LatLngRad& p1, double az, double distance,
                            double& lat, double& lng) {
  if (distance < kEpsilon) {
    lat = p1.lat;
    lng = p1.lng;
    return;
  }
  az = pos_angle(az);
  if (az < kEpsilon || fabs(az - M_PI) < kEpsilon) {
    lat = az < kEpsilon ? p1.lat + distance : p1.lat - distance;
    if (fabs(lat - M_PI_2) < kEpsilon) {
      lat = M_PI_2;
      lng = 0.0;
    } else if (fabs(lat + M_PI_2) < kEpsilon) {
      lat = -M_PI_2;
      lng = 0.0;
    } else {
      lng = constrain_lng(p1.lng);
    }
  } else {
    double sinlat = sin(p1.lat) * cos(distance) +
                    cos(p1.lat) * sin(distance) * cos(az);
    if (sinlat > 1.0) sinlat = 1.0;
    if (sinlat < -1.0) sinlat = -1.0;
    lat = asin(sinlat);
    if (fabs(lat - M_PI_2) < kEpsilon) {
      lat = M_PI_2;
      lng = 0.0;
    } else if (fabs(lat + M_PI_2) < kEpsilon) {
      lat = -M_PI_2;
      lng = 0.0;
    } else {
      double sinlng = sin(az) * sin(distance) / cos(lat);
      double coslng = (cos(distance) - sin(p1.lat) * sin(lat)) / cos(p1.lat) / cos(lat);
      if (sinlng > 1.0) sinlng = 1.0;
      if (sinlng < -1.0) sinlng = -1.0;
      if (coslng > 1.0) coslng = 1.0;
      if (coslng < -1.0) coslng = -1.0;
      lng = constrain_lng(p1.lng + atan2(sinlng, coslng));
    }
  }
}

// _adjustOverageClassII, substrate = 0 (cell centers only).
struct FaceIJK {
  int face;
  CoordIJK coord;
};
inline bool adjust_overage_class2(FaceIJK& f, int res, bool pent_leading4) {
  CoordIJK& ijk = f.coord;
  const int max_dim = k_maxDimByCIIres[res];
  if (ijk.i + ijk.j + ijk.k <= max_dim) return false;  // NO_OVERAGE
  const FaceOrient* orient;
  if (ijk.k > 0) {
    if (ijk.j > 0) {  // jk quadrant
      orient = &kFaceNeighbors[f.face][3 /*JK*/];
    } else {  // ik quadrant
      orient = &kFaceNeighbors[f.face][2 /*KI*/];
      if (pent_leading4) {
        // translate origin to pentagon center, rotate cw, translate back
        CoordIJK tmp{ijk.i - max_dim, ijk.j, ijk.k};
        ijk_rotate60cw(tmp);
        ijk = {tmp.i + max_dim, tmp.j, tmp.k};
      }
    }
  } else {  // ij quadrant
    orient = &kFaceNeighbors[f.face][1 /*IJ*/];
  }
  f.face = orient->face;
  for (int i = 0; i < orient->ccwRot60; ++i) ijk_rotate60ccw(ijk);
  const int unit_scale = k_unitScaleByCIIres[res];
  ijk.i += orient->translate.i * unit_scale;
  ijk.j += orient->translate.j * unit_scale;
  ijk.k += orient->translate.k * unit_scale;
  normalize(ijk);
  return true;  // NEW_FACE
}

// _h3ToFaceIjk (+ WithInitializedFijk), verbatim logic.
inline void h3_to_face_ijk(std::uint64_t h, int res, FaceIJK& f) {
  int bc = static_cast<int>((h >> kBcOffset) & 127u);
  if (kIsPentagon[bc] && leading_nonzero_digit(h, res) == 5) {
    h = rotate60cw_index(h, res);
  }
  f.face = kHomeFace[bc];
  f.coord = {kHomeIJK[bc][0], kHomeIJK[bc][1], kHomeIJK[bc][2]};

  const bool possible_overage =
      kIsPentagon[bc] ||
      !(res == 0 || (f.coord.i == 0 && f.coord.j == 0 && f.coord.k == 0));
  for (int r = 1; r <= res; ++r) {
    if (r % 2 == 1) down_ap7(f.coord);  // Class III
    else down_ap7r(f.coord);
    neighbor_digit(f.coord, digit_at(h, r));
  }
  if (!possible_overage) return;

  const CoordIJK orig = f.coord;
  int adj_res = res;
  if (res % 2 == 1) {  // Class III: drop into the finer Class II grid
    down_ap7r(f.coord);
    ++adj_res;
  }
  const bool pent_leading4 = kIsPentagon[bc] && leading_nonzero_digit(h, res) == 4;
  if (adjust_overage_class2(f, adj_res, pent_leading4)) {
    if (kIsPentagon[bc]) {
      while (adjust_overage_class2(f, adj_res, false)) continue;
    }
    if (adj_res != res) up_ap7r(f.coord);
  } else if (adj_res != res) {
    f.coord = orig;
  }
}

}  // namespace detail

template <int Res>
class Grid {
  static_assert(Res >= 0 && Res <= 15, "H3 resolution");

 public:
  // lat/lng in RADIANS (like H3's internal LatLng).
  static std::uint64_t cell(double lat, double lng) {
    // ---- face + gnomonic hex2d (H3 _geoToHex2d, verbatim math) ----
    const double clat = cos(lat);
    const Vec3 p{cos(lng) * clat, sin(lng) * clat, sin(lat)};

    // Nearest of the 20 icosahedron face centres. This used to sit behind an
    // opt-in single-face cache, on the theory that consecutive points in a
    // track share a face and the search could be skipped. The theory was right
    // and the optimisation was worthless: the cache hit 99.9999% of the time
    // and bought between -1.1% and +2.8% depending on function and -O level.
    // Twenty unrollable distance computations vectorise into almost nothing,
    // and the branch predictor already exploits the locality the cache was
    // built to exploit. Deleting it removed the library's only mutable state,
    // which is why there is now nothing here to make thread-unsafe.
    int face = 0;
    double sqd = 5.0;
    for (int f = 0; f < 20; ++f) {
      const Vec3& c = kFaceCenterPoint[f];
      const double dx = c.x - p.x, dy = c.y - p.y, dz = c.z - p.z;
      const double d = dx * dx + dy * dy + dz * dz;
      if (d < sqd) {
        face = f;
        sqd = d;
      }
    }

    double r = acos(1 - sqd / 2);
    double vx = 0.0, vy = 0.0;
    if (r >= static_cast<double>(kEpsilon)) {
      double theta = detail::pos_angle(
          kFaceAxesAzRadsCII[face][0] -
          detail::pos_angle(detail::geo_azimuth(kFaceCenterGeo[face], lat, lng)));
      // Long double, folded exactly as H3 folds it: the constants stay long
      // double and each statement truncates back to double. This costs ~20%
      // via x87 and it is not optional.
      //
      // An earlier version truncated these three constants to double first and
      // claimed the differential suite had measured that bit-identical. It had
      // not. The suite samples random points and marches boundaries in 1 cm
      // steps, while the divergence band around a cell vertex is a few
      // nanometres wide, so nothing in it ever probed one. Feeding
      // cellToBoundary vertices straight back in exposed 42,253 mismatches out
      // of 192,727 probes (21.9%), at every resolution except 0. Every
      // mismatch was an adjacent cell. See the vertex test in tests/.
      //
      // If you want the faster arithmetic, that is what cell_fast() is for:
      // it is 4.1x rather than 2.2x and is explicitly outside the contract.
      // cell() is the one that must be exact.
      if constexpr (Res % 2 == 1) {  // Class III
        theta = detail::pos_angle(static_cast<double>(theta - kAp7RotRads));
      }
      r = tan(r);
      r = static_cast<double>(r / kRes0UGnomonic);
      for (int i = 0; i < Res; ++i) r = static_cast<double>(r * kSqrt7);
      vx = r * cos(theta);
      vy = r * sin(theta);
    }

    // ---- cube round + digit walk (_hex2dToCoordIJK + _faceIjkToH3) ----
    CoordIJK ijk;
    detail::hex2d_to_ijk(vx, vy, ijk);
    return assemble(face, ijk);
  }

 private:
  // The integer back end shared by cell() and cell_fast(): aperture-7 digit
  // walk, base-cell lookup, pentagon rotations.
  static std::uint64_t assemble(int face, CoordIJK ijk) {
    std::uint64_t h = kH3Init;
    h = (h & ~(15ull << kModeOffset)) | (static_cast<std::uint64_t>(kCellMode) << kModeOffset);
    h = (h & ~(15ull << kResOffset)) | (static_cast<std::uint64_t>(Res) << kResOffset);

    for (int rr = Res - 1; rr >= 0; --rr) {
      const CoordIJK last = ijk;
      CoordIJK center;
      if ((rr + 1) % 2 == 1) {  // res rr+1 is Class III
        detail::up_ap7(ijk);
        center = ijk;
        detail::down_ap7(center);
      } else {
        detail::up_ap7r(ijk);
        center = ijk;
        detail::down_ap7r(center);
      }
      CoordIJK diff{last.i - center.i, last.j - center.j, last.k - center.k};
      h = detail::set_digit(h, rr + 1, detail::unit_ijk_to_digit(diff));
    }

    if (ijk.i > 2 || ijk.j > 2 || ijk.k > 2) return 0;  // out of range

    const BaseCellRot& bcr = kFaceIjkBaseCells[face][ijk.i][ijk.j][ijk.k];
    h = (h & ~(127ull << kBcOffset)) | (static_cast<std::uint64_t>(bcr.baseCell) << kBcOffset);

    if (kIsPentagon[bcr.baseCell]) {
      if (detail::leading_nonzero_digit(h, Res) == kKAxesDigit) {
        const bool cw = kCwOffsetPent[bcr.baseCell][0] == face ||
                        kCwOffsetPent[bcr.baseCell][1] == face;
        h = cw ? detail::rotate60cw_index(h, Res) : detail::rotate60ccw_index(h, Res);
      }
      for (int i = 0; i < bcr.ccwRot60; ++i) h = detail::rotate_pent60ccw_index(h, Res);
    } else {
      for (int i = 0; i < bcr.ccwRot60; ++i) h = detail::rotate60ccw_index(h, Res);
    }
    return h;
  }

  // Scaled tangent axes per face for cell_fast: e1 along the face's CII
  // azimuth (minus the Class III rotation at odd Res: theta' = (az0-rot)-az),
  // e2 at azimuth-90deg, both scaled by sqrt7^Res / RES0_U_GNOMONIC. The
  // divergence test compares cell_fast against H3 over millions of points —
  // an axis-convention slip fails there at ~100% rate, a correct build
  // diverges only within float-noise of cell boundaries.
  struct FastAxes {
    Vec3 e1[20], e2[20];
    FastAxes() {
      long double scale = 1.0L / kRes0UGnomonic;
      for (int i = 0; i < Res; ++i) scale *= kSqrt7;
      for (int f = 0; f < 20; ++f) {
        const LatLngRad c = kFaceCenterGeo[f];
        double az0 = kFaceAxesAzRadsCII[f][0];
        if constexpr (Res % 2 == 1) az0 = az0 - static_cast<double>(kAp7RotRads);
        const double slat = sin(c.lat), clat_ = cos(c.lat);
        const double slng = sin(c.lng), clng = cos(c.lng);
        const Vec3 east{-slng, clng, 0.0};
        const Vec3 north{-slat * clng, -slat * slng, clat_};
        auto tangent = [&](double az) {
          return Vec3{north.x * cos(az) + east.x * sin(az),
                      north.y * cos(az) + east.y * sin(az),
                      north.z * cos(az) + east.z * sin(az)};
        };
        const Vec3 t1 = tangent(az0);
        const Vec3 t2 = tangent(az0 - M_PI_2);
        const double s = static_cast<double>(scale);
        e1[f] = {t1.x * s, t1.y * s, t1.z * s};
        e2[f] = {t2.x * s, t2.y * s, t2.z * s};
      }
    }
  };
  static const FastAxes& fast_axes() {
    static const FastAxes ax;
    return ax;
  }

 public:
  // Ring-1 (gridDisk k=1), replicating H3's facade exactly INCLUDING output
  // order: the unsafe spiral (origin, then J,JK,K,IK,I,IJ moves after the
  // initial I step), with the safe hash-set fallback (maxIdx=7, origin%7
  // probing) when a pentagon or its distortion area is met. out[7]; empty
  // slots are 0. Order matters: the matcher's candidate sort tie-breaks can
  // depend on cell visit order, and byte-identity is the contract.
  static void ring1(std::uint64_t origin, std::uint64_t out[7]) {
    for (int i = 0; i < 7; ++i) out[i] = 0;
    if (!unsafe_ring1(origin, out)) {
      for (int i = 0; i < 7; ++i) out[i] = 0;
      safe_disk1(origin, out, 1);
    }
  }

  // Cell center (cellToLatLng), radians. Verbatim _h3ToFaceIjk +
  // _faceIjkToGeo (inverse gnomonic through the face's polar frame).
  static void center(std::uint64_t cell, double& lat, double& lng) {
    detail::FaceIJK f;
    detail::h3_to_face_ijk(cell, Res, f);
    double vx, vy;
    detail::ijk_to_hex2d(f.coord, vx, vy);
    double r = sqrt(vx * vx + vy * vy);
    if (r < static_cast<double>(kEpsilon)) {
      lat = kFaceCenterGeo[f.face].lat;
      lng = kFaceCenterGeo[f.face].lng;
      return;
    }
    // Mixed precision VERBATIM: H3 folds the long-double constants into the
    // expression and truncates the long-double result back to double each
    // statement — replicating that is what bit-exactness means here.
    double theta = atan2(vy, vx);
    for (int i = 0; i < Res; ++i) r = static_cast<double>(r / kSqrt7);
    r = static_cast<double>(r * kRes0UGnomonic);
    r = atan(r);
    if constexpr (Res % 2 == 1) {  // Class III
      theta = detail::pos_angle(static_cast<double>(theta + kAp7RotRads));
    }
    theta = detail::pos_angle(kFaceAxesAzRadsCII[f.face][0] - theta);
    detail::geo_az_distance(kFaceCenterGeo[f.face], theta, r, lat, lng);
  }

  // OPT-IN fast path: pure vector gnomonic projection — hex2d comes from
  // three dot products and one divide against precomputed per-face axes
  // (azimuth frame + Class III rotation + resolution scale folded in),
  // replacing the exact path's acos/tan/atan2/sincos chain. NOT bit-identical
  // to H3: last-ulp hex2d differences can flip points that sit within
  // float-noise of a cell boundary (measured divergence rate in the test
  // suite; use only where a one-in-1e8 neighbor-cell assignment is
  // acceptable). The integer digit walk is shared with the exact path.
  static std::uint64_t cell_fast(double lat, double lng) {
    const double clat = cos(lat);
    const Vec3 p{cos(lng) * clat, sin(lng) * clat, sin(lat)};

    int face = 0;
    double sqd = 5.0;
    for (int f = 0; f < 20; ++f) {
      const Vec3& c = kFaceCenterPoint[f];
      const double dx = c.x - p.x, dy = c.y - p.y, dz = c.z - p.z;
      const double d = dx * dx + dy * dy + dz * dz;
      if (d < sqd) {
        face = f;
        sqd = d;
      }
    }

    const FastAxes& ax = fast_axes();
    const Vec3& c = kFaceCenterPoint[face];
    const double w = p.x * c.x + p.y * c.y + p.z * c.z;  // cos(angular dist)
    const Vec3& e1 = ax.e1[face];
    const Vec3& e2 = ax.e2[face];
    const double vx = (p.x * e1.x + p.y * e1.y + p.z * e1.z) / w;
    const double vy = (p.x * e2.x + p.y * e2.y + p.z * e2.z) / w;

    CoordIJK ijk;
    detail::hex2d_to_ijk(vx, vy, ijk);
    return assemble(face, ijk);
  }

  // Degree convenience (the matcher works in degrees). The long-double
  // constant and mixed multiply replicate H3's degsToRads bit-for-bit.
  static std::uint64_t cell_deg(double lat_deg, double lng_deg) {
    constexpr long double kPi180 = 0.0174532925199432957692369076848861271111L;
    return cell(static_cast<double>(lat_deg * kPi180),
                static_cast<double>(lng_deg * kPi180));
  }

 private:
  static bool unsafe_ring1(std::uint64_t origin, std::uint64_t out[7]) {
    // gridDiskDistancesUnsafe, k=1. DIRECTIONS = {J,JK,K,IK,I,IJ} = {2,3,1,5,4,6};
    // NEXT_RING_DIRECTION = I = 4.
    constexpr int kDirections[6] = {2, 3, 1, 5, 4, 6};
    int idx = 0;
    out[idx++] = origin;
    if (detail::is_pentagon_index(origin, Res)) return false;
    int rotations = 0;
    // enter ring 1 (not emitted)
    origin = detail::neighbor_rotations(origin, 4, rotations, Res);
    if (origin == 0) return false;
    if (detail::is_pentagon_index(origin, Res)) return false;
    for (int direction = 0; direction < 6; ++direction) {
      origin = detail::neighbor_rotations(origin, kDirections[direction],
                                          rotations, Res);
      if (origin == 0) return false;
      out[idx++] = origin;
      if (detail::is_pentagon_index(origin, Res)) return false;
    }
    return true;
  }

  // _gridDiskDistancesInternal with maxIdx = 7 (the k=1 safe fallback): a
  // 7-slot open-addressed hash set keyed by origin % 7, recursing to
  // neighbors with distance bookkeeping.
  static void safe_disk1(std::uint64_t origin, std::uint64_t out[7], int k,
                         int cur_k = 0, int* distances = nullptr) {
    int local_dist[7] = {0, 0, 0, 0, 0, 0, 0};
    if (!distances) distances = local_dist;
    int off = static_cast<int>(origin % 7);
    while (out[off] != 0 && out[off] != origin) off = (off + 1) % 7;
    if (out[off] == origin && distances[off] <= cur_k) return;
    out[off] = origin;
    distances[off] = cur_k;
    if (cur_k >= k) return;
    constexpr int kDirections[6] = {2, 3, 1, 5, 4, 6};
    for (int i = 0; i < 6; ++i) {
      int rotations = 0;
      const std::uint64_t nb =
          detail::neighbor_rotations(origin, kDirections[i], rotations, Res);
      if (nb == 0) continue;  // E_PENTAGON: expected when traversing off one
      safe_disk1(nb, out, k, cur_k + 1, distances);
    }
  }

};

}  // namespace nanoh3
