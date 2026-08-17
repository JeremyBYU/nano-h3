// nanoh3 vs H3 latLngToCell, in the two regimes that matter:
//  - uniform: independent points, no locality (nanoh3's floor — the face
//    cache rarely helps, the integer digit walk is the whole win);
//  - trace: a 15 m-step random walk (the matcher's regime — face cache ~100%).
// Points are precomputed; each op is exactly one conversion.
#include <h3api.h>
#include <nanobench.h>

#include <cmath>
#include <random>
#include <vector>

#include "nanoh3/nanoh3.hpp"

int main() {
  std::mt19937_64 rng(20260814);
  std::uniform_real_distribution<double> u(0.0, 1.0);

  std::vector<LatLng> uniform;
  for (int i = 0; i < 65536; ++i) {
    uniform.push_back({std::asin(2.0 * u(rng) - 1.0), (2.0 * u(rng) - 1.0) * M_PI});
  }
  std::vector<LatLng> trace;
  {
    constexpr double kD2R = M_PI / 180.0;
    double lat = 45.52, lng = -122.68;  // Portland
    std::normal_distribution<double> step(0.0, 15.0 / 111000.0);
    for (int i = 0; i < 65536; ++i) {
      lat += step(rng);
      lng += step(rng);
      trace.push_back({lat * kD2R, lng * kD2R});
    }
  }

  std::vector<LatLng> oneface;
  {
    constexpr double kD2R = M_PI / 180.0;
    std::uniform_real_distribution<double> dlat(44.0, 47.0), dlng(-124.0, -120.0);
    for (int i = 0; i < 65536; ++i) {
      oneface.push_back({dlat(rng) * kD2R, dlng(rng) * kD2R});
    }
  }

  ankerl::nanobench::Bench b;
  b.title("lat/lng -> res-11 cell").unit("pt").warmup(3).minEpochIterations(50);

  std::size_t k = 0;
  b.run("h3 latLngToCell, uniform", [&] {
    H3Index out = 0;
    latLngToCell(&uniform[k++ & 65535], 11, &out);
    ankerl::nanobench::doNotOptimizeAway(out);
  });
  k = 0;
  b.run("nanoh3, uniform, no cache", [&] {
    const auto& p = uniform[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell(p.lat, p.lng));
  });
  k = 0;
  nanoh3::Cache cu;
  b.run("nanoh3, uniform, cache", [&] {
    const auto& p = uniform[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell(p.lat, p.lng, &cu));
  });
  k = 0;
  b.run("h3 latLngToCell, trace walk", [&] {
    H3Index out = 0;
    latLngToCell(&trace[k++ & 65535], 11, &out);
    ankerl::nanobench::doNotOptimizeAway(out);
  });
  k = 0;
  nanoh3::Cache ct;
  b.run("nanoh3, trace walk, cache", [&] {
    const auto& p = trace[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell(p.lat, p.lng, &ct));
  });
  // The row whose absence hid the fact that Cache is worth ~0%. Without a
  // trace/no-cache pair there is no way to tell whether the trace speedup comes
  // from the cache or from the hardware exploiting the same locality for free.
  k = 0;
  b.run("nanoh3, trace walk, no cache", [&] {
    const auto& p = trace[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell(p.lat, p.lng));
  });
  // Scattered points confined to ONE face: separates "same face" locality from
  // "consecutive points" locality. Both are free; neither needs a Cache.
  k = 0;
  b.run("nanoh3, one face, scattered, no cache", [&] {
    const auto& p = oneface[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell(p.lat, p.lng));
  });
  k = 0;
  nanoh3::Cache cf;
  b.run("nanoh3 cell_fast, trace walk, cache", [&] {
    const auto& p = trace[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell_fast(p.lat, p.lng, &cf));
  });
  k = 0;
  nanoh3::Cache cfu;
  b.run("nanoh3 cell_fast, uniform, cache", [&] {
    const auto& p = uniform[k++ & 65535];
    ankerl::nanobench::doNotOptimizeAway(nanoh3::Grid<11>::cell_fast(p.lat, p.lng, &cfu));
  });
  return 0;
}
