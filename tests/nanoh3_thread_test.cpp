// Concurrency. Since the face cache was removed there is exactly ONE piece of
// shared mutable state left in the library: the function-local static inside
// fast_axes(), initialised on first use and read-only thereafter.
//
// C++11 guarantees that initialisation is thread-safe. "Guaranteed by the
// standard" and "true in this build" are different claims, so every thread is
// released from a spin barrier to make its first cell_fast() call as close to
// simultaneous as possible, which is the only moment that static is ever
// written.
//
// Everything else is a pure function of its arguments, which is the whole
// reason this file is short. It used to also police per-thread cache ownership;
// deleting the cache deleted that obligation.
//
// No H3 function is called inside a thread. Expected values are precomputed
// single-threaded, so under TSan any report is unambiguously about this
// library rather than about h3, which is linked uninstrumented.
#include <doctest/doctest.h>
#include <h3api.h>

#include <atomic>
#include <cmath>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

#include "nanoh3/nanoh3.hpp"

namespace {

struct Expected {
  double lat, lng;
  std::uint64_t c11, c7;
  double clat, clng;
  std::uint64_t ring[7];
};

}  // namespace

TEST_CASE("nanoh3: 16 threads, pure entry points, simultaneous lazy-static init") {
  const unsigned kThreads = 16;
  const int kPoints = 20'000;

  std::vector<Expected> ex(kPoints);
  {
    std::mt19937_64 rng(20260817);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (auto& e : ex) {
      e.lat = std::asin(2.0 * u(rng) - 1.0);
      e.lng = (2.0 * u(rng) - 1.0) * M_PI;
      LatLng g{e.lat, e.lng};
      H3Index a = 0, b = 0;
      REQUIRE(latLngToCell(&g, 11, &a) == E_SUCCESS);
      REQUIRE(latLngToCell(&g, 7, &b) == E_SUCCESS);
      e.c11 = a;
      e.c7 = b;
      LatLng c;
      REQUIRE(cellToLatLng(a, &c) == E_SUCCESS);
      e.clat = c.lat;
      e.clng = c.lng;
      H3Index r[7] = {};
      REQUIRE(gridDisk(a, 1, r) == E_SUCCESS);
      for (int i = 0; i < 7; ++i) e.ring[i] = r[i];
    }
  }

  std::atomic<unsigned> at_gate{0};
  std::atomic<bool> go{false};
  std::atomic<long> mismatches{0};
  std::atomic<long> iterations{0};

  auto work = [&](unsigned id) {
    // Spin barrier rather than std::barrier: this library is C++17.
    at_gate.fetch_add(1, std::memory_order_acq_rel);
    while (!go.load(std::memory_order_acquire)) { /* spin */ }

    // First touch of fast_axes()'s static, as simultaneous as we can arrange.
    for (int i = 0; i < 256; ++i) {
      const auto& e = ex[static_cast<std::size_t>(i) % ex.size()];
      nanoh3::Grid<11>::cell_fast(e.lat, e.lng);
    }

    long bad = 0;
    for (std::size_t n = 0; n < ex.size(); ++n) {
      const Expected& e = ex[(n + id * 997) % ex.size()];  // decorrelate threads
      if (nanoh3::Grid<11>::cell(e.lat, e.lng) != e.c11) ++bad;
      if (nanoh3::Grid<7>::cell(e.lat, e.lng) != e.c7) ++bad;
      double lat = 0, lng = 0;
      nanoh3::Grid<11>::center(e.c11, lat, lng);
      if (lat != e.clat || lng != e.clng) ++bad;
      std::uint64_t r[7] = {};
      nanoh3::Grid<11>::ring1(e.c11, r);
      for (int i = 0; i < 7; ++i) if (r[i] != e.ring[i]) ++bad;
      iterations.fetch_add(1, std::memory_order_relaxed);
    }
    mismatches.fetch_add(bad, std::memory_order_relaxed);
  };

  std::vector<std::thread> threads;
  threads.reserve(kThreads);
  for (unsigned i = 0; i < kThreads; ++i) threads.emplace_back(work, i);
  while (at_gate.load(std::memory_order_acquire) < kThreads) { /* spin */ }
  go.store(true, std::memory_order_release);
  for (auto& t : threads) t.join();

  MESSAGE("threads: ", kThreads, ", iterations: ", iterations.load());
  REQUIRE(mismatches.load() == 0);
}

// TSan was confirmed non-vacuous on this harness while the face cache still
// existed: making the per-thread Cache static produced 4 data-race reports on
// nanoh3.hpp's cache->face accesses, against 0 with private caches. With the
// cache gone there is no longer any way to provoke a race from the public API,
// which is the point.
