// Concurrency. The library documents "give each thread its own Cache"; this
// checks that claim instead of asserting it.
//
// There are exactly two pieces of shared mutable state:
//
//   1. Cache, which the caller owns. Per-thread is safe. Sharing one is a real
//      data race, and TSan reports it on nanoh3.hpp's cache->face accesses if
//      you try. That is not simulated here, because a test suite should not
//      contain a deliberate race; see the comment at the bottom for how to
//      reproduce it.
//   2. the function-local static inside fast_axes(), initialised on first use.
//      C++11 guarantees that initialisation is thread-safe, but "guaranteed by
//      the standard" and "true in this build" are different claims, so every
//      thread is released from a spin barrier to make its first cell_fast()
//      call as close to simultaneous as possible.
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

TEST_CASE("nanoh3: 16 threads, private caches, simultaneous lazy-static init") {
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
    nanoh3::Cache fast_cache;
    for (int i = 0; i < 256; ++i) {
      const auto& e = ex[static_cast<std::size_t>(i) % ex.size()];
      nanoh3::Grid<11>::cell_fast(e.lat, e.lng, &fast_cache);
    }

    nanoh3::Cache mine;  // private, as documented
    long bad = 0;
    for (std::size_t n = 0; n < ex.size(); ++n) {
      const Expected& e = ex[(n + id * 997) % ex.size()];  // decorrelate threads
      if (nanoh3::Grid<11>::cell(e.lat, e.lng, &mine) != e.c11) ++bad;
      if (nanoh3::Grid<11>::cell(e.lat, e.lng) != e.c11) ++bad;
      if (nanoh3::Grid<7>::cell(e.lat, e.lng, &mine) != e.c7) ++bad;
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

// To confirm TSan is actually watching this file rather than silently inert,
// make `mine` above static and rebuild with -fsanitize=thread. It reports a
// data race on nanoh3.hpp's cache->face read and write. That was verified when
// this test was written: 4 reports with a shared cache, 0 with private caches.
