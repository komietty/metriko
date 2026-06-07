#ifndef METRIKO_TEST_CHECK_H
#define METRIKO_TEST_CHECK_H
#include <iostream>

// assert と違い NDEBUG でも無効化されない。失敗時に main から 1 を返す（CTest が失敗判定）。
#define CHECK(x) do { if (!(x)) { \
    std::cerr << "FAIL: " #x "  @ " << __FILE__ << ":" << __LINE__ << "\n"; \
    return 1; } } while (0)

#endif
