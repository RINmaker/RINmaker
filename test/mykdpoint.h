#pragma once

#include <chrono>
#include <random>
#include <algorithm>
#include <gtest/gtest.h>

#include "utils/spatial/kdtree.h"
#include "utils/spatial/kdpoint.h"

template<size_t K>
class MyKDPoint : public kdpoint<K>
{
public:
    explicit MyKDPoint(std::array<double, K> const& pos) : kdpoint<K>(pos) {}
    MyKDPoint(MyKDPoint const& other) : kdpoint<K>(other) {}
};