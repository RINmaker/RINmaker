#pragma once

#include "../src/kdtree/kdpoint.h"

template<size_t K>
class MyKDPoint : public kdpoint<K>
{
public:
    explicit MyKDPoint(std::array<double, K> const& pos) : kdpoint<K>(pos) {}
    MyKDPoint(MyKDPoint const& other) : kdpoint<K>(other) {}
};