#pragma once

#include <array>
#include <string>

#include "geometry.h"

// element of spatial
template<size_t K>
class kdpoint
{
protected:
    std::array<double, K> _position;

    explicit kdpoint(std::array<double, K> const& pos) : _position(pos)
    {}

    kdpoint(kdpoint const& other) : _position(other._position)
    {}

public:
    kdpoint() = delete;

    const double& operator[](size_t const axis) const
    { return _position[axis % K]; }

    explicit operator std::array<double, K> const&() const
    { return _position; }

    kdpoint<K>& operator=(kdpoint<K> const& rhs)
    {
        _position = (std::array<double, K>) rhs;
        return *this;
    }

    // coordinate-wise subtraction
    kdpoint<K> operator-(kdpoint<K> const& rhs) const
    {
        return kdpoint<K>(geom::difference<K>(_position, rhs._position));
    }

    // NOTE operator- would not be ideal, because it usually means coordinate-wise subtraction
    double distance(kdpoint<K> const& other) const
    {
        return geom::distance<K>(_position, other._position);
    }

    explicit operator std::string() const
    {
        std::string output;
        for (size_t i = 0; i < K; ++i)
        {
            output += std::to_string((*this)[i]);
            if (i + 1 < K)
            {
                output += ";";
            }
        }
        return output;
    }
};