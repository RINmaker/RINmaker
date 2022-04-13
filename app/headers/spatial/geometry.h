#pragma once

#include <vector>
#include <array>
#include <cmath>

namespace geom
{
constexpr double PI_GRECO = 3.14159265358979323846;

// Centroid of points
template <size_t K>
std::array<double, K> centroid(std::vector<std::array<double, K> const> const& points)
{
    std::array<double, K> sum({});

    for (auto const& p : points)
        for (size_t i = 0; i < K; ++i)
            sum[i] += p[i];

    double num = points.size();
    for (size_t i = 0; i < K; ++i)
        sum[i] /= num;

    return sum;
}

// Vector sum: v+w
template <size_t K>
std::array<double, K> sum(std::array<double, K> const& v, std::array<double, K> const& w)
{
    std::array<double, K> sum({});
    for (size_t i = 0; i < K; ++i)
        sum[i] = v[i] + w[i];

    return sum;
}

// Vector difference: v-w
template <size_t K>
std::array<double, K> difference(std::array<double, K> const& v, std::array<double, K> const& w)
{
    std::array<double, K> diff({});
    for (size_t i = 0; i < K; ++i)
        diff[i] = v[i] - w[i];

    return diff;
}

// Scalar product between vectors
template <size_t K>
double dot(std::array<double, K> const& v, std::array<double, K> const& w)
{
    double sum = 0;
    for (size_t i = 0; i < K; ++i)
        sum += v[i] * w[i];

    return sum;
}

// Modulus of a vector
template <size_t K>
double magnitude(std::array<double, K> const& v)
{
    return sqrt(dot(v, v));
}

// Normalize a vector
template <size_t K>
std::array<double, K> normalize(std::array<double, K> const& v)
{
    std::array<double, K> u(v);

    double m = magnitude(v);
    for (size_t i = 0; i < K; ++i)
        u[i] /= m;

    return u;
}

// Distance between two points
template <size_t K>
double distance(std::array<double, K> const& a, std::array<double, K> const& b)
{
    return magnitude(difference(a, b));
}

// Angle between two vectors in degrees
template <size_t K>
double angle(std::array<double, K> const& v, std::array<double, K> const& w)
{
    //TODO segnala errore se uno dei sue vettori è nullo
    return acos(dot(v, w) / (magnitude(v) * magnitude(w))) * (180.0 / PI_GRECO);
}

// Angle between the directions of two vectors in degrees ° (Returns the smaller angle)
template <size_t K>
double d_angle(std::array<double, K> const& v, std::array<double, K > const& w)
{
    double a = angle(v, w);
    return a > 90 ? 180 - a : a;
}

// The vector product is mathematically defined only between vectors in R^3
inline static std::array<double, 3> cross(std::array<double, 3> const& v, std::array<double, 3> const& w)
{
    return std::array<double, 3>({
        (v[1] * w[2]) - (v[2] * w[1]),
        (v[2] * w[0]) - (v[0] * w[2]),
        (v[0] * w[1]) - (v[1] * w[0])
        });
}
}
