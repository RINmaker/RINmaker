#pragma once

#include <algorithm>
#include <vector>

#include "kdpoint.h"

template<class T, size_t K>
class kdtree {
public:
    class node;

private:
    // root node
    node *root;

public:
    // shallow copy O(1)
    kdtree &operator=(kdtree &&other) noexcept;

    // shallow copy O(1)
    kdtree(kdtree &&) noexcept;

    // balanced spatial O(n log^2 n)
    explicit kdtree(std::vector<T const *> &network);

    // empty spatial
    kdtree() : root(nullptr) {};

    // destruct everything O(n)
    ~kdtree();

    // range search in a (2*range)-edged cube O(log n)
    std::vector<T const *> range_search(kdpoint<K> const &test, double range) const;

    node *get_root() const {
        return root;
    }
};

template<class T, size_t K>
class kdtree<T, K>::node {
public:
    // what im holding
    T const *element;

    // left child
    node *left;

    // right child
    node *right;

    // balanced build O(n log^2 n)
    node(std::vector<T const *> &network, size_t first, size_t last, size_t depth);

    // destruct everything O(n)
    ~node();

    // range search in a (2*range)-edged cube O(log n)
    void range_search(kdpoint<K> const &test, double range, size_t depth, std::vector<T const *> &neighbors) const;
};

template<class T, size_t K>
kdtree<T, K> &kdtree<T, K>::operator=(kdtree<T, K> &&other) noexcept {
    if (this != &other) {
        root = other.root;
        other.root = nullptr;
    }

    return *this;
}

template<class T, size_t K>
kdtree<T, K>::kdtree(kdtree<T, K> &&other) noexcept : root(other.root) {
    other.root = nullptr;
}

template<class T, size_t K>
kdtree<T, K>::kdtree(std::vector<T const *> &network) {
    root = network.empty() ? nullptr : new node(network, 0, network.size() - 1, 0);
}

template<class T, size_t K>
kdtree<T, K>::~kdtree() {
    delete root;
}

template<class T, size_t K>
std::vector<T const *> kdtree<T, K>::range_search(kdpoint<K> const &test, double const range) const {
    std::vector<T const *> neighbors;

    if (root != nullptr) {
        root->range_search(test, range, 0, neighbors);
    }

    return neighbors;
}

template<class T, size_t K>
kdtree<T, K>::node::node(std::vector<T const *> &network, size_t const first, size_t const last, size_t const depth) {
    std::stable_sort(
            network.begin() + first, network.begin() + last + 1, [depth](kdpoint<K> const *a, kdpoint<K> const *b) {
                return (*a)[depth] < (*b)[depth];
            }
    );

    size_t const median = (first + last) / 2;

    element = network[median];

    left = median - first == 0 ? nullptr : new node(network, first, median - 1, depth + 1);
    right = last - median == 0 ? nullptr : new node(network, median + 1, last, depth + 1);
}

template<class T, size_t K>
kdtree<T, K>::node::~node() {
    delete left;
    delete right;
}

template<class T, size_t K>
void kdtree<T, K>::node::range_search(kdpoint<K> const &test, double const range, size_t const depth,
                                      std::vector<T const *> &neighbors) const {
    if ((*element)[depth] < test[depth] - range) {
        if (right != nullptr) {
            right->range_search(test, range, depth + 1, neighbors);
        }
    } else if ((*element)[depth] > test[depth] + range) {
        if (left != nullptr) {
            left->range_search(test, range, depth + 1, neighbors);
        }
    } else {
        if (left != nullptr) {
            left->range_search(test, range, depth + 1, neighbors);
        }

        if (right != nullptr) {
            right->range_search(test, range, depth + 1, neighbors);
        }

        if (element->distance(test) <= range) {
            neighbors.push_back(element);
        }
    }
}