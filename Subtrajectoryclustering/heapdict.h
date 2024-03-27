#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <memory>
#ifndef HEAPDICT_H
#define HEAPDICT_H
template <typename K, typename V>
class HeapDict {
public:
    HeapDict();

    void clear();

    void insert(const K& key, const V& value);

    void erase(const K& key);

    V operator[](const K& key) const;

    size_t size() const;

    std::pair<K, V> popItem();

    std::pair<K, V> peekItem() const;

private:
    struct Wrapper {
        V value;
        K key;
        size_t position;

        Wrapper(const V& v, const K& k, size_t pos) : value(v), key(k), position(pos) {}
    };

    std::vector<std::shared_ptr<Wrapper>> heap;
    std::unordered_map<K, std::shared_ptr<Wrapper>> dictionary;

    size_t parent(size_t i) const;

    size_t left(size_t i) const;

    size_t right(size_t i) const;

    void decreaseKey(size_t i);

    void minHeapify(size_t i);

    void swap(size_t i, size_t j);
};

#endif  HEAPDICT_H
