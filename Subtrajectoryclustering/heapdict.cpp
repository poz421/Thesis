#include "heapdict.h"

template <typename K, typename V>
HeapDict<K,V>::HeapDict() {}

    template <typename K, typename V>
    void HeapDict<K, V>::clear() {
        heap.clear();
        dictionary.clear();
    }

    template <typename K, typename V>
    void HeapDict<K, V>::insert(const K& key, const V& value) {
        if (dictionary.find(key) != dictionary.end()) {
            erase(key);
        }

        auto wrapper = std::make_shared<Wrapper>(value, key, heap.size());
        dictionary[key] = wrapper;
        heap.push_back(wrapper);
        decreaseKey(heap.size() - 1);
    }

    template <typename K, typename V>
    void HeapDict<K, V>::erase(const K& key) {
        auto wrapper = dictionary[key];
        while (wrapper->position) {
            auto parentPos = parent(wrapper->position);
            auto parent = heap[parentPos];
            swap(wrapper->position, parent->position);
        }
        popItem();
    }

    template <typename K, typename V>
    V HeapDict<K, V>::operator[](const K& key) const {
        return dictionary.at(key)->value;
    }

    template <typename K, typename V>
    size_t HeapDict<K, V>::size() const {
        return dictionary.size();
    }

    template <typename K, typename V>
    std::pair<K, V> HeapDict<K, V>::popItem() {
        if (heap.size() == 1) {
            auto wrapper = heap[0];
            heap.pop_back();
            dictionary.erase(wrapper->key);
            return std::make_pair(wrapper->key, wrapper->value);
        }
        else {
            auto wrapper = heap[0];
            heap[0] = heap.back();
            heap[0]->position = 0;
            heap.pop_back();
            minHeapify(0);
            dictionary.erase(wrapper->key);
            return std::make_pair(wrapper->key, wrapper->value);
        }
    }

    template <typename K, typename V>
    std::pair<K, V> HeapDict<K, V>::peekItem() const {
        return std::make_pair(heap[0]->key, heap[0]->value);
    }
   
    template <typename K, typename V>
    size_t HeapDict<K, V>::parent(size_t i) const {
        return (i - 1) >> 1;
    }

    template <typename K, typename V>
    size_t HeapDict<K, V>::left(size_t i) const {
        return (i << 1) + 1;
    }

    template <typename K, typename V>
    size_t HeapDict<K, V>::right(size_t i) const {
        return (i + 1) << 1;
    }

    template <typename K, typename V>
    void HeapDict<K, V>::decreaseKey(size_t i) {
        while (i > 0) {
            size_t parentIndex = parent(i);
            if (heap[parentIndex]->value < heap[i]->value) {
                break;
            }
            swap(i, parentIndex);
            i = parentIndex;
        }
    }

    template <typename K, typename V>
    void HeapDict<K, V>::minHeapify(size_t i) {
        size_t l = left(i);
        size_t r = right(i);
        size_t n = heap.size();

        size_t lowest = 0;
        if (l < n and heap[l]->value < heap[i]->value)
            lowest = l;
        else
            lowest = i;
        if (r < n and heap[r]->value < heap[lowest]->value)
            lowest = r;

        if (lowest != i) {
            swap(i, lowest);
            minHeapify(lowest);
        }
    }

    template <typename K, typename V>
    void HeapDict<K, V>::swap(size_t i, size_t j) {
        std::swap(heap[i], heap[j]);
        heap[i]->position = i;
        heap[j]->position = j;
    }

