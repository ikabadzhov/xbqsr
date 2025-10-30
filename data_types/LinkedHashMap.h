#ifndef GENOMICPIPELINE_LINKEDHASHMAP_H
#define GENOMICPIPELINE_LINKEDHASHMAP_H

#include <vector>
#include <unordered_map>
#include <stdexcept>

template <typename K, typename V>
class LinkedHashMap {
private:
    std::vector<K> keyOrder;
    std::unordered_map<K, V> values;
    std::unordered_map<K, bool> is_removed;

public:
    // Inserts a key-value pair into the map
    void put(const K& key, const V& value) {
        if (!containsKey(key)) {
            keyOrder.push_back(key);
        }
        values[key] = value;
        is_removed[key] = false;
    }

    const V& get(const K& key) const {
        if (!containsKey(key)) {
            throw std::out_of_range("Key not found");
        }
        return values.at(key);
    }

    // Removes the key-value pair for a given key
    void remove(const K& key) {
        if (containsKey(key)) {
            keyOrder.erase(std::remove(keyOrder.begin(), keyOrder.end(), key), keyOrder.end());
            values.erase(key);
            is_removed.erase(key);
        }
    }

    // Mark as removed the key-value pair for a given key
    void soft_remove(const K& key) {
        if (containsKey(key)) {
            is_removed[key] = true;
        }
    }

    // Checks if the map contains a given key
    bool containsKey(const K& key) const {
        if ((is_removed.contains(key) && is_removed.at(key)) ||  values.find(key) == values.end()) {
            return false;
        }
        return true;
    }

    // Iterates over and applies a function to each key-value pair
    template <typename Func>
    void for_each(Func f) const {
        for (const auto& key : keyOrder) {
            if (!is_removed.at(key)) {
                f(key, values.at(key));
            }
        }
    }

    // Returns the size of the map
    [[nodiscard]] size_t size() const {
        return values.size();
    }

    // Checks if the map is empty
    [[nodiscard]] bool isEmpty() const {
        return values.empty();
    }

    // More functions can be added to further mimic Java's LinkedHashMap
};

#endif //GENOMICPIPELINE_LINKEDHASHMAP_H
