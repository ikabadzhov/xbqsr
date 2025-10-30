#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <cstdint>

class StringDictionary {
public:
    static constexpr int32_t EMPTY_STRING     = -1;
    static constexpr int32_t STRING_NOT_FOUND = -2;

    int32_t add(const std::string& s) {
        if (s.empty()) return EMPTY_STRING;
        auto it = s2id_.find(s);
        if (it != s2id_.end()) return it->second;
        auto id = static_cast<int32_t>(id2s_.size());
        s2id_.emplace(s, id);
        id2s_.push_back(s);
        return id;
    }
    int32_t getId(const std::string& s) const {
        if (s.empty()) return EMPTY_STRING;
        auto it = s2id_.find(s);
        return it == s2id_.end() ? STRING_NOT_FOUND : it->second;
    }
    const std::string& getString(int32_t id) const {
        if (id < 0 || static_cast<size_t>(id) >= id2s_.size())
            throw std::out_of_range("StringDictionary: id out of range");
        return id2s_[static_cast<size_t>(id)];
    }
    int32_t size() const noexcept { return static_cast<int32_t>(id2s_.size()); }

    void clear() {
        s2id_.clear();
        id2s_.clear();
    }

private:
    std::unordered_map<std::string,int32_t> s2id_;
    std::vector<std::string>                id2s_;
};