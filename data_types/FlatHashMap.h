#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>
#include <new>
#include <cassert>

// ========================= Fast string hash =========================

struct fast_str_hash {
    std::size_t operator()(const std::string& s) const noexcept {
        return hash_bytes(reinterpret_cast<const unsigned char*>(s.data()), s.size());
    }
    std::size_t operator()(std::string_view s) const noexcept {
        return hash_bytes(reinterpret_cast<const unsigned char*>(s.data()), s.size());
    }

private:
    static std::size_t hash_bytes(const unsigned char* p, std::size_t n) noexcept {
        // FNV-1a 64-bit (simple, stable, good enough)
        std::uint64_t h = 1469598103934665603ull;
        for (std::size_t i = 0; i < n; ++i) {
            h ^= p[i];
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

// ========================= flat_hash_map (simplified) =========================

namespace fast {

// Very small SwissTable-style map:
// - open addressing, linear probing
// - ctrl[i] == 0 => empty
// - ctrl[i] == tag (1..255) => occupied, tag = high 8 bits of hash (never 0)
// - no erase, no rehash-on-delete, no fancy batching

template<class Key, class T, class Hash = std::hash<Key>>
class flat_hash_map {
    struct Entry {
        Key key;
        T   val;
    };

    static constexpr std::uint8_t EMPTY = 0;
    static constexpr float MAX_LOAD = 0.5f;

    std::uint8_t* ctrl_   = nullptr;
    Entry*        slots_  = nullptr;
    std::size_t   cap_    = 0;   // always power of two when > 0
    std::size_t   mask_   = 0;   // cap_ - 1
    std::size_t   size_   = 0;
    Hash          hasher_{};

public:
    flat_hash_map() = default;

    explicit flat_hash_map(std::size_t n, const Hash& h = Hash{})
        : hasher_(h)
    {
        if (n) {
            reserve(n);
        }
    }

    ~flat_hash_map() {
        destroy_storage();
    }

    flat_hash_map(const flat_hash_map&) = delete;
    flat_hash_map& operator=(const flat_hash_map&) = delete;

    flat_hash_map(flat_hash_map&& other) noexcept {
        steal(other);
    }

    flat_hash_map& operator=(flat_hash_map&& other) noexcept {
        if (this != &other) {
            destroy_storage();
            steal(other);
        }
        return *this;
    }

    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }

    void reserve(std::size_t n) {
        const std::size_t want = static_cast<std::size_t>(n / MAX_LOAD) + 1;
        if (want <= cap_) return;
        rehash(next_pow2(want));
    }

    // Insert-or-assign.
    template<class K, class... Args>
    T& emplace_or_assign(K&& key, Args&&... args) {
        ensure_capacity_for_one_more();

        std::size_t h = hasher_(key);
        std::uint8_t tag = make_tag(h);
        std::size_t idx = h & mask_;

        for (;;) {
            std::uint8_t c = ctrl_[idx];

            if (c == EMPTY) {
                // insert new
                ctrl_[idx] = tag;
                ::new (&slots_[idx].key) Key(std::forward<K>(key));
                ::new (&slots_[idx].val) T(std::forward<Args>(args)...);
                ++size_;
                return slots_[idx].val;
            }

            if (c == tag && slots_[idx].key == key) {
                // overwrite
                slots_[idx].val = T(std::forward<Args>(args)...);
                return slots_[idx].val;
            }

            idx = (idx + 1) & mask_;
        }
    }

    // Convenience insert for your benchmark: emplace(Key,T)
    std::pair<T*, bool> insert(const Key& key, const T& val) {
        ensure_capacity_for_one_more();

        std::size_t h = hasher_(key);
        std::uint8_t tag = make_tag(h);
        std::size_t idx = h & mask_;

        for (;;) {
            std::uint8_t c = ctrl_[idx];

            if (c == EMPTY) {
                ctrl_[idx] = tag;
                ::new (&slots_[idx].key) Key(key);
                ::new (&slots_[idx].val) T(val);
                ++size_;
                return { &slots_[idx].val, true };
            }

            if (c == tag && slots_[idx].key == key) {
                return { &slots_[idx].val, false }; // already there
            }

            idx = (idx + 1) & mask_;
        }
    }

    // Find: returns pointer to value or nullptr
    T* find(const Key& key) {
        if (!cap_) return nullptr;

        std::size_t h = hasher_(key);
        std::uint8_t tag = make_tag(h);
        std::size_t idx = h & mask_;

        for (;;) {
            std::uint8_t c = ctrl_[idx];

            if (c == EMPTY) {
                return nullptr;
            }

            if (c == tag && slots_[idx].key == key) {
                return &slots_[idx].val;
            }

            idx = (idx + 1) & mask_;
        }
    }

    const T* find(const Key& key) const {
        return const_cast<flat_hash_map*>(this)->find(key);
    }

private:
    static std::size_t next_pow2(std::size_t x) {
        if (x < 8) return 8;
        x--;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
#if INTPTR_MAX == INT64_MAX
        x |= x >> 32;
#endif
        return x + 1;
    }

    static std::uint8_t make_tag(std::size_t h) {
        std::uint8_t t = static_cast<std::uint8_t>(h >> (sizeof(std::size_t)*8 - 8));
        return t ? t : 1; // never 0 = reserved for empty
    }

    void ensure_capacity_for_one_more() {
        if (cap_ == 0) {
            rehash(8);
        } else if (size_ + 1 > static_cast<std::size_t>(cap_ * MAX_LOAD)) {
            rehash(cap_ * 2);
        }
    }

    void rehash(std::size_t new_cap) {
        assert((new_cap & (new_cap - 1)) == 0 && "capacity must be power of two");

        std::uint8_t* new_ctrl = new std::uint8_t[new_cap];
        Entry*        new_slots = static_cast<Entry*>(::operator new[](new_cap * sizeof(Entry)));

        for (std::size_t i = 0; i < new_cap; ++i)
            new_ctrl[i] = EMPTY;

        const std::size_t old_cap = cap_;
        std::uint8_t* old_ctrl = ctrl_;
        Entry* old_slots = slots_;

        cap_  = new_cap;
        mask_ = new_cap - 1;
        ctrl_ = new_ctrl;
        slots_= new_slots;
        size_ = 0;

        if (old_slots) {
            for (std::size_t i = 0; i < old_cap; ++i) {
                if (old_ctrl[i] != EMPTY) {
                    // move k,v into new table
                    Key& k = old_slots[i].key;
                    T&   v = old_slots[i].val;
                    emplace_or_assign(std::move(k), std::move(v));
                    k.~Key();
                    v.~T();
                }
            }
            delete[] old_ctrl;
            ::operator delete[](old_slots);
        }
    }

    void destroy_storage() {
        if (!slots_) return;
        for (std::size_t i = 0; i < cap_; ++i) {
            if (ctrl_[i] != EMPTY) {
                slots_[i].key.~Key();
                slots_[i].val.~T();
            }
        }
        delete[] ctrl_;
        ::operator delete[](slots_);
        ctrl_ = nullptr;
        slots_ = nullptr;
        cap_ = mask_ = size_ = 0;
    }

    void steal(flat_hash_map& o) noexcept {
        ctrl_ = o.ctrl_;
        slots_ = o.slots_;
        cap_ = o.cap_;
        mask_ = o.mask_;
        size_ = o.size_;
        hasher_ = std::move(o.hasher_);

        o.ctrl_ = nullptr;
        o.slots_ = nullptr;
        o.cap_ = o.mask_ = o.size_ = 0;
    }
};

} // namespace fast

//using fast_string_map = fast::flat_hash_map<std::string,int,fast_str_hash>;
