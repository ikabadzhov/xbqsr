struct transparent_string_hash {
    using is_transparent = void;  // enables hetero lookup

    std::size_t operator()(std::string_view s) const noexcept {
        return std::hash<std::string_view>{}(s);
    }
    std::size_t operator()(const std::string& s) const noexcept {
        return std::hash<std::string>{}(s);
    }
    std::size_t operator()(const char* s) const noexcept {
        return std::hash<std::string_view>{}(std::string_view{s});
    }
};

struct transparent_string_equal {
    using is_transparent = void;  // enables hetero lookup

    bool operator()(std::string_view a,
                    std::string_view b) const noexcept {
        return a == b;
    }
};