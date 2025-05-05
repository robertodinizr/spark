
#ifndef ENUM_BIT_OPS_H
#define ENUM_BIT_OPS_H

#define ENUM_CLASS_BIT_OPS(enum_class_name, under_type)                                \
    inline constexpr enum_class_name operator&(enum_class_name x, enum_class_name y) { \
        return static_cast<enum_class_name>(static_cast<under_type>(x) &               \
                                            static_cast<under_type>(y));               \
    }                                                                                  \
                                                                                       \
    inline constexpr enum_class_name operator|(enum_class_name x, enum_class_name y) { \
        return static_cast<enum_class_name>(static_cast<under_type>(x) |               \
                                            static_cast<under_type>(y));               \
    }                                                                                  \
                                                                                       \
    inline constexpr enum_class_name operator^(enum_class_name x, enum_class_name y) { \
        return static_cast<enum_class_name>(static_cast<under_type>(x) ^               \
                                            static_cast<under_type>(y));               \
    }                                                                                  \
                                                                                       \
    inline constexpr enum_class_name operator~(enum_class_name x) {                    \
        return static_cast<enum_class_name>(~static_cast<under_type>(x));              \
    }                                                                                  \
                                                                                       \
    inline enum_class_name& operator&=(enum_class_name& x, enum_class_name y) {        \
        x = x & y;                                                                     \
        return x;                                                                      \
    }                                                                                  \
                                                                                       \
    inline enum_class_name& operator|=(enum_class_name& x, enum_class_name y) {        \
        x = x | y;                                                                     \
        return x;                                                                      \
    }                                                                                  \
                                                                                       \
    inline enum_class_name& operator^=(enum_class_name& x, enum_class_name y) {        \
        x = x ^ y;                                                                     \
        return x;                                                                      \
    }

#endif  // ENUM_BIT_OPS_H
