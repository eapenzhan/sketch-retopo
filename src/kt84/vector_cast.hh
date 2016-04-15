#pragma once
namespace kt84 {
    namespace internal {
        // to achieve partial specialization on function template
        template <int N, typename T_to, typename T_from>
        struct VectorCast {
            static T_to cast(const T_from& v) {
                T_to result;
                for (int i = 0; i < N; ++i)
                    result[i] = v[i];
                return result;
            }
        };
        template <typename T_to, typename T_from>
        struct VectorCast<2, T_to, T_from> {
            static T_to cast(const T_from& v) { return { v[0], v[1] }; }
        };
        template <typename T_to, typename T_from>
        struct VectorCast<3, T_to, T_from> {
            static T_to cast(const T_from& v) { return { v[0], v[1], v[2] }; }
        };
        template <typename T_to, typename T_from>
        struct VectorCast<4, T_to, T_from> {
            static T_to cast(const T_from& v) { return { v[0], v[1], v[2], v[3] }; }
        };
    }
    template <int N, typename T_to, typename T_from>
    inline T_to vector_cast(const T_from& v) {
        return internal::VectorCast<N, T_to, T_from>::cast(v);
    }
}
