#include "dbg.hh"
#include <cmath>
#include <Eigen/LU>

namespace dbg {

#define KT84DBG_CPP_MAKE_EIGEN_MATRIX_FUNC_SUB1(Type, Size)\
    Type&                           at  (Eigen::Matrix<Type, Size, Size>& m, int i, int j) { return m(i, j); }\
    Eigen::Matrix<Type, Size, 1>    row (Eigen::Matrix<Type, Size, Size>& m, int i, Eigen::Matrix<Type, Size, 1>& dst) { return dst = m.row(i); }\
    Eigen::Matrix<Type, Size, 1>    col (Eigen::Matrix<Type, Size, Size>& m, int i, Eigen::Matrix<Type, Size, 1>& dst) { return dst = m.col(i); }\
    Type                            det (Eigen::Matrix<Type, Size, Size>& m) { return m.determinant(); }\
    Eigen::Matrix<Type, Size, Size> add (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst) { return dst = m1 + m2; }\
    Eigen::Matrix<Type, Size, Size> sub (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst) { return dst = m1 - m2; }\
    Eigen::Matrix<Type, Size, Size> mul (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst) { return dst = m1 * m2; }\
    Eigen::Matrix<Type, Size, Size> mul (Type s, Eigen::Matrix<Type, Size, Size>& m, Eigen::Matrix<Type, Size, Size>& dst) { return dst = s * m; }\
    Eigen::Matrix<Type, Size, Size> wsum(Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Type w1, Type w2, Eigen::Matrix<Type, Size, Size>& dst) { return dst = w1 * m1 + w2 * m2; }\
    Eigen::Matrix<Type, Size, Size> wsum(Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& m3, Type w1, Type w2, Type w3, Eigen::Matrix<Type, Size, Size>& dst) { return dst = w1 * m1 + w2 * m2 + w3 * m3; }\
    Eigen::Matrix<Type, Size, Size> rand(Eigen::Matrix<Type, Size, Size>& m, int rows, int cols) { return m = Eigen::Matrix<Type, Size, Size>::Random(rows, cols); }\
    Eigen::Matrix<Type, Size, Size> zero(Eigen::Matrix<Type, Size, Size>& m, int rows, int cols) { return m = Eigen::Matrix<Type, Size, Size>::Zero  (rows, cols); }\
    Eigen::Matrix<Type, Size, Size> rand(Eigen::Matrix<Type, Size, Size>& m) { return m = Eigen::Matrix<Type, Size, Size>::Random(Size, Size); }\
    Eigen::Matrix<Type, Size, Size> zero(Eigen::Matrix<Type, Size, Size>& m) { return m = Eigen::Matrix<Type, Size, Size>::Zero  (Size, Size); }\

#define KT84DBG_CPP_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, Rows, Cols)\
    Type&                           at   (Eigen::Matrix<Type, Rows, Cols>& v, int i) { return v[i]; }\
    Type                            dot  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2) { return v1.dot(v2); }\
    double                          norm(Eigen::Matrix<Type, Rows, Cols>& v) { return std::sqrt(v.dot(v)); }\
    Eigen::Matrix<Type, Rows, Cols> normalize(Eigen::Matrix<Type, Rows, Cols>& v) { return v = v / norm(v); }\
    Eigen::Matrix<Type, Rows, Cols> cross(Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst) {\
        dst.resize(3);\
        dst[0] = v1[1] * v2[2] - v1[2] * v2[1];\
        dst[1] = v1[2] * v2[0] - v1[0] * v2[2];\
        dst[2] = v1[0] * v2[1] - v1[1] * v2[0];\
        return dst;\
    }\
    Eigen::Matrix<Type, Rows, Cols> add  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst) { return dst = v1 + v2; }\
    Eigen::Matrix<Type, Rows, Cols> sub  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst) { return dst = v1 - v2; }\
    Eigen::Matrix<Type, Rows, Cols> mul  (Type s, Eigen::Matrix<Type, Rows, Cols>& v, Eigen::Matrix<Type, Rows, Cols>& dst) { return dst = s * v; }\
    /* TODO: matrix-vector multiplication */\
    Eigen::Matrix<Type, Rows, Cols> wsum (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Type w1, Type w2, Eigen::Matrix<Type, Rows, Cols>& dst) { return dst = w1 * v1 + w2 * v2; }\
    Eigen::Matrix<Type, Rows, Cols> wsum (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& v3, Type w1, Type w2, Type w3, Eigen::Matrix<Type, Rows, Cols>& dst) { return dst = w1 * v1 + w2 * v2 + w3 * v3; }\
    Eigen::Matrix<Type, Rows, Cols> rand(Eigen::Matrix<Type, Rows, Cols>& v, int size) { return v = Eigen::Matrix<Type, Rows, Cols>::Random(size); }\
    Eigen::Matrix<Type, Rows, Cols> zero(Eigen::Matrix<Type, Rows, Cols>& v, int size) { return v = Eigen::Matrix<Type, Rows, Cols>::Zero  (size); }\
    Eigen::Matrix<Type, Rows, Cols> rand(Eigen::Matrix<Type, Rows, Cols>& v) { return v = Eigen::Matrix<Type, Rows, Cols>::Random(Rows * Cols); }\
    Eigen::Matrix<Type, Rows, Cols> zero(Eigen::Matrix<Type, Rows, Cols>& v) { return v = Eigen::Matrix<Type, Rows, Cols>::Zero  (Rows * Cols); }\

#define KT84DBG_CPP_MAKE_EIGEN_VECTOR_FUNC_SUB2(Type, Size)\
    KT84DBG_CPP_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, Size, 1)\
    KT84DBG_CPP_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, 1, Size)\

#define KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB1(Type, Size)\
    KT84DBG_CPP_MAKE_EIGEN_MATRIX_FUNC_SUB1(Type, Size)\
    KT84DBG_CPP_MAKE_EIGEN_VECTOR_FUNC_SUB2(Type, Size)\

#define KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB2(Type)\
    KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB1(Type,  2)\
    KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB1(Type,  3)\
    KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB1(Type,  4)\
    KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB1(Type, -1)\

#define KT84DBG_CPP_MAKE_STL_FUNC_SUB1(Type)\
    Type& at(std::vector<Type>& c, int i) { return c[i]; }\

#define KT84DBG_CPP_MAKE_VARIABLE_SUB1(Type, TypeSuffix, Size, SizeSuffix)\
    Eigen::Matrix<Type, Size, Size> mat##SizeSuffix##TypeSuffix [pool_size];\
    Eigen::Matrix<Type, Size, 1>    vec##SizeSuffix##TypeSuffix [pool_size];\
    Eigen::Matrix<Type, 1, Size>    cvec##SizeSuffix##TypeSuffix[pool_size];\

#define KT84DBG_CPP_MAKE_VARIABLE_SUB2(Type, TypeSuffix)\
    Type x##TypeSuffix[pool_size];\
    KT84DBG_CPP_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  2, 2)\
    KT84DBG_CPP_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  3, 3)\
    KT84DBG_CPP_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  4, 4)\
    KT84DBG_CPP_MAKE_VARIABLE_SUB1(Type, TypeSuffix, -1, x)\

#define KT84DBG_CPP_MAKE_SUB1(Type, TypeSuffix)\
    KT84DBG_CPP_MAKE_EIGEN_FUNC_SUB2(Type)\
    KT84DBG_CPP_MAKE_STL_FUNC_SUB1(Type)\
    KT84DBG_CPP_MAKE_VARIABLE_SUB2(Type, TypeSuffix)\

KT84DBG_CPP_MAKE_SUB1(int   , i)
KT84DBG_CPP_MAKE_SUB1(double, d)
KT84DBG_CPP_MAKE_SUB1(float , f)

#define KT84DBG_CPP_EIGEN_COPY_SUB1(TypeSrc, TypeDst, Size)\
    void copy(Eigen::Matrix<TypeSrc, Size, Size>& src, Eigen::Matrix<TypeDst, Size, Size>& dst) { dst = src.cast<TypeDst>(); }\
    void copy(Eigen::Matrix<TypeSrc, Size, 1   >& src, Eigen::Matrix<TypeDst, Size, 1   >& dst) { dst = src.cast<TypeDst>(); }\
    void copy(Eigen::Matrix<TypeSrc, 1   , Size>& src, Eigen::Matrix<TypeDst, 1   , Size>& dst) { dst = src.cast<TypeDst>(); }\

#define KT84DBG_CPP_EIGEN_COPY_SUB2(TypeSrc, TypeDst)\
    KT84DBG_CPP_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  2)\
    KT84DBG_CPP_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  3)\
    KT84DBG_CPP_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  4)\
    KT84DBG_CPP_EIGEN_COPY_SUB1(TypeSrc, TypeDst, -1)\

KT84DBG_CPP_EIGEN_COPY_SUB2(int   , int   )
KT84DBG_CPP_EIGEN_COPY_SUB2(int   , double)
KT84DBG_CPP_EIGEN_COPY_SUB2(int   , float )
KT84DBG_CPP_EIGEN_COPY_SUB2(double, int   )
KT84DBG_CPP_EIGEN_COPY_SUB2(double, double)
KT84DBG_CPP_EIGEN_COPY_SUB2(double, float )
KT84DBG_CPP_EIGEN_COPY_SUB2(float , int   )
KT84DBG_CPP_EIGEN_COPY_SUB2(float , double)
KT84DBG_CPP_EIGEN_COPY_SUB2(float , float )

}
