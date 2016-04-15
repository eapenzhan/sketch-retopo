#pragma once
#include <Eigen/Core>
#include <vector>

// stupid but useful stuff for debugging with visual studio's immediate window

namespace dbg {

#define KT84DBG_H_MAKE_EIGEN_MATRIX_FUNC_SUB1(Type, Size)\
    Type&                           at  (Eigen::Matrix<Type, Size, Size>& m, int i, int j);\
    Eigen::Matrix<Type, Size, 1>    row (Eigen::Matrix<Type, Size, Size>& m, int i, Eigen::Matrix<Type, Size, 1>& dst);\
    Eigen::Matrix<Type, Size, 1>    col (Eigen::Matrix<Type, Size, Size>& m, int i, Eigen::Matrix<Type, Size, 1>& dst);\
    Type                            det (Eigen::Matrix<Type, Size, Size>& m);\
    Eigen::Matrix<Type, Size, Size> add (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> sub (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> mul (Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> mul (Type s, Eigen::Matrix<Type, Size, Size>& m, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> wsum(Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Type w1, Type w2, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> wsum(Eigen::Matrix<Type, Size, Size>& m1, Eigen::Matrix<Type, Size, Size>& m2, Eigen::Matrix<Type, Size, Size>& m3, Type w1, Type w2, Type w3, Eigen::Matrix<Type, Size, Size>& dst);\
    Eigen::Matrix<Type, Size, Size> rand(Eigen::Matrix<Type, Size, Size>& m, int rows, int cols);\
    Eigen::Matrix<Type, Size, Size> zero(Eigen::Matrix<Type, Size, Size>& m, int rows, int cols);\
    Eigen::Matrix<Type, Size, Size> rand(Eigen::Matrix<Type, Size, Size>& m);\
    Eigen::Matrix<Type, Size, Size> zero(Eigen::Matrix<Type, Size, Size>& m);\

#define KT84DBG_H_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, Rows, Cols)\
    Type&                           at   (Eigen::Matrix<Type, Rows, Cols>& v, int i);\
    Type                            dot  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2);\
    double                          norm(Eigen::Matrix<Type, Rows, Cols>& v);\
    Eigen::Matrix<Type, Rows, Cols> normalize(Eigen::Matrix<Type, Rows, Cols>& v);\
    Eigen::Matrix<Type, Rows, Cols> cross(Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst);\
    Eigen::Matrix<Type, Rows, Cols> add  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst);\
    Eigen::Matrix<Type, Rows, Cols> sub  (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& dst);\
    Eigen::Matrix<Type, Rows, Cols> mul  (Type s, Eigen::Matrix<Type, Rows, Cols>& v, Eigen::Matrix<Type, Rows, Cols>& dst);\
    /* TODO: matrix-vector multiplication */\
    Eigen::Matrix<Type, Rows, Cols> wsum (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Type w1, Type w2, Eigen::Matrix<Type, Rows, Cols>& dst);\
    Eigen::Matrix<Type, Rows, Cols> wsum (Eigen::Matrix<Type, Rows, Cols>& v1, Eigen::Matrix<Type, Rows, Cols>& v2, Eigen::Matrix<Type, Rows, Cols>& v3, Type w1, Type w2, Type w3, Eigen::Matrix<Type, Rows, Cols>& dst);\
    Eigen::Matrix<Type, Rows, Cols> rand(Eigen::Matrix<Type, Rows, Cols>& v, int size);\
    Eigen::Matrix<Type, Rows, Cols> zero(Eigen::Matrix<Type, Rows, Cols>& v, int size);\
    Eigen::Matrix<Type, Rows, Cols> rand(Eigen::Matrix<Type, Rows, Cols>& v);\
    Eigen::Matrix<Type, Rows, Cols> zero(Eigen::Matrix<Type, Rows, Cols>& v);\

#define KT84DBG_H_MAKE_EIGEN_VECTOR_FUNC_SUB2(Type, Size)\
    KT84DBG_H_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, Size, 1)\
    KT84DBG_H_MAKE_EIGEN_VECTOR_FUNC_SUB1(Type, 1, Size)\

#define KT84DBG_H_MAKE_EIGEN_FUNC_SUB1(Type, Size)\
    KT84DBG_H_MAKE_EIGEN_MATRIX_FUNC_SUB1(Type, Size)\
    KT84DBG_H_MAKE_EIGEN_VECTOR_FUNC_SUB2(Type, Size)\

#define KT84DBG_H_MAKE_EIGEN_FUNC_SUB2(Type)\
    KT84DBG_H_MAKE_EIGEN_FUNC_SUB1(Type,  2)\
    KT84DBG_H_MAKE_EIGEN_FUNC_SUB1(Type,  3)\
    KT84DBG_H_MAKE_EIGEN_FUNC_SUB1(Type,  4)\
    KT84DBG_H_MAKE_EIGEN_FUNC_SUB1(Type, -1)\

#define KT84DBG_H_MAKE_STL_FUNC_SUB1(Type)\
    Type& at(std::vector<Type>& c, int i);\

const int pool_size = 100;

#define KT84DBG_H_MAKE_VARIABLE_SUB1(Type, TypeSuffix, Size, SizeSuffix)\
    extern Eigen::Matrix<Type, Size, Size> mat##SizeSuffix##TypeSuffix [pool_size];\
    extern Eigen::Matrix<Type, Size, 1>    vec##SizeSuffix##TypeSuffix [pool_size];\
    extern Eigen::Matrix<Type, 1, Size>    cvec##SizeSuffix##TypeSuffix[pool_size];\

#define KT84DBG_H_MAKE_VARIABLE_SUB2(Type, TypeSuffix)\
    extern Type x##TypeSuffix[pool_size];\
    KT84DBG_H_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  2, 2)\
    KT84DBG_H_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  3, 3)\
    KT84DBG_H_MAKE_VARIABLE_SUB1(Type, TypeSuffix,  4, 4)\
    KT84DBG_H_MAKE_VARIABLE_SUB1(Type, TypeSuffix, -1, x)\

#define KT84DBG_H_MAKE_SUB1(Type, TypeSuffix)\
    KT84DBG_H_MAKE_EIGEN_FUNC_SUB2(Type)\
    KT84DBG_H_MAKE_STL_FUNC_SUB1(Type)\
    KT84DBG_H_MAKE_VARIABLE_SUB2(Type, TypeSuffix)\

KT84DBG_H_MAKE_SUB1(int   , i)
KT84DBG_H_MAKE_SUB1(double, d)
KT84DBG_H_MAKE_SUB1(float , f)

#define KT84DBG_H_EIGEN_COPY_SUB1(TypeSrc, TypeDst, Size)\
    void copy(Eigen::Matrix<TypeSrc, Size, Size>& src, Eigen::Matrix<TypeDst, Size, Size>& dst);\
    void copy(Eigen::Matrix<TypeSrc, Size, 1   >& src, Eigen::Matrix<TypeDst, Size, 1   >& dst);\
    void copy(Eigen::Matrix<TypeSrc, 1   , Size>& src, Eigen::Matrix<TypeDst, 1   , Size>& dst);\

#define KT84DBG_H_EIGEN_COPY_SUB2(TypeSrc, TypeDst)\
    KT84DBG_H_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  2)\
    KT84DBG_H_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  3)\
    KT84DBG_H_EIGEN_COPY_SUB1(TypeSrc, TypeDst,  4)\
    KT84DBG_H_EIGEN_COPY_SUB1(TypeSrc, TypeDst, -1)\

KT84DBG_H_EIGEN_COPY_SUB2(int   , int   )
KT84DBG_H_EIGEN_COPY_SUB2(int   , double)
KT84DBG_H_EIGEN_COPY_SUB2(int   , float )
KT84DBG_H_EIGEN_COPY_SUB2(double, int   )
KT84DBG_H_EIGEN_COPY_SUB2(double, double)
KT84DBG_H_EIGEN_COPY_SUB2(double, float )
KT84DBG_H_EIGEN_COPY_SUB2(float , int   )
KT84DBG_H_EIGEN_COPY_SUB2(float , double)
KT84DBG_H_EIGEN_COPY_SUB2(float , float )

}

/*

if you're compiling dbg.cpp with your project instead of linking against *.lib,
make sure that these functions are included in the binary by specifying the linker option:
    Linker->Optimization->References->No (/OPT:NOREF)

example code:

int main() {
    MatrixXd V = MatrixXd::Random(100, 3);
    auto vec = util::make_vector<int>(5, 1, 2, 3, 4, 5);
    dbg::at(V, 0, 0);           // <-- dummy call to any function in dbg::, otherwise the Immediate Window won't find symbols, weird...
    return 0;
}


we can work in the immediate window like this:

dbg::at(V, 0, 2)
-0.12662129581591230
dbg::at(vec, 3)
4
dbg::row(V, 0, dbg::vecxd[0])
[3] (dynamic row vector)
    [size]: 3
    [0]: -0.99749748222296819
    [1]: -0.70500808740501109
    [2]: -0.12662129581591230
    [Raw View]: 0x005df930 {...}
dbg::row(V, 1, dbg::vecxd[1])
[3] (dynamic row vector)
    [size]: 3
    [0]: 0.12717062898648024
    [1]: 0.89916684469130526
    [2]: -0.53514206366161077
    [Raw View]: 0x005df930 {...}
dbg::row(V, 2, dbg::vecxd[2])
[3] (dynamic row vector)
    [size]: 3
    [0]: -0.61339152195806756
    [1]: -0.71684926908169810
    [2]: 0.15573595385601369
    [Raw View]: 0x005df930 {...}
dbg::sub(dbg::vecxd[1], dbg::vecxd[0], dbg::vecxd[3])
[3] (dynamic row vector)
    [size]: 3
    [0]: 1.1246681112094485
    [1]: 1.6041749320963163
    [2]: -0.40852076784569846
    [Raw View]: 0x005df930 {...}
dbg::sub(dbg::vecxd[2], dbg::vecxd[0], dbg::vecxd[4])
[3] (dynamic row vector)
    [size]: 3
    [0]: 0.38410596026490063
    [1]: -0.011841181676687018
    [2]: 0.28235724967192599
    [Raw View]: 0x005df930 {...}
dbg::cross(dbg::vecxd[3], dbg::vecxd[4], dbg::vecxd[5])
[3] (dynamic row vector)
    [size]: 3
    [0]: 0.44811305318860389
    [1]: -0.47447345649634626
    [2]: -0.62949055215654492
    [Raw View]: 0x005df930 {...}
dbg::normalize(dbg::vecxd[5])
[3] (dynamic row vector)
    [size]: 3
    [0]: 0.49419897060183438
    [1]: -0.52327039373186424
    [2]: -0.69423012935177841
    [Raw View]: 0x005df930 {...}

*/
