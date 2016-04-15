#pragma once
#include <Eigen/Geometry>
#include "../../vector_cast.hh"
#include "../DerivedPtrHolder.hh"

namespace kt84 {

template <class TMeshBase, class TMesh>
struct BoundingBox : public DerivedPtrHolder<TMesh, BoundingBox<TMeshBase, TMesh>> {
    static const int N = TMeshBase::Point::size_;
    Eigen::AlignedBox<double, N> boundingBox;

    void boundingBox_compute() {
        TMesh* mesh = DerivedPtrHolder<TMesh, BoundingBox<TMeshBase, TMesh>>::derived_ptr;
        boundingBox.setEmpty();
        for (auto v : mesh->vertices())
            boundingBox.extend(vector_cast<N, Eigen::Matrix<double, N, 1>>(mesh->point(v)));
    }
    double boundingBox_diagonal_norm() const {
        return boundingBox.diagonal().norm();
    }
};

}
