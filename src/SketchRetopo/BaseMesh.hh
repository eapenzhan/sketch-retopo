#pragma once

#include <Eigen/Core>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <kt84/openmesh/base/ExpMap.hh>
#include <kt84/openmesh/base/CotanWeight.hh>
#include <kt84/openmesh/base/CenterOfMass.hh>
#include <kt84/openmesh/base/BoundingBox.hh>
#include <kt84/openmesh/base/LaplaceIterative.hh>
#include <kt84/openmesh/base/NormalVariation.hh>
#include <kt84/openmesh/base/NormalCurvature.hh>
#include <kt84/openmesh/base/Gradient.hh>
#include <kt84/openmesh/base/Utility.hh>
#ifndef NDEBUG
#   include <kt84/openmesh/base/DebugInfo.hh>
#endif
#include <kt84/geometry/PointNormal.hh>
#include <kt84/graphics/DisplayList.hh>

struct BaseMeshTraits : public OpenMesh::DefaultTraits {
    typedef OpenMesh::Vec3d Point;
    typedef OpenMesh::Vec3d Normal;
    typedef OpenMesh::Vec2d TexCoord2D;
    
    VertexTraits
        , public kt84::ExpMap_VertexTraits
        , public kt84::LaplaceIterative_VertexTraits<3>
        , public kt84::NormalCurvature_VertexTraits
        , public kt84::NormalVariation_VertexTraits
        , public kt84::Gradient_VertexTraits
    {
        Eigen::Vector3d normal_smoothed;
        Eigen::Vector2d k1k2;
        double          snake_feature;
        double          snake_energy;
        Eigen::Vector3d snake_energy_gradient;
        Eigen::Vector2d harmonic_uv;                // UV computed with harmonic parameterization in SketchRetopo::compute_patch_interior_pn_harmonic(). Not to be confused with the texture UVs
        
        VertexT()
            : normal_smoothed(Eigen::Vector3d::Zero())
            , k1k2(Eigen::Vector2d::Zero())
            , snake_feature()
            , snake_energy()
            , snake_energy_gradient(Eigen::Vector3d::Zero())
        {}
    };
    
    FaceTraits
        , public kt84::ExpMap_FaceTraits
        , public kt84::CotanWeight_FaceTraits
        , public kt84::Gradient_FaceTraits<3>
    {
        kt84::BaryCoordZero snake_feature_gradient;
        kt84::BaryCoordZero snake_energy_gradient;
        bool floodfill_flag;
        bool is_hidden;
        FaceT()
            : floodfill_flag()
            , is_hidden()
        {}
    };
    
    EdgeTraits
        , public kt84::CotanWeight_EdgeTraits
    {};
    
    HalfedgeTraits
        , public kt84::LaplaceIterative_HalfedgeTraits
        , public kt84::CotanWeight_HalfedgeTraits
    {};
};

typedef OpenMesh::TriMesh_ArrayKernelT<BaseMeshTraits> BaseMeshBase;

struct BaseMesh
    : public BaseMeshBase
    , public kt84::ExpMap          <BaseMeshBase, BaseMesh>
    , public kt84::CotanWeight     <BaseMeshBase, BaseMesh>
    , public kt84::CenterOfMass    <BaseMeshBase, BaseMesh>
    , public kt84::BoundingBox     <BaseMeshBase, BaseMesh>
    , public kt84::LaplaceIterative<BaseMeshBase, BaseMesh, 3>
    , public kt84::NormalCurvature <BaseMeshBase, BaseMesh>
    , public kt84::NormalVariation <BaseMeshBase, BaseMesh>
    , public kt84::Gradient        <BaseMeshBase, BaseMesh>
    , public kt84::Utility         <BaseMeshBase, BaseMesh>
#ifndef NDEBUG
    , public kt84::DebugInfo       <BaseMeshBase, BaseMesh>
#endif
{
    kt84::DisplayList displist_fill;
    kt84::DisplayList displist_line;
    int normalSmoothIter;
    int featureSmoothIter;
    double snake_energy_max;
    double snake_energy_min;
    bool use_feature_as_energy;
    bool negate_energy;
    bool has_texture;
    
    enum class SnakeFeatureType {
        MEAN_CURVATURE = 0,
        NORMAL_VARIATION,
        NORMAL_CURVATURE
    } snake_feature_type;
    
    BaseMesh();
    void init();
    void smooth_normal();
    void smooth_expmap();
    void compute_k1k2();
    void compute_snake_energy();
    void render_fill();
    void render_line();
    void invalidate_displist() {
        displist_fill.invalidate();
        displist_line.invalidate();
    }
    kt84::PointNormal get_pn(VHandle v) const;
};
