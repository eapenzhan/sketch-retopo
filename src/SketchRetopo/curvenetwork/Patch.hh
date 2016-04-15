#pragma once
#include "ElemPtrT.hh"
#include <patchgen/PatchParam.hh>
#include <patchgen/PatchVertexTraits.hh>
#include <utility>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <kt84/openmesh/base/LaplaceDirect.hh>
#include <kt84/openmesh/base/LaplaceIterative.hh>
#include <kt84/openmesh/base/Utility.hh>
#include <kt84/openmesh/base/DebugInfo.hh>
#include <kt84/graphics/DisplayList.hh>
#include <kt84/geometry/PolylineT.hh>

namespace curvenetwork {
    struct PatchTraits : public OpenMesh::DefaultTraits {
        typedef OpenMesh::Vec3d Point;              // these are basically unused
        typedef OpenMesh::Vec3d Normal;
        
        VertexTraits
            , public patchgen::PatchVertexTraits
            , public kt84::LaplaceDirect_VertexTraits   <6>
            , public kt84::LaplaceIterative_VertexTraits<6>
        {
            kt84::PointNormal     pn;
            kt84::PointNormal     pn_temp;      // temporary data used by smoothing
            Eigen::Vector2d harmonic_uv;        // Canonical UVs used for SketchRetopo::compute_patch_interior_pn_harmonic()
            
            VertexT()
                : pn     (kt84::PointNormal::Zero())
                , pn_temp(kt84::PointNormal::Zero())
            {}
        };
        
        HalfedgeTraits
            , public kt84::LaplaceDirect_HalfedgeTraits
            , public kt84::LaplaceIterative_HalfedgeTraits
        {
            Halfchain* halfchain;
            int        index_wrt_halfchain;
            
            HalfedgeT()
                : halfchain()
                , index_wrt_halfchain(-1)
            {}
        };
    };
    
    typedef OpenMesh::PolyMesh_ArrayKernelT<PatchTraits> PatchBase;
    
    struct Patch
        : public PatchBase
        , public kt84::LaplaceDirect   <PatchBase, Patch, 6>
        , public kt84::LaplaceIterative<PatchBase, Patch, 6>
        , public kt84::Utility         <PatchBase, Patch>
#ifndef NDEBUG
        , public kt84::DebugInfo       <PatchBase, Patch>
#endif
    {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;
        
        // connectivity info
        HalfchainPtr halfchain;
        
        // flag indicating the failure of patch generation because of impossible combination of num_subdiv of its surrounding edgechains. rendered in different color with no edge lines.
        bool is_failure;
        
        bool is_hidden;
        
        // display list
        kt84::DisplayList displist_phong_fill;
        kt84::DisplayList displist_phong_line;
        kt84::DisplayList displist_phong_line_boundary;
        kt84::DisplayList displist_flat_fill ;
        kt84::DisplayList displist_flat_line ;
        
        static Patch imaginary_patch_symmetry;              // imaginary patch to represent curves on symmetry plane      (id: -2)
        static Patch imaginary_patch_boundary;              // imaginary patch to represent curves on basemesh boundaries (id: -3)
        static double quadSize;                             // size of individual quads, used to determine Edgechain::num_subdiv. (don't know if it should really be declared here...)
        static bool   use_even_num_subdiv;
        static bool   prefer_rect_proc3;
        static const int failure_patch_num_subdiv;
        
        explicit Patch(int id_ = -1);
        void clear();
        
        // connectivity-related functions
        int num_corners() const;
        HHandle opposite_boundary_halfedge(HHandle boundary_halfedge) const;
        Patch*  opposite_patch(HHandle boundary_halfedge) const;
        curvenetwork::Vertex* vertex_patch_to_curveNetwork(VHandle v) const;
        VHandle vertex_curveNetwork_to_patch(curvenetwork::Vertex* v) const;
        
        // other utility functions
        VHandle patch_corner         (int corner_index) const;
        int     num_subdiv           (int side_index) const;
        int     num_free_halfchains  (int side_index) const;
        double  side_length          (int side_index) const;
        void    change_num_subdiv    (int side_index, int target_num_subdiv);
        std::pair<HHandle, HHandle> get_side_halfedges(int side_index) const;
        kt84::Polyline_PointNormal  get_side_polyline (int side_index) const;
        double distance(const Eigen::Vector3d& p) const;
        
        // patch topology generation
        patchgen::PatchParam param;
        void add_padding();
        void generate_topology(bool init_param = true);
        
        void set_halfedgeData ();
        void laplaceSolver_init();
        void set_boundary_pn();

        bool is_imaginary_symmetry() const;
        bool is_imaginary_boundary() const;
        bool is_imaginary() const;
        
        // rendering
        void render_phong_fill();
        void render_phong_line();
        void render_phong_line_boundary();
        void render_flat_fill();
        void render_flat_line();
        void render_irregular_vertices() const;
        void invalidate_displist();
        
        // static utility
        static bool is_ordinary(const Patch* patch);
    };
}
