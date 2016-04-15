#pragma once
#include <list>
#include <kt84/graphics/DisplayList.hh>
#include <kt84/geometry/PolylineT.hh>
#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Patch.hh"

namespace curvenetwork {
    struct Core {
        // elements stored using std::list
        std::list<Vertex   > vertices;
        std::list<Halfedge > halfedges;
        std::list<Halfchain> halfchains;
        std::list<Edgechain> edgechains;
        std::list<Patch    > patches;
        
        // display list for elements with relatively complex rendering scheme
        kt84::DisplayList displist_corner;
        
        void add_curve_loop(kt84::Polyline_PointNormal& curve);                                                     // non-const argument is intentional here...
        void add_curve_open(kt84::Polyline_PointNormal& curve, double dist_snap, bool corner_at_openend = true);
        void erase_curve(const kt84::Polyline_PointNormal& scribble);
        void erase_curve_sub(Halfedge* h);
        
        void flip_corner_valance2(Vertex* v);
        void flip_corner_valanceN(Vertex* v, Eigen::Vector3d reference_point);
    
        void generate_halfchains();
    
        enum SnapType {
            SNAPTYPE_OPENEND = 0,
            SNAPTYPE_CORNER,
            SNAPTYPE_SIDE,
        };
        bool calc_snap(SnapType snapType, Vertex*& v_snap, Vertex* v_other, Halfedge*& h_snap, double dist_snap, bool corner_at_openend);       // h_snap: halfedfge pointing toward v_snap (computed if snapping actually happens)
    
        bool calc_intersect(Vertex* v0, Vertex* v1, Halfedge*& h_v0_u, Halfedge*& h_u_v1);                // h_v0_u, h_u_v1: computed if snapping actually happens
    
        enum class FlagOp {
            SUBST  = 0,
            AND   ,
            OR    ,
        };
        void set_flag_vertices  (FlagOp op, int arg);
        void set_flag_halfedges (FlagOp op, int arg);
        void set_flag_halfchains(FlagOp op, int arg);
        void set_flag_edgechains(FlagOp op, int arg);
        void set_flag_patches   (FlagOp op, int arg);
        
        VertexPtr    new_vertex   ();
        HalfedgePtr  new_halfedge ();
        HalfchainPtr new_halfchain();
        EdgechainPtr new_edgechain();
        PatchPtr     new_patch    ();
    
        void validate_all_ptr();
    
        void garbage_collect();
    
        Core& operator=(const Core& rhs);
    
        void clear();
        bool empty() const;
    
        // rendering
        void render_edgechains() const;
        void render_edgechains_num_subdiv() const;
        void render_corners();
        void render_endpoints() const;
        void invalidate_displist();
    };
}
