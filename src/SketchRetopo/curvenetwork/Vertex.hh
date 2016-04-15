#pragma once
#include "ElemPtrT.hh"
#include <kt84/geometry/PointNormal.hh>

namespace curvenetwork {
    struct Vertex {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;               // flag for general purpose
        
        // connectivity info
        HalfedgePtr halfedge;     // outgoing halfedge
        
        // additional data
        kt84::PointNormal pn;
        bool is_hidden;
        
        explicit Vertex(int id_ = -1);
        bool is_openend () const;
        bool is_corner  () const;
        bool is_boundary() const;
        int  valence    () const;
        bool on_symmetry_plane() const;
    };
}
