#pragma once
#include "ElemPtrT.hh"
#include <Eigen/Core>

namespace curvenetwork {
    struct Halfedge {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;
        
        // connectivity info
        VertexPtr    vertex   ;         // vertex to which this halfedge points
        HalfedgePtr  next     ;
        HalfedgePtr  prev     ;
        HalfedgePtr  opposite ;
        HalfchainPtr halfchain;
        
        // additional data
        bool is_corner;         // true if this and next halfedge forms a corner
        Patch* imaginary_patch;
        
        explicit Halfedge(int id_ = -1);
        
        // connectivity-related functions
        VertexPtr from_vertex() const;
        PatchPtr  patch      () const;
        
        // other utility functions
        Eigen::Vector3d toVector3d() const;
    };
}
