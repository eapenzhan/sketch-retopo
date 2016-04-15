#pragma once
#include "ElemPtrT.hh"
#include <kt84/geometry/PolylineT.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>

namespace curvenetwork {
    struct Halfchain {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;
        
        // connectivity info
        HalfedgePtr  halfedge_front;
        HalfedgePtr  halfedge_back;
        PatchPtr     patch;
        EdgechainPtr edgechain;
        
        explicit Halfchain(int id_ = -1);
        
        // connectivity-related functions
        HalfchainPtr next        () const;
        HalfchainPtr prev        () const;
        HalfchainPtr opposite    () const;
        bool         is_corner   () const;
        int          num_subdiv  () const;
        VertexPtr    vertex_front() const;
        VertexPtr    vertex_back () const;
        
        // other utility functions
        bool is_loop() const;
        kt84::Polyline_PointNormal toPolyline() const;
        double length() const;
        void resample() const;
        kt84::PointNormal pn_at(double t) const;
        OpenMesh::PolyConnectivity::VHandle patch_vertex_at(int index) const;       // index should be within [0, num_subdiv()]
        bool is_boundary() const;
    };
}
