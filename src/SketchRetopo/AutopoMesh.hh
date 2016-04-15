#pragma once
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <kt84/openmesh/base/Utility.hh>
#include <kt84/openmesh/base/DebugInfo.hh>
#include "curvenetwork/decl.hh"

struct AutopoMeshTraits : public OpenMesh::DefaultTraits {
    typedef OpenMesh::Vec3d Point;
    typedef OpenMesh::Vec3d Normal;
    
    VertexTraits {
        curvenetwork::Vertex* cn_vertex;
        VertexT()
            : cn_vertex(nullptr)
        {}
    };
    
    EdgeTraits {
        bool on_separatrix;
        EdgeT()
            : on_separatrix(false)
        {}
    };
    
    HalfedgeTraits {
        curvenetwork::Halfedge* cn_halfedge;
        bool is_processed;
        HalfedgeT()
            : cn_halfedge(nullptr)
            , is_processed(false)
        {}
    };
};

typedef OpenMesh::PolyMesh_ArrayKernelT<AutopoMeshTraits> AutopoMeshBase;

struct AutopoMesh
    : public AutopoMeshBase
    , public kt84::Utility  <AutopoMeshBase, AutopoMesh>
#ifndef NDEBUG
    , public kt84::DebugInfo<AutopoMeshBase, AutopoMesh>
#endif
{
    
    AutopoMesh() {
        request_face_normals();
        request_vertex_normals();
    }
    void init() {
        update_normals();
#ifndef NDEBUG
        debugInfo_get(true, true, true);
#endif
    }
    bool is_corner   (VHandle v) const;
    bool is_regular(VHandle v) const;
    bool on_separatrix(VHandle v) const;
    HHandle edgeflow_next(HHandle h) const;
    HHandle edgeflow_prev(HHandle h) const;
    void trace_separatrix();
};
