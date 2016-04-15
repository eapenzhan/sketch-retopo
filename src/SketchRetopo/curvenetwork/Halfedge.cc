#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"

curvenetwork::Halfedge::Halfedge(int id_)
    : id        (id_)
    , is_deleted()
    , flag      ()
    , is_corner(false)
    , imaginary_patch(nullptr)
{}
curvenetwork::VertexPtr curvenetwork::Halfedge::from_vertex() const {
    return opposite->vertex;
}
curvenetwork::PatchPtr  curvenetwork::Halfedge::patch() const {
    return halfchain->patch;
}
Eigen::Vector3d curvenetwork::Halfedge::toVector3d() const {
    return (vertex->pn - from_vertex()->pn).head(3);
}
