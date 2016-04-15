#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Patch.hh"
#include "Circulator.hh"
using namespace kt84;

curvenetwork::Halfchain::Halfchain(int id_)
    : id        (id_)
    , is_deleted()
    , flag      ()
{}

curvenetwork::HalfchainPtr curvenetwork::Halfchain::next() const {
    if (halfedge_back ->next)
        return halfedge_back ->next->halfchain;
    return HalfchainPtr();
}

curvenetwork::HalfchainPtr curvenetwork::Halfchain::prev() const {
    if (halfedge_front->prev)
        return halfedge_front->prev->halfchain;
    return HalfchainPtr();
}

curvenetwork::HalfchainPtr curvenetwork::Halfchain::opposite() const {
    return halfedge_front->opposite->halfchain;
}

bool curvenetwork::Halfchain::is_corner() const {
    return halfedge_back->is_corner;
}

int curvenetwork::Halfchain::num_subdiv() const {
    return edgechain->num_subdiv;
}

curvenetwork::VertexPtr curvenetwork::Halfchain::vertex_front() const {
    return halfedge_front->from_vertex();
}

curvenetwork::VertexPtr curvenetwork::Halfchain::vertex_back () const {
    return halfedge_back ->vertex;
}

bool curvenetwork::Halfchain::is_loop() const {
    return halfedge_front == halfedge_back->next;
}

Polyline_PointNormal curvenetwork::Halfchain::toPolyline() const {
    Polyline_PointNormal result;
    result.push_back(halfedge_front->from_vertex()->pn);
    for (ConstCHIter h(this); h; ++h) {
        result.push_back(h->vertex->pn);
    }
    if (is_loop()) {
        result.pop_back();
        result.is_loop = true;
    }
    return result;
}

double curvenetwork::Halfchain::length() const {
    return toPolyline().length();
}

void curvenetwork::Halfchain::resample() const {
    auto polyline = toPolyline();
    polyline.resample();
    int i = 1;
    for (auto h = ConstCHIter(this); &*h != halfedge_back; ++i, ++h)
        h->vertex->pn = polyline[i];
}

PointNormal curvenetwork::Halfchain::pn_at(double t) const {
    return toPolyline().point_at(t);
}

OpenMesh::PolyConnectivity::VHandle curvenetwork::Halfchain::patch_vertex_at(int index) const {
    if (!Patch::is_ordinary(patch) || index < 0 || num_subdiv() < index)
        return Patch::VHandle();
    
    auto h_start = patch->prev_halfedge_handle(patch->halfedge_handle(patch->patch_corner(0)));
    // look for patch side corresponding to this halfchain
    for (auto h = h_start; ; ) {
        if (patch->data(h).halfchain == this) {
            // found! get the corresponding boundary vertex
            for (int i = 0; i < index; ++i)
                h = patch->prev_halfedge_handle(h);
            return patch->to_vertex_handle(h);
        }

        h = patch->prev_halfedge_handle(h);
        if (h == h_start)
            break;
    }
    // should not reach here...
    return Patch::VHandle();
}
bool curvenetwork::Halfchain::is_boundary() const {
    return !Patch::is_ordinary(patch) && Patch::is_ordinary(opposite()->patch);
}
