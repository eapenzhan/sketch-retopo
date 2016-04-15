#include "helper.hh"
#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Patch.hh"
#include "Circulator.hh"

void curvenetwork::set_next_prev(Halfedge* h_pred, Halfedge* h_succ) {
    if (h_pred) h_pred->next = h_succ;
    if (h_succ) h_succ->prev = h_pred;
}
void curvenetwork::set_vertex_halfedge(Vertex* v0, Vertex* v1, Halfedge* h01, Halfedge* h10) {
    v0->halfedge = h01;
    v1->halfedge = h10;
    h01->vertex = v1;
    h10->vertex = v0;
    h01->opposite = h10;
    h10->opposite = h01;
}
void curvenetwork::trace_halfchains(Halfchain* in_c, Halfchain*& out_c_front, Halfchain*& out_c_back) {
    auto c_forward  = in_c;
    auto c_backward = in_c;
    for (int i = 0; ; ++i) {
        if (i % 2 == 0 && c_forward) {
            if (!c_forward->next())
                out_c_back = c_forward;
            c_forward = c_forward->next();
        } else if (i % 2 == 1 && c_backward) {
            if (!c_backward->prev())
                out_c_front = c_backward;
            c_backward = c_backward->prev();
        }
        if (c_forward == c_backward)
            break;
    }
    if (c_forward) {
        // loop detected
        out_c_front = in_c;
        out_c_back  = in_c->prev();
    }
}
void curvenetwork::delete_patch(Patch* p) {
    if (!p || p->is_deleted || p->is_imaginary()) return;     // never delete imaginary patch
    p->is_deleted = true;
    
    // clear pointer from halfchains to this face
    for (PCIter c(p); c; ++c) {
        c->patch = nullptr;
        
        if (c->opposite() && !Patch::is_ordinary(c->opposite()->patch))
            c->edgechain->num_subdiv = 0;
    }
}
void curvenetwork::delete_halfchain(Halfchain* c) {
    if (!c || c->is_deleted) return;
    c->is_deleted = true;
    
    delete_patch(c->patch);
    
    // clear pointer from halfedges to this halfchain
    for (CHIter h(c); h; ++h)
        h->halfchain = 0;
    
    // delete its edgechain if isolated
    auto e = c->edgechain;
    if (e->halfchain[0]->is_deleted && e->halfchain[1]->is_deleted)
        e->is_deleted = true;
}
void curvenetwork::delete_halfedge(Halfedge* h) {
    if (!h || h->is_deleted) return;
    h->is_deleted = true;
    
    delete_halfchain(h->halfchain);
    
    auto v = h->from_vertex();
    // check if v is isolated
    bool is_isolated = true;
    for (VOHIter voh(v); voh; ++voh) {
        if (voh->is_deleted) continue;
        v->halfedge = &*voh;
        is_isolated = false;
        break;
    }
    if (is_isolated)
        // v is isolated --> delete
        v->is_deleted = true;
}
