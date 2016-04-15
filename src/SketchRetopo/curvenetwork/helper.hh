#pragma once
#include "decl.hh"

namespace curvenetwork {
    // helper functions
    void set_next_prev(Halfedge* h_pred, Halfedge* h_succ);
    void set_vertex_halfedge(Vertex* v0, Vertex*  v1, Halfedge* h01, Halfedge* h10);
    void trace_halfchains(Halfchain* in_c, Halfchain*& out_c_front, Halfchain*& out_c_back);
    void delete_patch    (Patch    * p);
    void delete_halfchain(Halfchain* c);
    void delete_halfedge (Halfedge * h);
}
