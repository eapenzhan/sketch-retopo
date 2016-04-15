#include "AutopoMesh.hh"

bool AutopoMesh::is_corner(VHandle v) const {
    int num_separatrices = 0;
    for (auto e : ve_range(v))
        num_separatrices += data(e).on_separatrix;
    return num_separatrices == valence(v);
}
bool AutopoMesh::is_regular(VHandle v) const {
    int valence_ = valence(v);
    return is_boundary(v) && valence_ == 3 || valence_ == 4;
}
bool AutopoMesh::on_separatrix(VHandle v) const {
    for (auto e : ve_range(v))
        if (data(e).on_separatrix)
            return true;
    return false;
}
AutopoMesh::HHandle AutopoMesh::edgeflow_next(HHandle h) const {
    if (is_boundary(h))
        return next_halfedge_handle(h);
    
    return next_halfedge_handle(opposite_halfedge_handle(next_halfedge_handle(h)));
}
AutopoMesh::HHandle AutopoMesh::edgeflow_prev(HHandle h) const {
    if (is_boundary(h))
        return prev_halfedge_handle(h);
    
    return prev_halfedge_handle(opposite_halfedge_handle(prev_halfedge_handle(h)));
}
void AutopoMesh::trace_separatrix() {
    for (auto v : vertices()) {
        if (is_regular(v))
            continue;
        
        for (auto h_start : voh_range(v)) {
            if (data(edge_handle(h_start)).on_separatrix)
                continue;
            
            for (auto h = h_start; ; ) {
                data(edge_handle(h)).on_separatrix = true;
                
                auto v2 = to_vertex_handle(h);
                if (!is_regular(v2) || !is_boundary(h) && is_boundary(v2))
                    break;
                
                h = edgeflow_next(h);
            }
        }
    }
    
    for (auto e : edges()) {
        if (is_boundary(e))
            data(e).on_separatrix = true;
    }
}
