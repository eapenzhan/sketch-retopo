#include "Core.hh"
#include "Circulator.hh"
#include "helper.hh"
#include <Eigen/LU>
#include <kt84/eigen_util.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

void curvenetwork::Core::add_curve_open(Polyline_PointNormal& curve, double dist_snap, bool corner_at_openend) {
    // clear flag
    set_flag_vertices (FlagOp::SUBST, 0);
    set_flag_halfedges(FlagOp::SUBST, 0);
    
    // add vertex
    vector<Vertex*> vertices_added;
    for (int i = 0; i < curve.size(); ++i) {
        auto v = new_vertex();
        v->flag |= 2;           // flag for newly generated elements
        v->pn = curve[i];
        vertices_added.push_back(v);
    }
    
    // add halfedge without next/prev info
    vector<Halfedge*> halfedges_forward;
    vector<Halfedge*> halfedges_backward;
    for (int i = 0; i < curve.size() - 1; ++i) {
        Vertex* v0 = vertices_added[i];
        Vertex* v1 = vertices_added[i + 1];
        Halfedge* h_snap = 0;
        Halfedge* h_v0_u = 0;
        Halfedge* h_u_v1 = 0;
        if (i == 0 &&
            (  calc_snap(SNAPTYPE_OPENEND, v0, v1, h_snap, dist_snap, corner_at_openend)
            || calc_snap(SNAPTYPE_CORNER , v0, v1, h_snap, dist_snap, corner_at_openend)
            || calc_snap(SNAPTYPE_SIDE   , v0, v1, h_snap, dist_snap, corner_at_openend)))
        {
            // snapping happend at curve starting point
            halfedges_forward .push_back(h_snap->opposite);
            halfedges_backward.push_back(h_snap);
        } else if (i == curve.size() - 2 &&
            (  calc_snap(SNAPTYPE_OPENEND, v1, v0, h_snap, dist_snap, corner_at_openend)
            || calc_snap(SNAPTYPE_CORNER , v1, v0, h_snap, dist_snap, corner_at_openend)
            || calc_snap(SNAPTYPE_SIDE   , v1, v0, h_snap, dist_snap, corner_at_openend)))
        {
            // snapping happend at curve ending point
            halfedges_forward .push_back(h_snap);
            halfedges_backward.push_back(h_snap->opposite);
        } else if (calc_intersect(v0, v1, h_v0_u, h_u_v1)) {
            // intersection detected
            halfedges_forward .push_back(h_v0_u);
            halfedges_forward .push_back(h_u_v1);
            halfedges_backward.push_back(h_v0_u->opposite);
            halfedges_backward.push_back(h_u_v1->opposite);
        } else {
            auto h01 = new_halfedge();
            auto h10 = new_halfedge();
            h01->flag = h10->flag = 2;        // flag for newly generated elements
            set_vertex_halfedge(v0, v1, h01, h10);
            halfedges_forward .push_back(h01);
            halfedges_backward.push_back(h10);
        }
    }
    
    // set halfedge prev/next info
    size_t n = halfedges_forward.size();
    for (size_t i = 0; i < n; ++i) {
        auto h0 = halfedges_forward[i];
        if (!h0->prev && i != 0    ) h0->prev = halfedges_forward[i - 1];
        if (!h0->next && i != n - 1) h0->next = halfedges_forward[i + 1];
        auto h1 = halfedges_backward[i];
        if (!h1->prev && i != n - 1) h1->prev = halfedges_backward[i + 1];
        if (!h1->next && i != 0    ) h1->next = halfedges_backward[i - 1];
    }
    
    generate_halfchains();
    
    garbage_collect();
}

void curvenetwork::Core::add_curve_loop(Polyline_PointNormal& curve) {
    assert(curve.is_loop);
    
    // clear flag
    set_flag_vertices (FlagOp::SUBST, 0);
    set_flag_halfedges(FlagOp::SUBST, 0);
    
    // add vertex
    vector<Vertex*> vertices_added;
    for (int i = 0; i < curve.size(); ++i) {
        auto v = new_vertex();
        v->flag = 2;                            // flag for newly generated elements
        v->pn = curve[i];
        vertices_added.push_back(v);
    }
    
    // add halfedge without next/prev info
    vector<Halfedge*> halfedges_forward;
    vector<Halfedge*> halfedges_backward;
    for (int i = 0; i < curve.size(); ++i) {
        Vertex* v0 = vertices_added[i];
        Vertex* v1 = vertices_added[i < curve.size() - 1 ? i + 1 : 0];
        Halfedge* h_v0_u = 0;
        Halfedge* h_u_v1 = 0;
        if (calc_intersect(v0, v1, h_v0_u, h_u_v1)) {
            halfedges_forward .push_back(h_v0_u);
            halfedges_forward .push_back(h_u_v1);
            halfedges_backward.push_back(h_v0_u->opposite);
            halfedges_backward.push_back(h_u_v1->opposite);
        } else {
            auto h01 = new_halfedge();
            auto h10 = new_halfedge();
            h01->flag = h10->flag = 2;        // flag for newly generated elements
            set_vertex_halfedge(v0, v1, h01, h10);
            halfedges_forward .push_back(h01);
            halfedges_backward.push_back(h10);
        }
    }
    
    // set halfedge prev/next info
    for (size_t i = 0; i < halfedges_forward.size(); ++i) {
        int i_next = i < halfedges_forward.size() - 1 ? i + 1 : 0;
        int i_prev = i == 0 ? halfedges_forward.size() - 1 : i - 1;
        auto h0 = halfedges_forward[i];
        if (!h0->prev) h0->prev = halfedges_forward[i_prev];
        if (!h0->next) h0->next = halfedges_forward[i_next];
        auto h1 = halfedges_backward[i];
        if (!h1->prev) h1->prev = halfedges_backward[i_next];
        if (!h1->next) h1->next = halfedges_backward[i_prev];
    }
    
    generate_halfchains();
    
    garbage_collect();
}

void curvenetwork::Core::erase_curve(const Polyline_PointNormal& scribble) {
    for (int i = 1; i < scribble.size(); ++i) {
        Vector3d v0 = scribble[i - 1].head(3);
        Vector3d v1 = scribble[i    ].head(3);
        
        for (auto h = halfedges.begin(); h != halfedges.end(); ++h) {
            if (h->is_deleted)
                continue;
            
            if (h->imaginary_patch || h->opposite && h->opposite->imaginary_patch)
                // never erase halfedges adjacent to imaginary patch
                continue;
            
            // calc intersection
            Vector3d w0 = h->opposite->vertex->pn.head(3);
            Vector3d w1 = h          ->vertex->pn.head(3);
            
            // check for intersection of bbox
            double rv = (v1 - v0).cwiseAbs().maxCoeff() * 0.5;
            double rw = (w1 - w0).cwiseAbs().maxCoeff() * 0.5;
            double rvw = (v0 + v1 - (w0 + w1)).cwiseAbs().minCoeff() * 0.5;
            if (rv + rw < rvw)
                continue;
            
            // (1 - s) * v0 + s * v1 = (1 - t) * w0 + t * w1;
            // (v1 - v0) * s - (w1 - w0) * t = w0 - v0
            // p = v1 - v0
            // q = w1 - w0
            // b = w0 - v0
            // A = [p, -q]
            // x = inv(A^T * A) * A^T * b
            Vector3d p = v1 - v0;
            Vector3d q = w1 - w0;
            Vector3d b = w0 - v0;
            
            if (p.squaredNorm() * q.squaredNorm() == p.dot(q) * p.dot(q))
                // these edges are colinear
                continue;
            
            Matrix<double, 3, 2> A;
            A << p, -q;
            Vector2d x = (A.transpose() * A).inverse() * A.transpose() * b;
            if (x[0] < 0 || x[1] < 0 || 1 < x[0] || 1 < x[1])
                // intersected outside either edges
                continue;
            
            double dist = (A * x - b).norm();
            if ((p.norm() + q.norm()) < dist)
                // intersections are far away
                continue;
            
            // intersection detected
            if (!h->halfchain)
                // this should not happen!
                break;
            
            erase_curve_sub(&*h);
            break;
        }
    }
}
void curvenetwork::Core::erase_curve_sub(Halfedge* h) {
    if (h->is_deleted)
        return;
    
    auto h_front = h->halfchain->halfedge_front;
    auto h_back  = h->halfchain->halfedge_back;
    
    // update next/prev info around endpoints of halfchain
    pair<Halfedge*, Halfedge*> prev_next_pair[2];
    int valence[2];
    Halfedge* h_endpoint[2] = { h_front->opposite, h_back };
    for (int i = 0; i < 2; ++i) {
        valence[i] = h_endpoint[i]->vertex->valence();
        prev_next_pair[i] = make_pair(h_endpoint[i]->opposite->prev, h_endpoint[i]->next);
    }
    
    // delete all incident halfedge/halfchain/face
    for (auto h1 = h_front; ; h1 = h1->next) {
        delete_halfedge(h1);
        delete_halfedge(h1->opposite);
        if (h1 == h_back) break;
    }
    
    for (int i = 0; i < 2; ++i) {
        auto h_prev = prev_next_pair[i].first;
        auto h_next = prev_next_pair[i].second;
        if (valence[i] == 2) {          // make open endpoint
            h_prev->next = nullptr;
            h_next->prev = nullptr;
            h_prev->is_corner = false;
        
        } else if (valence[i] > 2) {
            set_next_prev(h_prev, h_next);
            int num_corners = 0;
            for (VIHIter vih(h_prev->vertex); vih; ++vih)
                num_corners += vih->is_corner ? 1 : 0;
            h_prev->is_corner = valence[i] == 3 || num_corners < 2;
        }
    }
    
    garbage_collect();
    
}

void curvenetwork::Core::flip_corner_valance2(curvenetwork::Vertex* v) {
    if (v->valence() != 2) return;
    //assert(v->valence() == 2);
    
    // flip is_corner
    bool is_corner = v->halfedge->opposite->is_corner;
    v->halfedge->opposite->is_corner = !is_corner;
    v->halfedge->prev    ->is_corner = !is_corner;
    
    // delete associated halfchains
    delete_halfchain(v->halfedge->halfchain);
    delete_halfchain(v->halfedge->opposite->halfchain);
    delete_halfchain(v->halfedge->opposite->next->halfchain);
    delete_halfchain(v->halfedge->prev->halfchain);
    
    generate_halfchains();
    
    garbage_collect();
}
void curvenetwork::Core::flip_corner_valanceN(Vertex* v, Vector3d reference_point) {
    if (v->valence() < 3) return;
    
    auto d_ref = (reference_point - v->pn.head(3)).normalized();
    
    Halfedge* h_max = 0;
    double dot_max = -numeric_limits<double>::max();
    for (VIHIter vih(v); vih; ++vih) {
        auto d0 = -vih       ->toVector3d().normalized();
        auto d1 =  vih->next->toVector3d().normalized();
        double dot = d0.dot(d_ref) + d1.dot(d_ref);
        if (dot_max < dot) {
            dot_max = dot;
            h_max = &*vih;
        }
    }
    
    // flip is_corner
    h_max->is_corner = !h_max->is_corner;
    
    // delete associated halfchains
    delete_halfchain(h_max      ->halfchain);
    delete_halfchain(h_max->next->halfchain);
    
    generate_halfchains();
    
    garbage_collect();
}

void curvenetwork::Core::generate_halfchains() {
    for (auto h = halfedges.begin(); h != halfedges.end(); ++h) {
        if (h->is_deleted) continue;
        
        if (h->halfchain && !h->halfchain->is_deleted)
            continue;
        
        auto c = new_halfchain();
        
        // trace halfedges forward
        c->halfedge_back = &*h;
        while (true) {
            c->halfedge_back->halfchain = c;
            auto v = c->halfedge_back->vertex;
            if (v->is_corner() || v->is_openend())
                break;
            c->halfedge_back = c->halfedge_back->next;
            if (c->halfedge_back == &*h)
                // loop with no corner
                break;
        }
        if (c->halfedge_back == &*h && !h->vertex->is_corner() && !h->vertex->is_openend()) {
            // this halfchain is loop
            c->halfedge_front = &*h;
            c->halfedge_back = h->prev;
        } else {
            // trace halfedges backward
            c->halfedge_front = &*h;
            while (true) {
                 c->halfedge_front->halfchain = c;
                 auto v = c->halfedge_front->opposite->vertex;
                 if (v->is_corner() || v->is_openend())
                     break;
                 c->halfedge_front = c->halfedge_front->prev;
            }
        }
        
        // reference to corresponding edgechain
        if (c->opposite()) {
            // its opposite side is already created --> point to an existing edgechain
            auto e = c->opposite()->edgechain;
            e->halfchain[e->halfchain[0] == c->opposite() ? 1 : 0] = c;
            c->edgechain = e;
        } else {
            // its opposite side is not yet created --> create a new edgechain, point to it
            auto e = new_edgechain();
            e->halfchain[0] = c;
            c->edgechain = e;
            
            // resample evenly
            //c->resample();
        }
    }
}

bool curvenetwork::Core::calc_snap(SnapType snapType, Vertex*& v_snap, Vertex* v_other, Halfedge*& h_snap, double dist_snap, bool corner_at_openend) {
    double dist_min = DBL_MAX;
    Vertex* v_min = nullptr;
    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
        if (v->is_deleted) continue;
        
        if (v->flag & 2)
            // flag at second lowest bit is on --> newly generated vertex
            continue;
        
        bool is_openend = v->is_openend();
        bool is_corner  = v->is_corner ();
        
        if (snapType == SNAPTYPE_OPENEND)   { if (!is_openend)              continue; } else
        if (snapType == SNAPTYPE_CORNER )   { if (is_openend || !is_corner) continue; } else
        /* (snapType == SNAPTYPE_SIDE   )*/ { if (is_openend ||  is_corner) continue; }
        
        double dist = pn_norm(v_snap->pn - v->pn);
        if (dist < dist_min) {
            dist_min = dist;
            v_min = &*v;
        }
    }
    
    if (dist_snap < dist_min)
        // no snapping
        return false;
        
    Halfedge* h_snap_next = 0;              // in the end should satisfy h_snap     ->next == h_snap_next
    Halfedge* h_snap_prev = 0;              // in the end should satisfy h_snap_prev->next == h_snap_opposite
    
    if (snapType == SNAPTYPE_CORNER || snapType == SNAPTYPE_SIDE) {
        // get star configuration around v_min
        Vector3d n = v_min->pn.tail(3);
        Vector3d d0;                // pivot orientation
        vector<double> star_angles;
        for (VOHIter h(v_min); h; ++h) {
            Vector3d d = (h->vertex->pn - v_min->pn).head(3);
            d -= n.dot(d) * n;      // project along normal
            d.normalize();
            if (&*h == v_min->halfedge) {
                d0 = d;
                continue;
            }
            double angle = atan2(d0.cross(d).dot(n), d0.dot(d));
            if (angle < 0)
                angle += 2 * M_PI;
            star_angles.push_back(angle);
        }
        
        // find halfedge that should be connected by h_snap
        Vector3d d1 = (v_other->pn - v_min->pn).head(3);
        d1 -= n.dot(d1) * n;
        d1.normalize();
        double angle = atan2(d0.cross(d1).dot(n), d0.dot(d1));
        if (angle < 0)
            angle += 2 * M_PI;
        
        size_t i = 0;
        for (VOHIter h(v_min); i < star_angles.size(); ++h, ++i) {
            if (star_angles[i] < angle) {
                h_snap_next = &*++h;
                break;
            }
        }
        if (!h_snap_next)
            h_snap_next = v_min->halfedge;
        
        h_snap_prev = h_snap_next->prev;
        
    } else {    // snapType == SNAPTYPE_OPENEND
        h_snap_next = v_min->halfedge;
        h_snap_prev = h_snap_next->opposite;
    }
    
    if (snapType == SNAPTYPE_CORNER) {
        if (h_snap_next->halfchain) delete_patch(h_snap_next->halfchain->patch);        // why this checking needed...?
        if (h_snap_prev->halfchain) delete_patch(h_snap_prev->halfchain->patch);
        
    } else if (snapType == SNAPTYPE_SIDE) {
        delete_halfchain(h_snap_next->halfchain);
        delete_halfchain(h_snap_next->opposite->halfchain);
    }
    
    // create new halfedge
    h_snap = new_halfedge();        auto h_snap_opposite = new_halfedge();
    
    // flag for newly generated elements
    h_snap->flag |= 2;              h_snap_opposite->flag |= 2;
    
    // set connectivity info
    set_vertex_halfedge(v_other, v_min, h_snap, h_snap_opposite);
    
    // set next/prev info
    set_next_prev(h_snap     , h_snap_next    );
    set_next_prev(h_snap_prev, h_snap_opposite);
    
    // set corner flag
    bool is_corner = snapType != SNAPTYPE_OPENEND || corner_at_openend;
    h_snap->is_corner = h_snap_prev->is_corner = is_corner;
    if (!is_corner) {
        delete_halfchain(h_snap_next          ->halfchain);
        delete_halfchain(h_snap_next->opposite->halfchain);
    }
    
    // delete unneeded v_snap
    v_snap->is_deleted = true;

    return true;
}

bool curvenetwork::Core::calc_intersect(Vertex* v0, Vertex* v1, Halfedge*& h_v0_u, Halfedge*& h_u_v1) {
    // clear visited flag (while keeping other flags)
    set_flag_halfedges(FlagOp::AND, ~1);
    
    for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
        auto h = &*i;
        
        if (h->is_deleted) continue;
        
        // check if this edge is already visited or newly generated
        if (h->flag & 1 || h->flag & 2)
            continue;
        // set visited flag
        h          ->flag |= 1;
        h->opposite->flag |= 1;
        
        // calc intersection
        Vertex* w0 = h->opposite->vertex;
        Vertex* w1 = h->vertex;
        
        // check for intersection of bbox
        double rv = (v1->pn - v0->pn).head(3).cwiseAbs().maxCoeff() * 0.5;
        double rw = (w1->pn - w0->pn).head(3).cwiseAbs().maxCoeff() * 0.5;
        double rvw = (v0->pn + v1->pn - (w0->pn + w1->pn)).head(3).cwiseAbs().minCoeff() * 0.5;
        if (rv + rw < rvw)
            continue;
        
        // (1 - s) * v0 + s * v1 = (1 - t) * w0 + t * w1;
        // (v1 - v0) * s - (w1 - w0) * t = w0 - v0
        // p = v1 - v0
        // q = w1 - w0
        // b = w0 - v0
        // A = [p, -q]
        // x = inv(A^T * A) * A^T * b
        Vector3d p = (v1->pn - v0->pn).head(3);
        Vector3d q = (w1->pn - w0->pn).head(3);
        Vector3d b = (w0->pn - v0->pn).head(3);
        
        if (p.squaredNorm() * q.squaredNorm() == p.dot(q) * p.dot(q))
            // these edges are colinear
            continue;
        
        Matrix<double, 3, 2> A;
        A << p, -q;
        Vector2d x = (A.transpose() * A).inverse() * A.transpose() * b;
        if (x[0] < 0 || x[1] < 0 || 1 < x[0] || 1 < x[1])
            // intersected outside either edges
            continue;
        
        double dist = (A * x - b).norm();
        if ((p.norm() + q.norm()) < dist)
            // intersections are far away
            continue;
        
        // intersection detected
        auto u = new_vertex();
        
        // delete halfchains
        delete_halfchain(h          ->halfchain);
        delete_halfchain(h->opposite->halfchain);
        
        // flag for newly generated elements
        u->flag = 2;
        
        u->pn = 0.5 * ((1 - x[0]) * v0->pn + x[0] * v1->pn + (1 - x[1]) * w0->pn + x[1] * w1->pn);
        pn_normalize(u->pn);
        // normal at intersection
        Vector3d n = u->pn.tail(3);
        
        // make sure (v0-->v1, w0-->w1, n) forms a left-hand axis configuration
        if (p.cross(q).dot(n) < 0) {
            swap(w0, w1);
            h = h->opposite;
        }
        
        // generate 8 halfedges
        //        v1
        //        ||
        //        ||
        //  w1====u====w0
        //        ||
        //        ||
        //        v0
             h_v0_u = new_halfedge();   auto h_u_v0 = new_halfedge();
        auto h_v1_u = new_halfedge();        h_u_v1 = new_halfedge();
        auto h_w0_u = new_halfedge();   auto h_u_w0 = new_halfedge();
        auto h_w1_u = new_halfedge();   auto h_u_w1 = new_halfedge();
        
        // flag for newly generated elements
        h_v0_u->flag = 2;       h_u_v0->flag = 2;
        h_v1_u->flag = 2;       h_u_v1->flag = 2;
        h_w0_u->flag = 2;       h_u_w0->flag = 2;
        h_w1_u->flag = 2;       h_u_w1->flag = 2;
        
        // set connectivity info
        set_vertex_halfedge(u, v0, h_u_v0, h_v0_u);
        set_vertex_halfedge(u, v1, h_u_v1, h_v1_u);
        set_vertex_halfedge(u, w0, h_u_w0, h_w0_u);
        set_vertex_halfedge(u, w1, h_u_w1, h_w1_u);
        
        // set next/prev around u
        set_next_prev(h_v0_u, h_u_w1);
        set_next_prev(h_w1_u, h_u_v1);
        set_next_prev(h_v1_u, h_u_w0);
        set_next_prev(h_w0_u, h_u_v0);
        
        // set next/prev connecting to w0, w1 from outside
        if (h->next) {
            set_next_prev(h_u_w1, h->next);
            set_next_prev(h->opposite->prev, h_w1_u);
        }
        if (h->prev) {
            set_next_prev(h->prev, h_w0_u);
            set_next_prev(h_u_w0, h->opposite->next);
        }
        
        // set corner flag
        h_v0_u->is_corner = true;
        h_v1_u->is_corner = true;
        h_w0_u->is_corner = true;
        h_w1_u->is_corner = true;
        
        // delete this halfedge (and its opposite)
        delete_halfedge(h);
        delete_halfedge(h->opposite);
        
        return true;
    }
    return false;
}

namespace {
    template <class T>
    void flagop_sub(curvenetwork::Core::FlagOp op, int arg, list<T>& elements) {
        for (auto e = elements.begin(); e != elements.end(); ++e) {
            if (e->is_deleted) continue;
            
            e->flag =
                op == curvenetwork::Core::FlagOp::SUBST  ? arg :
                op == curvenetwork::Core::FlagOp::AND    ? e->flag & arg :
                op == curvenetwork::Core::FlagOp::OR     ? e->flag | arg :
                0;
        }
    }
}

void curvenetwork::Core::set_flag_vertices  (FlagOp op, int arg) { flagop_sub(op, arg, vertices  ); }
void curvenetwork::Core::set_flag_halfedges (FlagOp op, int arg) { flagop_sub(op, arg, halfedges ); }
void curvenetwork::Core::set_flag_halfchains(FlagOp op, int arg) { flagop_sub(op, arg, halfchains); }
void curvenetwork::Core::set_flag_edgechains(FlagOp op, int arg) { flagop_sub(op, arg, edgechains); }
void curvenetwork::Core::set_flag_patches   (FlagOp op, int arg) { flagop_sub(op, arg, patches   ); }

curvenetwork::VertexPtr    curvenetwork::Core::new_vertex   () { vertices  .push_back(Vertex   (vertices  .empty() ? 0 : vertices  .back().id + 1)); return VertexPtr   (vertices  .back()); }
curvenetwork::HalfedgePtr  curvenetwork::Core::new_halfedge () { halfedges .push_back(Halfedge (halfedges .empty() ? 0 : halfedges .back().id + 1)); return HalfedgePtr (halfedges .back()); }
curvenetwork::HalfchainPtr curvenetwork::Core::new_halfchain() { halfchains.push_back(Halfchain(halfchains.empty() ? 0 : halfchains.back().id + 1)); return HalfchainPtr(halfchains.back()); }
curvenetwork::EdgechainPtr curvenetwork::Core::new_edgechain() { edgechains.push_back(Edgechain(edgechains.empty() ? 0 : edgechains.back().id + 1)); return EdgechainPtr(edgechains.back()); }
curvenetwork::PatchPtr     curvenetwork::Core::new_patch    () { patches   .push_back(Patch    (patches   .empty() ? 0 : patches   .back().id + 1)); return PatchPtr    (patches   .back()); }

namespace {
    template <class T>
    map<int, T*> validate_all_ptr_get_map(list<T>& elem_list) {
        map<int, T*> result;
        for (auto elem = elem_list.begin(); elem != elem_list.end(); ++elem)
            result.insert(make_pair(elem->id, &*elem));
        return result;
    }
    template <class T>
    void validate_all_ptr_sub(curvenetwork::ElemPtrT<T>& elem_ptr, const map<int, T*>& elem_map) {
        auto found = elem_map.find(elem_ptr.id);
        elem_ptr.ptr = found == elem_map.end() ? nullptr : found->second;
    }
}

void curvenetwork::Core::validate_all_ptr() {
    auto vertex_map    = validate_all_ptr_get_map(vertices  );
    auto halfedge_map  = validate_all_ptr_get_map(halfedges );
    auto halfchain_map = validate_all_ptr_get_map(halfchains);
    auto edgechain_map = validate_all_ptr_get_map(edgechains);
    auto patch_map     = validate_all_ptr_get_map(patches   );
    
    for (auto& v : vertices) {
        validate_all_ptr_sub(v.halfedge, halfedge_map);
    }
    
    for (auto& h : halfedges) {
        validate_all_ptr_sub(h.vertex   , vertex_map);
        validate_all_ptr_sub(h.next     , halfedge_map);
        validate_all_ptr_sub(h.prev     , halfedge_map);
        validate_all_ptr_sub(h.opposite , halfedge_map);
        validate_all_ptr_sub(h.halfchain, halfchain_map);
    }
    
    for (auto& c : halfchains) {
        validate_all_ptr_sub(c.halfedge_front, halfedge_map );
        validate_all_ptr_sub(c.halfedge_back , halfedge_map );
        validate_all_ptr_sub(c.patch         , patch_map    );
        validate_all_ptr_sub(c.edgechain     , edgechain_map);
        
        if (c.patch.id == -2) c.patch = &Patch::imaginary_patch_symmetry;
        if (c.patch.id == -3) c.patch = &Patch::imaginary_patch_boundary;
    }
    
    for (auto& e : edgechains) {
        validate_all_ptr_sub(e.halfchain[0], halfchain_map);
        validate_all_ptr_sub(e.halfchain[1], halfchain_map);
    }
    
    for (auto& p : patches) {
        validate_all_ptr_sub(p.halfchain, halfchain_map);
        p.set_halfedgeData();
    }
}

void curvenetwork::Core::garbage_collect() {
    vertices  .remove_if([] (Vertex   & v) { return v.is_deleted; });
    halfedges .remove_if([] (Halfedge & h) { return h.is_deleted; });
    halfchains.remove_if([] (Halfchain& c) { return c.is_deleted; });
    edgechains.remove_if([] (Edgechain& e) { return e.is_deleted; });
    patches   .remove_if([] (Patch    & p) { return p.is_deleted; });
}

curvenetwork::Core& curvenetwork::Core::operator=(const Core& rhs) {
    vertices   = rhs.vertices  ;
    halfedges  = rhs.halfedges ;
    halfchains = rhs.halfchains;
    edgechains = rhs.edgechains;
    patches    = rhs.patches   ;
    
    invalidate_displist();
    return *this;
}

void curvenetwork::Core::clear() {
    vertices  .clear();
    halfedges .clear();
    halfchains.clear();
    edgechains.clear();
    patches   .clear();
    
    invalidate_displist();
}
bool curvenetwork::Core::empty() const {
    return vertices.empty();
}
void curvenetwork::Core::invalidate_displist() {
    displist_corner.invalidate();
}
void curvenetwork::Core::render_edgechains() const {
    glLineWidth(2);
    glBegin(GL_LINES);
    for (auto& e : edgechains) {
        if (e.is_deleted) continue;
        
        if (e.num_subdiv == 0)
            glColor3d(1, 0.5, 0);
        else
            glColor3d(1, 0.8, 0);
        
        for (curvenetwork::CHIter h(e.halfchain[0].ptr); h; ++h) {
            auto v0 = h->vertex;
            auto v1 = h->from_vertex();
            if (v0->is_hidden || v1->is_hidden) continue;
            glVertex3dv(&v0->pn[0]);
            glVertex3dv(&v1->pn[0]);
        }
    }
    glEnd();
}
void curvenetwork::Core::render_edgechains_num_subdiv() const {
    glPointSize(10);
    glBegin(GL_POINTS);
    for (auto& e : edgechains) {
        if (e.is_deleted || Patch::is_ordinary(e.halfchain[0]->patch) || Patch::is_ordinary(e.halfchain[1]->patch)) continue;
        
        bool is_hidden = false;
        for (curvenetwork::CHIter h(e.halfchain[0].ptr); h; ++h) {
            is_hidden = h->vertex->is_hidden;
            if (is_hidden) break;
        }
        if (is_hidden) continue;
        
        for (int i = 1; i < e.num_subdiv; ++i) {
            double t = i / static_cast<double>(e.num_subdiv);
            auto pn = e.halfchain[0]->pn_at(t);
            glVertex3dv(&pn[0]);
        }
    }
    glEnd();
}
void curvenetwork::Core::render_corners() {
    displist_corner.render([&] () {
        for (auto& v : vertices) {
            if (v.is_deleted || v.is_hidden || !v.is_corner())
                continue;
            
            Vector3d n = v.pn.tail(3);
            for (curvenetwork::VIHIter h(&v); h; ++h) {
                if (h->is_deleted || !h->next || h->next->is_deleted) continue;
                if (h->is_corner || h->imaginary_patch)
                    glColor4d(1, 0, 0.5, 0.5);
                else
                    glColor4d(0, 0.5, 1, 0.5);
                
                Vector3d d0 = -h      ->toVector3d();
                Vector3d d1 =  h->next->toVector3d();
                eigen_util::orthonormalize(n, d0);
                eigen_util::orthonormalize(n, d1);
                
                double angle = util::acos_clamped(d0.dot(d1));
                if (d1.cross(d0).dot(n) < 0)
                    angle = 2 * util::pi() - angle;
                
                int num_fans = static_cast<int>(angle / (util::pi() / 16)) + 1;
                AngleAxisd rot(angle / num_fans, n);
                Vector3d d_current = d1;
                
                glBegin(GL_TRIANGLE_FAN);
                glTexCoord3d(0, 0, 0);
                glNormal3dv(&v.pn[3]);
                glVertex3dv(&v.pn[0]);
                for (int i = 0; i <= num_fans; ++i, d_current = rot * d_current) {
                    glTexCoord3d(d_current);
                    glVertex3dv(&v.pn[0]);
                }
                glEnd();
            }
        }
    });
}
void curvenetwork::Core::render_endpoints() const {
    glPointSize(7);
    glColor3d(0, 0.5, 0);
    glBegin(GL_POINTS);
    for (auto& v : vertices) {
        if (v.is_deleted || v.is_hidden || !v.is_openend())
            continue;
        glVertex3dv(&v.pn[0]);
    }
    glEnd();
}