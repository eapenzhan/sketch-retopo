#include "../SketchRetopo.hh"
#include <kt84/graphics/graphics_util.hh>
#include <boost/range/algorithm.hpp>
#include <kt84/container_util.hh>
#include <kt84/math/RBFKernel.hh>
#include <kt84/MinSelector.hh>
#include <kt84/ScopeExit.hh>
#include "../curvenetwork/Circulator.hh"
#include "../curvenetwork/helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateDeformCurve::StateDeformCurve()
    : State("deform curve", Vector3d(0.3, 0.8, 0.2))
{}
void StateDeformCurve::init() {
    mode            = Mode::None;
    selected_vertex = nullptr;
    core.common_sketch_curve       .clear();
    core.common_affected_edgechains.clear();
}
void StateDeformCurve::display() {
    if (mode == Mode::OverSketch) {
        // oversketched curve
        glBegin(GL_LINE_STRIP);
        glColor3d(state_color);
        for (auto& pn : core.common_sketch_curve)
            glVertex3dv(&pn[0]);
        glEnd();
    }
}
void StateDeformCurve::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button != Button::LEFT || !core.common_pn_mouse)
        return;
    
    // look for closest corner vertex
    MinSelector<curvenetwork::Vertex*> v_corner;
    MinSelector<curvenetwork::Vertex*> v_side;
    for (auto& v : core.curvenetwork.vertices) {
        if (core.configSaved.symmetric && v.on_symmetry_plane())
            continue;
        
        double dist = pn_norm(v.pn - *core.common_pn_mouse);
        if (core.configTemp.snapSize() < dist)
            continue;
        (v.is_corner() || v.is_openend() ? v_corner : v_side).update(dist, &v);
    }
    
    if (v_corner.value) {
        selected_vertex = v_corner.value;
    
    } else if (v_side.value) {
        selected_vertex = v_side.value;
    
    } else {
        init();
        return;
    }
    
    mode = ctrl_pressed ? Mode::OverSketch : (selected_vertex->is_corner() || selected_vertex->is_openend()) ? Mode::DragCorner : Mode::DragSide;
    
    if (mode == Mode::DragCorner || mode == Mode::DragSide) {
        for (curvenetwork::VOCIter c(selected_vertex); c; ++c)
            core.common_affected_edgechains.push_back(c->edgechain);
        
        prev_pos = selected_vertex->pn.head(3);
    }
    
    core.memento_store();
}

void StateDeformCurve::mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (mode == Mode::None)
        return;
    
    if (mode == Mode::OverSketch) {
        if (!core.common_pn_mouse) {
            init();
            return;
        }
        
        // look for edgechain, where one of its endpoints is selected_vertex the other is snapped to by *core.common_pn_mouse, and which is closest to core.common_sketch_curve
        MinSelector<curvenetwork::Vertex*> other_vertex;
        for (auto& v : core.curvenetwork.vertices) {
            if (core.configSaved.symmetric && v.on_symmetry_plane())
                continue;
            
            if (&v == selected_vertex)
                continue;
            
            double dist = pn_norm(v.pn - *core.common_pn_mouse);
            other_vertex.update(dist, &v);
        }
        
        if (other_vertex.score > core.configTemp.snapSize()) {
            // no close vertex found
            init();
            return;
        }
        
        // find a sequence of halfedges connecting between the two vertices, without meetin any other corner vertex
        curvenetwork::Halfedge* h_found = nullptr;
        int sequence_length;
        for (curvenetwork::VOHIter h(selected_vertex); h; ++h) {
            sequence_length = 1;
            for (auto h2 = &*h; ; h2 = h2->next, ++sequence_length) {
                auto v = h2->vertex;
                if (v == other_vertex.value) {
                    h_found = &*h;
                    break;
                }
                if (v->is_corner() || v->is_openend())
                    break;
            }
            if (h_found)
                break;
        }
        if (!h_found) {
            init();
            return;
        }
        
        auto h = h_found;
        for (int i = 1; i < sequence_length; ++i, h = h->next) {
            double t = i / static_cast<double>(sequence_length);
            h->vertex->pn = core.common_sketch_curve.point_at(t);
        }
        
        core.common_affected_edgechains.push_back(h->halfchain->edgechain);
    }
    
    // collect affected patches
    vector<curvenetwork::Patch*> affected_patches;
    container_util::remove_duplicate(core.common_affected_edgechains);
    for (auto e : core.common_affected_edgechains) {
        for (int i = 0; i < 2; ++i) {
            auto patch = e->halfchain[i]->patch;
            if (!patch || patch->is_imaginary())
                continue;
            affected_patches.push_back(patch);
        }
    }
    container_util::remove_duplicate(affected_patches);
    for (auto patch : affected_patches)
        if (!patch->is_deleted) patch->clear();         // NOTE: patch parameter is kept intact!
    
    auto scope_exit = make_ScopeExit([&](){
        // update patch vertex positions
        for (auto patch : affected_patches) {
            if (patch->is_deleted) continue;
            patch->generate_topology(false);
            core.compute_patch_interior_pn(patch);
            patch->invalidate_displist();
        }
        core.curvenetwork.invalidate_displist();
        
        core.common_affected_edgechains.clear();
        core.curvenetwork.garbage_collect();
        
        init();
    });
    
    // no snapping when selected_vertex is neither corner nor openend
    if (!selected_vertex->is_corner() && !selected_vertex->is_openend()) return;
    
    // snap selected vertex to corner/side
    MinSelector<curvenetwork::Vertex*> v_openend;
    MinSelector<curvenetwork::Vertex*> v_corner;
    MinSelector<curvenetwork::Vertex*> v_corner_adjacent;
    MinSelector<curvenetwork::Vertex*> v_side;
    vector<curvenetwork::Edgechain*> e_adjacent;
    for (auto c = curvenetwork::VOCIter(selected_vertex); c; ++c)
        e_adjacent.push_back(c->edgechain);
    for (auto& v : core.curvenetwork.vertices) {
        if (&v == selected_vertex) continue;
        
        double dist = pn_norm(selected_vertex->pn - v.pn);
        if (dist > core.configTemp.snapSize()) continue;
        
        if (v.is_openend()) {
            if (v.halfedge->halfchain->vertex_back() == selected_vertex) continue;      // degenerate!
            v_openend.update(dist, &v);
        
        } else if (v.is_corner()) {
            if (selected_vertex->is_openend() && selected_vertex->halfedge->halfchain->vertex_back() == &v) continue;       // degenerate!
            bool is_adjacent = false;
            for (auto c = curvenetwork::VOCIter(selected_vertex); c; ++c) {
                if (c->halfedge_back->vertex == &v) {
                    is_adjacent = true;
                    break;
                }
            }
            (is_adjacent ? v_corner_adjacent : v_corner).update(dist, &v);
        
        } else if (boost::range::find(e_adjacent, v.halfedge->halfchain->edgechain) == e_adjacent.end()) {
            v_side.update(dist, &v);
        }
    }
    
    auto v_snapped = v_openend.value         ? v_openend.value         :
                     v_corner_adjacent.value ? v_corner_adjacent.value :
                     v_corner.value          ? v_corner.value          :
                     v_side.value;
    if (!v_snapped) return;
    
    // offset to vertices on edgechanins emanating from selected_vertex
    Vector3d offset_value = (v_snapped->pn - selected_vertex->pn).head(3);
    auto offset_func = [&] (PointNormal& pn, double w) {
        pn.head(3) += w * offset_value;
        core.project(pn);
    };
    for (curvenetwork::VOCIter c(selected_vertex); c; ++c) {
        int n = c->toPolyline().size() - 1;
        auto h = c->halfedge_front;
        for (int i = 1; i < n; ++i, h = h->next) {
            // linear weight
            double w = 1.0 - i / static_cast<double>(n);
            offset_func(h->vertex->pn, w);
        }
    }
    
    if (v_snapped == v_corner_adjacent.value) {
        // find halfchain from selected_vertex to v_snapped
        curvenetwork::Halfchain* c_deleted = nullptr;
        for (auto c = curvenetwork::VOCIter(selected_vertex); c; ++c) {
            if (c->vertex_back() == v_snapped) {
                c_deleted = &*c;
                break;
            }
        }
        // reconnect vertex, delete unused one
        for (auto h = curvenetwork::VIHIter(selected_vertex); h; ++h)
            h->vertex = v_snapped;
        selected_vertex->is_deleted = true;
        // store halfedges that will be connected
        auto h0_prev_next = make_pair(c_deleted            ->halfedge_front->prev, c_deleted            ->halfedge_back->next);
        auto h1_prev_next = make_pair(c_deleted->opposite()->halfedge_front->prev, c_deleted->opposite()->halfedge_back->next);
        // delete halfedges
        auto h_deleted_front = c_deleted->halfedge_front;
        auto h_deleted_back  = c_deleted->halfedge_back;
        for (auto h = h_deleted_front; ; h = h->next) {
            curvenetwork::delete_halfedge(h->opposite);
            curvenetwork::delete_halfedge(h);
            if (h == h_deleted_back) break;
        }
        // connect halfedges
        curvenetwork::set_next_prev(h0_prev_next.first, h0_prev_next.second);
        curvenetwork::set_next_prev(h1_prev_next.first, h1_prev_next.second);
    
    } else {
        // get boundary halfedges around snapped vertices
        auto get_booundary_halfedges = [] (curvenetwork::Vertex* v, curvenetwork::Halfedge*& h_prev, curvenetwork::Halfedge*& h_next) {
            h_prev = h_next = nullptr;
            for (auto h = curvenetwork::VIHIter(v); h; ++h) {
                if (h->patch()) continue;
                if (h_prev) return false;
                h_prev = &*h;
            }
            for (auto h = curvenetwork::VOHIter(v); h; ++h) {
                if (h->patch()) continue;
                if (h_next) return false;
                h_next = &*h;
            }
            return true;
        };
        curvenetwork::Halfedge* hsrc_prev, *hsrc_next;
        curvenetwork::Halfedge* htgt_prev, *htgt_next;
        if (!get_booundary_halfedges(selected_vertex, hsrc_prev, hsrc_next)) return;
        if (!get_booundary_halfedges(v_snapped      , htgt_prev, htgt_next)) return;
        
        // reconnect vertex, delete unused one
        vector<curvenetwork::Halfedge*> hsrc_vih;
        for (auto h = curvenetwork::VIHIter(selected_vertex); h; ++h)
            hsrc_vih.push_back(&*h);
        for (auto h : hsrc_vih)
            h->vertex = v_snapped;
        selected_vertex->is_deleted = true;
        
        if (v_snapped == v_side.value) {
            // delete halfchain
            curvenetwork::delete_halfchain(htgt_prev->halfchain);
            curvenetwork::delete_halfchain(htgt_prev->opposite->halfchain);
        }
        
        curvenetwork::set_next_prev(hsrc_prev, htgt_next);
        curvenetwork::set_next_prev(htgt_prev, hsrc_next);
        hsrc_prev->is_corner = true;
        htgt_prev->is_corner = true;
        
        if (v_snapped == v_side.value) core.curvenetwork.generate_halfchains();
        
        auto csrc_prev = hsrc_prev->halfchain;
        auto csrc_next = hsrc_next->halfchain;
        auto ctgt_prev = htgt_prev->halfchain;
        auto ctgt_next = htgt_next->halfchain;
        
        // remove any redundant curves
        if (!csrc_next->is_deleted && csrc_next->vertex_back() == ctgt_prev->vertex_front()) {
            // set next/prev at endpoints of deleted halfchain
            curvenetwork::set_next_prev(htgt_prev, hsrc_next->opposite->next);
            curvenetwork::set_next_prev(csrc_next->halfedge_back->opposite->prev, csrc_next->halfedge_back->next);
            auto h_deleted_front = csrc_next->halfedge_front;
            auto h_deleted_back  = csrc_next->halfedge_back;
            for (auto h = h_deleted_front; ; h = h->next) {
                curvenetwork::delete_halfedge(h->opposite);
                curvenetwork::delete_halfedge(h);
                if (h == h_deleted_back) break;
            }
        }
        if (!csrc_prev->is_deleted && csrc_prev->vertex_front() == ctgt_next->vertex_back()) {
            // set next/prev at endpoints of deleted halfchain
            curvenetwork::set_next_prev(hsrc_prev->opposite->prev, htgt_next);
            curvenetwork::set_next_prev(csrc_prev->halfedge_front->prev, csrc_prev->halfedge_front->opposite->next);
            auto h_deleted_front = csrc_prev->halfedge_front;
            auto h_deleted_back  = csrc_prev->halfedge_back;
            for (auto h = h_deleted_front; ; h = h->next) {
                curvenetwork::delete_halfedge(h->opposite);
                curvenetwork::delete_halfedge(h);
                if (h == h_deleted_back) break;
            }
        }
    }
    
    core.generate_patches();
}
void StateDeformCurve::mouse_move(int mouse_x, int mouse_y) {
    if (!core.common_pn_mouse || mode == Mode::None)
        return;
    
    if (mode == Mode::DragCorner || mode == Mode::DragSide) {
        // offset value and functor
        Vector3d offset_value = core.common_pn_mouse->head(3) - prev_pos;
        auto offset_func = [&] (PointNormal& pn, double w) {
            pn.head(3) += w * offset_value;
            core.project(pn);
        };
        
        // apply offset to vertices connected to selected_vertex
        if (mode == Mode::DragCorner) {
            for (curvenetwork::VOCIter c(selected_vertex); c; ++c) {
                int n = c->toPolyline().size() - 1;
                auto h = c->halfedge_front;
                for (int i = 1; i < n; ++i, h = h->next) {
                    // linear weight
                    double w = 1.0 - i / static_cast<double>(n);
                    offset_func(h->vertex->pn, w);
                }
            }
            
        } else {
            auto h_start = curvenetwork::VOHIter(selected_vertex);
            for (int i = 0; i < 2; ++i, ++h_start) {
                int num_halfedges = 1;
                for (auto h = &*h_start; ; ) {
                    if (h->vertex->is_corner() || h->vertex->is_openend())
                        break;
                    ++num_halfedges;
                    h = h->next;
                    if (h == &*h_start)
                        return;
                }
                auto h = &*h_start;
                for (int j = 1; j < num_halfedges; ++j) {
                    // smooth weight
                    double r = j / static_cast<double>(num_halfedges);
                    double w = RBFKernel_Wendland()(j / static_cast<double>(num_halfedges));
                    offset_func(h->vertex->pn, w);
                    h = h->next;
                }
            }
        }
        // apply offset to selected_vertex
        offset_func(selected_vertex->pn, 1);
        
        prev_pos = core.common_pn_mouse->head(3);
    }
    
    if (mode == Mode::DragCorner)
        core.curvenetwork.invalidate_displist();
}
