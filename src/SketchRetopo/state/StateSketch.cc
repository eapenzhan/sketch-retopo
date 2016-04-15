#include "../SketchRetopo.hh"
#include <kt84/container_util.hh>
#include <kt84/graphics/graphics_util.hh>
#include "../curvenetwork/helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::container_util;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateSketch::StateSketch()
    : State("basic sketching", Vector3d(1, 0.5, 0.1))
{}

void StateSketch::init() {
    core.common_sketch_curve.clear();
}

void StateSketch::display() {
    // current curve
    glBegin(GL_LINE_STRIP);
    glColor3d(state_color);
    for (auto& pn : core.common_sketch_curve)
        glVertex3dv(&pn[0]);
    glEnd();
}

void StateSketch::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (!core.common_pn_mouse) return;
    
    if (core.common_key_pressed['z'])       // fix weird problem around clicked area
    {
        // collect halfedges crossing the snapping sphere outward
        vector<pair<double, curvenetwork::Halfedge*>> crossing_halfedges;       // pair<angle, halfedge> (angle is used later for sorting)
        for (auto& h : core.curvenetwork.halfedges) {
            double dist0 = pn_norm(*core.common_pn_mouse - h.from_vertex()->pn);
            double dist1 = pn_norm(*core.common_pn_mouse - h.vertex       ->pn);
            double dist_threshold = core.configTemp.snapSize();
            if (dist0 < dist_threshold && dist_threshold < dist1)
                crossing_halfedges.push_back(make_pair(0.0, &h));
            // erase halfedges that are completely inside snapping radius
            if (dist0 < dist_threshold && dist1 < dist_threshold)
                curvenetwork::delete_halfedge(&h);
        }
    
        if (crossing_halfedges.empty()) return;
        core.memento_store();
    
        // make new vertex at mouse pos
        auto v_center = core.curvenetwork.new_vertex();
        v_center->pn = *core.common_pn_mouse;
    
        // sort halfedge according to angle with respect to the center normal and first halfedge
        Vector3d d0 = (crossing_halfedges[0].second->from_vertex()->pn - v_center->pn).head(3);
        Vector3d n = v_center->pn.tail(3);
        Vector3d d1 = n.cross(d0);
        for (auto& p : crossing_halfedges) {
            Vector3d d = (p.second->from_vertex()->pn - v_center->pn).head(3);
            double angle = atan2(d1.dot(d), d0.dot(d));
            if (angle < 0) angle += 2 * util::pi();
            p.first = angle;
        }
        boost::range::sort(crossing_halfedges);
    
        // make new halfedges
        vector<curvenetwork::Halfedge*> new_halfedges;
        for (auto& p : crossing_halfedges) {
            auto h0 = core.curvenetwork.new_halfedge();
            auto h1 = core.curvenetwork.new_halfedge();
            new_halfedges.push_back(h0);
            curvenetwork::set_vertex_halfedge(v_center, p.second->from_vertex(), h0, h1);
            curvenetwork::set_next_prev(h0, p.second);
            curvenetwork::set_next_prev(p.second->opposite, h1);
            h1->is_corner = true;
        }
    
        // set next/prev between new halfedges
        for (int i = 0; i < crossing_halfedges.size(); ++i) {
            auto h0 = new_halfedges[i];
            auto h1 = at_mod(new_halfedges, i + 1);
            curvenetwork::set_next_prev(h1->opposite, h0);
        }
    
        core.curvenetwork.garbage_collect();
        core.curvenetwork.generate_halfchains();
        core.generate_patches();
        core.curvenetwork.invalidate_displist();
    }
}
void StateSketch::mouse_up(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (core.common_sketch_curve.size() < 3) {
        init();
        return;
    }
    
    core.memento_store();
    
    if (core.common_key_pressed['c'])           // erase parts of curve network crossing with scribble
    {
        core.curvenetwork.erase_curve(core.common_sketch_curve);
    }
    else if (core.common_key_pressed['x'])      // manually generate patch if the starting point of the stroke is bounded by loop of halfedges one of which is snapped to by the ending point of the stroke
    {
        // TODO

    } else {
        auto is_negative = [] (PointNormal& pn) { return pn.x() < 0; };
        bool has_negative = find_if(core.common_sketch_curve.begin(), core.common_sketch_curve.end(), is_negative) != core.common_sketch_curve.end();
        
        if (core.configSaved.symmetric && has_negative) {
            // in symmetry mode, remove all points on the negative x side from the curve
            if (core.common_sketch_curve.is_loop) {
                // make sure that back and front are on the negative and positive x sides, respectively
                auto pos = core.common_sketch_curve.begin();
                for ( ; pos != core.common_sketch_curve.end(); ++pos) {
                    if (is_negative(*pos) && !is_negative(*(pos + 1)))
                        break;
                }
                rotate(core.common_sketch_curve.begin(), pos, core.common_sketch_curve.end());
            }
            
            // split curve into pieces of open curves
            vector<Polyline_PointNormal> split_curves;
            bool is_new = true;
            for (auto& pn : core.common_sketch_curve) {
                if (is_negative(pn)) {
                    is_new = true;
                    continue;
                }
                if (is_new)
                    split_curves.push_back(Polyline_PointNormal());
                is_new = false;
                split_curves.back().push_back(pn);
            }
            
            // add split curves (which are all open strokes)
            for (auto& split_curve : split_curves)
                core.add_curve_open_delegate(split_curve, !ctrl_pressed);
            
        } else if (core.common_sketch_curve.is_loop) {
            // add closed loop
            core.curvenetwork.add_curve_loop(core.common_sketch_curve);
            
        } else {
            // add open stroke
            core.add_curve_open_delegate(core.common_sketch_curve, !ctrl_pressed);
        }
    }
    
    core.generate_patches();
    
    init();
}
void StateSketch::keyboard(unsigned char key, int x, int y) {
    if (key == 'g') {
        // erase all open-ended curves
        vector<curvenetwork::Edgechain*> erased_curves;
        for (auto& v : core.curvenetwork.vertices) {
            if (v.is_openend() && !v.on_symmetry_plane())
                erased_curves.push_back(v.halfedge->halfchain->edgechain);
        }
        
        container_util::remove_duplicate(erased_curves);
        if (erased_curves.empty())
            return;
        
        core.memento_store();
        
        vector<curvenetwork::Halfedge*> erased_halfedges;
        for (auto e : erased_curves)
            erased_halfedges.push_back(e->halfchain[0]->halfedge_front);
        
        for (auto h : erased_halfedges)
            core.curvenetwork.erase_curve_sub(h);
        
        core.generate_patches();
    }
}
