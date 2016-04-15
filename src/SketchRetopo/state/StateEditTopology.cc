#include "../SketchRetopo.hh"
#include <kt84/graphics/graphics_util.hh>
#include <kt84/eigen_util.hh>
#include <kt84/loop_util.hh>
#include <kt84/container_util.hh>
#include <kt84/MinSelector.hh>
#include <patchgen/decl.hh>
#include "../curvenetwork/Circulator.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateEditTopology::StateEditTopology()
    : State("edit topology", Vector3d(1, 0, 0))
{}

void StateEditTopology::init() {
    current_patch = nullptr;
    selected_variable = 0;
    core.curvenetwork.set_flag_edgechains(curvenetwork::Core::FlagOp::SUBST, 0);
}

void StateEditTopology::display() {
    glEnable(GL_STENCIL_TEST);
    glDepthMask(GL_FALSE);
    
    // draw currently selected patch
    if (current_patch) {
        glStencilFunc(GL_GREATER, ++core.stencil_index, ~0);
        glColor4d(state_color, 0.3);
        current_patch->render_phong_fill();
        
        // draw arrows indicating adjustable variables
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_STENCIL_TEST);
        glLineWidth(2);
        auto draw_doublesided_arrow = [&] (Vector3d p0, Vector3d p1, Vector3d n) {
            double arrow_size = (0.5 * p0 + 0.5 * p1 - core.camera->get_eye()).norm() * 0.01;
            Vector3d d = p1 - p0;
            p0 += 0.1 * d;
            p1 -= 0.1 * d;
            double r = d.norm();
            arrow_size = min<double>(arrow_size, r * 0.3);
            d *= arrow_size / r;
            glBegin(GL_LINES);
            glVertex3d(p0 + d);
            glVertex3d(p1 - d);
            glEnd();
            Vector3d e = 0.4 * n.cross(d);
            glBegin(GL_TRIANGLES);
            glVertex3d(p0);    glVertex3d(p0 + d - e);    glVertex3d(p0 + d + e);
            glVertex3d(p1);    glVertex3d(p1 - d + e);    glVertex3d(p1 - d - e);
            glEnd();
        };
        auto& param = current_patch->param;
        int num_sides = param.get_num_sides();
        // draw arrows representing padding parameters only when pattern permits variable adjustment
        int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
        if ( num_sides != 4 && num_sides < num_variables) {
            for (int i = 0; i < num_sides; ++i) {
                int j = param.permutation[i];
                auto polyline = current_patch->get_side_polyline(j);
                PointNormal pn = polyline.point_at(0.5);
                Vector3d n = pn.tail(3);
                Vector3d d = (polyline.point_at(0.51) - pn).head(3);
                d = n.cross(d).normalized() * average_edge_length * 0.5;
                Vector3d p0 = pn.head(3) - d;
                Vector3d p1 = pn.head(3) + d;
                glColor3d(i == selected_variable ? Vector3d(1, 0, 0) : Vector3d(0.5, 0, 0));
                draw_doublesided_arrow(p0, p1, n);
            }
        }
        const Vector3d color_table[4] = {
            Vector3d(0.6, 0.4, 0),
            Vector3d(0, 0.5, 0),
            Vector3d(0, 0.5, 0.6),
            Vector3d(0.5, 0, 0.5),
        };
        for (int i = num_sides; i < num_variables; ++i) {
            auto& variable_indicators = patchgen::get_variable_indicators(num_sides, param.pattern_id)[i - num_sides];
            for (auto& tag_pair : variable_indicators) {
                auto v0 = patchgen::find_tagged_vertex(*current_patch, tag_pair.first);
                auto v1 = patchgen::find_tagged_vertex(*current_patch, tag_pair.second);
                assert(v0.is_valid() && v1.is_valid());
                glColor3d((i == selected_variable ? 2 : 1) * color_table[i - num_sides]);
                PointNormal& pn0 = current_patch->data(v0).pn;
                PointNormal& pn1 = current_patch->data(v1).pn;
                draw_doublesided_arrow(pn0.head(3), pn1.head(3), (pn0 + pn1).tail(3).normalized());
            }
        }
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_STENCIL_TEST);
    }
    
    // draw currently selected halfchains
    glStencilFunc(GL_GREATER, ++core.stencil_index, ~0);
    glLineWidth(5);
    glColor4d(state_color, 0.5);
    
    for (auto& e : core.curvenetwork.edgechains) {
        if (!e.flag) continue;
        
        auto c = e.halfchain[0];
        glBegin(GL_LINE_STRIP);
        glVertex3dv(&c->halfedge_front->from_vertex()->pn[0]);
        for (curvenetwork::CHIter h(c); h; ++h) {
            auto& pn = h->vertex->pn;
            glVertex3d(pn.head(3) + pn.tail(3) * core.configSaved.projOffset * 0.1);
        }
        glEnd();
    }
    
    // visualize num_subdiv with points
    glStencilFunc(GL_GREATER, ++core.stencil_index, ~0);
    glColor3d(0.1, 0.2, 1.0);
    core.curvenetwork.render_edgechains_num_subdiv();
    
    glDisable(GL_STENCIL_TEST);
    glDepthMask(GL_TRUE);
}
void StateEditTopology::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (!core.common_pn_mouse) {
        init();
        return;
    }
    Vector3d mouse_pos = core.common_pn_mouse->head(3);
    
    if (current_patch) {
        auto& param = current_patch->param;
        int num_sides = param.get_num_sides();
        if (num_sides == 4) {
            // 4-sided patch -> edit topology by either dragging the singularity around, or adjusting param.x
            if (current_patch->distance(mouse_pos) < core.configTemp.snapSize()) {
                p_prev = mouse_pos;
                core.memento_store();
            } else {
                init();
            }
        } else {
            // select current patch's parameter
            MinSelector<int> selector;
            for (int i = 0; i < num_sides; ++i) {
                int j = param.permutation[i];
                auto polyline = current_patch->get_side_polyline(j);
                PointNormal pn = polyline.point_at(0.5);
                Vector3d n = pn.tail(3);
                Vector3d d = (polyline.point_at(0.51) - pn).head(3);
                d = n.cross(d).normalized() * average_edge_length * 0.5;
                Vector3d p0 = pn.head(3) - d;
                Vector3d p1 = pn.head(3) + d;
                auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
                selector.update(*dist, i);
            }
            int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
            for (int i = num_sides; i < num_variables; ++i) {
                auto& variable_indicators = patchgen::get_variable_indicators(num_sides, param.pattern_id)[i - num_sides];
                for (auto& tag_pair : variable_indicators) {
                    auto v0 = patchgen::find_tagged_vertex(*current_patch, tag_pair.first);
                    auto v1 = patchgen::find_tagged_vertex(*current_patch, tag_pair.second);
                    assert(v0.is_valid() && v1.is_valid());
                    Vector3d p0 = current_patch->data(v0).pn.head(3);
                    Vector3d p1 = current_patch->data(v1).pn.head(3);
                    auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
                    selector.update(*dist, i);
                }
            }
            if (selector.score < core.configTemp.snapSize()) {
                selected_variable = selector.value;
            } else {
                init();
            }
        }
        return;
    }
    
    // select closest edgechain
    MinSelector<curvenetwork::Edgechain*> e_closest;
    for (auto& v : core.curvenetwork.vertices) {
        if (v.is_corner()) continue;
        double dist = pn_norm(v.pn - *core.common_pn_mouse);
        if (dist > core.configTemp.snapSize()) continue;
        e_closest.update(dist, v.halfedge->halfchain->edgechain);
    }
    if (e_closest.value) {
        if (shift_pressed) {
            // toggle selection/de-selection
            e_closest.value->flag = e_closest.value->flag ? 0 : 1;
            
        } else if (ctrl_pressed) {
            // select sequence of halfchains on the same quad strip
            auto patch = e_closest.value->halfchain[curvenetwork::Patch::is_ordinary(e_closest.value->halfchain[0]->patch) ? 0 : 1]->patch;
            if (!curvenetwork::Patch::is_ordinary(patch))
                return;
            
            // look for patch halfedge closest to the mouse
            MinSelector<curvenetwork::Patch::HHandle> h_selector;
            for (auto h : patch->halfedges()) {
                if (!patch->is_boundary(h))
                    continue;
                
                auto v0v1 = patch->util_halfedge_to_vertex_pair(h);
                Vector3d p0 = patch->data(v0v1.first ).pn.head(3);
                Vector3d p1 = patch->data(v0v1.second).pn.head(3);
                
                double dist = *eigen_util::distance_to_line<Vector3d>(p0, p1, mouse_pos, true);
                h_selector.update(dist, h);
            }
            auto patch_halfedge_pairs = core.edgeLoop_walk(patch, h_selector.value);
            
            for (auto& patch_halfedge : patch_halfedge_pairs) {
                auto patch = patch_halfedge.first;
                auto h     = patch_halfedge.second;
                
                if (!patch->is_boundary(h))
                    continue;
                
                patch->data(h).halfchain->edgechain->flag = 1;
            }
            
        } else {
            // select currently clicked edgechain only
            init();
            e_closest.value->flag = 1;
        }
        return;
    }
    
    init();
    
    // select current patch
    MinSelector<curvenetwork::Patch*> patch_selector;
    for (auto& patch : core.curvenetwork.patches) {
        int num_sides  = patch.param.get_num_sides();
        int pattern_id = patch.param.pattern_id;
        
        // skip if 4-sided regular grid
        if (num_sides == 4 && pattern_id == 0)
            continue;
        
        // bounding box test
        AlignedBox3d bbox;
        for (auto v : patch.vertices())
            bbox.extend(patch.data(v).pn.head(3));
        Vector3d snap_margin = Vector3d::Constant(core.configTemp.snapSize());
        bbox.extend(bbox.max() + snap_margin);
        bbox.extend(bbox.min() - snap_margin);
        if (!bbox.contains(mouse_pos))
            continue;
        
        patch_selector.update(patch.distance(mouse_pos), &patch);
    }
    if (patch_selector.score < core.configTemp.snapSize()) {
        current_patch = patch_selector.value;
        
        selected_variable = current_patch->param.get_num_sides();
        
        average_edge_length = 0;
        for (auto e : current_patch->edges()) {
            auto v = current_patch->util_edge_to_vertex_pair(e);
            Vector3d p0 = current_patch->data(v.first ).pn.head(3);
            Vector3d p1 = current_patch->data(v.second).pn.head(3);
            average_edge_length += (p1 - p0).norm();
        }
        average_edge_length /= current_patch->n_edges();
    }
}
void StateEditTopology::mouse_move(int mouse_x, int mouse_y) {
    if (!current_patch || !core.common_pn_mouse || !core.common_button_left_pressed) return;
    
    auto& param = current_patch->param;
    if (param.get_num_sides() != 4) return;
    
    Vector3d d = core.common_pn_mouse->head(3) - p_prev;
    if (d.norm() < average_edge_length) return;
    d.normalize();
    
    Vector3d p[4], e[4];
    for (int i = 0; i < 4; ++i)
        p[i] = current_patch->data(current_patch->patch_corner(i)).pn.head(3);
    e[0] = (p[3] - p[0] + p[2] - p[1]).normalized();
    e[1] = (p[0] - p[1] + p[3] - p[2]).normalized();
    e[2] = -e[0];
    e[3] = -e[1];
    MinSelector<int> selector;
    for (int i = 0; i < 4; ++i)
        selector.update(-d.dot(e[i]), i);
    
    if (adjust_parameter(param, param.permutation[selector.value], param.permutation.is_flipped())) {         // super confusing... is this correct?
        current_patch->generate_topology(false);
        core.compute_patch_interior_pn(current_patch);
    }
    p_prev = core.common_pn_mouse->head(3);
}
void StateEditTopology::keyboard(unsigned char key, int x, int y) {
    if (current_patch) {
        if (!(key == 'z' || key == 'x' || key == 'c' || key == 'v' || key == 'b' || key == 'B')) return;
        core.memento_store();
        bool do_update = false;
        if (key == 'z' || key == 'x' || key == 'c' || key == 'v')
            do_update = adjust_parameter(current_patch->param, selected_variable, key == 'c' || key == 'v');
        if (key == 'b' || key == 'B')
            do_update = switch_pattern(current_patch->param, key == 'b');
        if (do_update) {
            current_patch->generate_topology(false);
            core.compute_patch_interior_pn(current_patch);
        }
        return;
    }
    
    if (!(key == 'z' || key == 'x' || key == 'c' || key == 'v')) return;
    core.memento_store();
    
    vector<curvenetwork::Patch*> updated_patches;
    for (auto& e : core.curvenetwork.edgechains) {
        if (!e.flag) continue;
        
        int num_subdiv_old = e.num_subdiv;
        
        e.num_subdiv +=
            key == 'z' ? -1 :
            key == 'v' ?  1 :
            key == 'x' ? -2 : 2;
        
        if (e.num_subdiv < 1) {
            e.num_subdiv = num_subdiv_old;
            continue;
        }
        
        for (int i = 0; i < 2; ++i) {
            auto p = e.halfchain[i]->patch;
            if (!p || p->is_imaginary())
                continue;
            
            updated_patches.push_back(p);
            p->clear();
        }
    }
    
    container_util::remove_duplicate(updated_patches);
    
    for (auto p : updated_patches) {
        p->generate_topology();
        core.compute_patch_interior_pn(p);
    }
}
