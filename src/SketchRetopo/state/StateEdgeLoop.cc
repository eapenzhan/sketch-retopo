#include "../SketchRetopo.hh"
#include <kt84/graphics/phong_tessellation.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateEdgeLoop::StateEdgeLoop()
    : State("edge loop", Vector3d(0, 0.5, 0))
{}

void StateEdgeLoop::init() {
    edgeLoop_preview.clear();
}

void StateEdgeLoop::display() {
    glEnable(GL_STENCIL_TEST);
    glDepthMask(GL_FALSE);
    glStencilFunc(GL_GREATER, ++core.stencil_index, ~0);
    
    glLineWidth(3);
    glColor3d(state_color);
    phong_tessellation::begin(phong_tessellation::Mode::LINE_STRIP);
    for (auto& pn : edgeLoop_preview)
        phong_tessellation::vertex(pn);
    phong_tessellation::end();
    
    glDisable(GL_STENCIL_TEST);
    glDepthMask(GL_TRUE);
}

void StateEdgeLoop::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (!core.common_pn_mouse)
        return;
    Vector3d p_mouse = core.common_pn_mouse->head(3);
    
    init();
    
    // look for boundary Patch::HHandle that is close enough
    curvenetwork::Patch::HHandle h_min;
    curvenetwork::Patch*         patch_min = nullptr;
    double                       dist_min  = util::dbl_max();
    for (auto& patch : core.curvenetwork.patches) {
        if (patch.is_failure)
            continue;
        
        for (auto h : patch.halfedges()) {
            if (!patch.is_boundary(h)) continue;
            
            auto v0v1 = patch.util_halfedge_to_vertex_pair(h);
            Vector3d p0 = patch.data(v0v1.first ).pn.head(3);
            Vector3d p1 = patch.data(v0v1.second).pn.head(3);
            
            double dist = *eigen_util::distance_to_line<Vector3d>(p0, p1, p_mouse, true);
            if (dist < dist_min && dist < core.configTemp.snapSize()) {
                h_min     = h;
                patch_min = &patch;
                dist_min  = dist;
            }
        }
    }
    
    if (!h_min.is_valid())
        return;
    
    Vector2d t;
    auto v0v1 = patch_min->util_halfedge_to_vertex_pair(h_min);
    Vector3d p0 = patch_min->data(v0v1.first ).pn.head(3);
    Vector3d p1 = patch_min->data(v0v1.second).pn.head(3);
    eigen_util::project_to_line(p0, p1, p_mouse, t);
    core.edgeLoop_insert(patch_min, h_min, t[1]);
}

void StateEdgeLoop::mouse_move(int mouse_x, int mouse_y) {
    if (!core.common_pn_mouse)
        return;
    Vector3d p_mouse = core.common_pn_mouse->head(3);
    
    init();
    
    // look for boundary Patch::HHandle that is close enough
    curvenetwork::Patch::HHandle h_min;
    curvenetwork::Patch*         patch_min = nullptr;
    double                       dist_min  = util::dbl_max();
    for (auto& patch : core.curvenetwork.patches) {
        if (patch.is_failure)
            continue;
        
        for (auto h : patch.halfedges()) {
            if (!patch.is_boundary(h)) continue;
            
            auto v0v1 = patch.util_halfedge_to_vertex_pair(h);
            Vector3d p0 = patch.data(v0v1.first ).pn.head(3);
            Vector3d p1 = patch.data(v0v1.second).pn.head(3);
            
            double dist = *eigen_util::distance_to_line<Vector3d>(p0, p1, p_mouse, true);
            if (dist < dist_min && dist < core.configTemp.snapSize()) {
                h_min     = h;
                patch_min = &patch;
                dist_min  = dist;
            }
        }
    }
    
    if (!h_min.is_valid())
        return;
    
    Vector2d t;
    auto v0v1 = patch_min->util_halfedge_to_vertex_pair(h_min);
    Vector3d p0 = patch_min->data(v0v1.first ).pn.head(3);
    Vector3d p1 = patch_min->data(v0v1.second).pn.head(3);
    eigen_util::project_to_line(p0, p1, p_mouse, t);
    edgeLoop_preview = core.edgeLoop_preview(patch_min, h_min, t[1]);
}
