#include "../SketchRetopo.hh"
#include <kt84/container_util.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateLaser::StateLaser()
    : State("laser sketching", Vector3d(0.5, 0, 1))
{}

void StateLaser::init() {
    pn_prev = pn_current = boost::none;
}
void StateLaser::display() {
    if (!pn_prev || !pn_current)
        return;
    
    glDisable(GL_DEPTH_TEST);
    glColor3d(state_color);
    glLineWidth(5);
    glBegin(GL_LINES);
    glVertex3dv(&(*pn_prev   )[0]);
    glVertex3dv(&(*pn_current)[0]);
    glEnd();
}
void StateLaser::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (!core.common_pn_mouse)
        return;
    pn_prev = *core.common_pn_mouse;
}
void StateLaser::mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (!pn_prev || !pn_current) {
        init();
        return;
    }
    
    // laser-cutting plane, depending on whether using perspective/ortho projection
    Vector3d plane_normal;
    if (core.configRender.use_ortho) {
        Vector3d d1 = core.camera->eye_to_center().normalized();
        Vector3d d2 = (*pn_current - *pn_prev).head(3);
        eigen_util::orthonormalize(d1, d2);
        plane_normal = d1.cross(d2);
        
    } else {
        Vector3d plane_e0 = pn_current->head(3) - core.camera->get_eye();
        Vector3d plane_e1 = pn_prev   ->head(3) - core.camera->get_eye();
        plane_normal = plane_e0.cross(plane_e1).normalized();
    }
    
    auto traced_curve = core.trace_on_plane(*pn_prev, plane_normal, pn_norm(*pn_current - *pn_prev) * 10);
    
    if (traced_curve.empty()) {
        init();
        return;
    }
    
    core.memento_store();
    
    if (core.configSaved.symmetric) {
        // check whether curve crosses symmetry plane -> if yes, trim it
        if (traced_curve.is_loop) {
            int rotate_pos = -1;
            for (int i = 0; i < traced_curve.size(); ++i) {
                int i_prev = i == 0 ? traced_curve.size() - 1 : i - 1;
                if (traced_curve[i_prev].x() < 0 && traced_curve[i].x() >= 0) {
                    rotate_pos = i;
                    break;
                }
            }
            if (rotate_pos != -1) {
                rotate(traced_curve.begin(), traced_curve.begin() + rotate_pos, traced_curve.end());
                traced_curve.is_loop = false;
            }
        }
        // erase all points with negative x
        container_util::remove_if(traced_curve, [] (PointNormal& pn) { return pn.x() < 0; });
    }
    
    if (traced_curve.is_loop)
        core.curvenetwork.add_curve_loop(traced_curve);
    else
        core.add_curve_open_delegate(traced_curve);
    
    core.generate_patches();

    init();
}
void StateLaser::mouse_move(int mouse_x, int mouse_y) {
    if (pn_prev && core.common_pn_mouse)
        pn_current = *core.common_pn_mouse;
}
