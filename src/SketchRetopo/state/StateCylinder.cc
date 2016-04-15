#include "../SketchRetopo.hh"
#include <kt84/loop_util.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateCylinder::StateCylinder()
    : State("cylinder sketching", Vector3d(0.5, 0.7, 0.3))
{}
void StateCylinder::init() {
    core.common_sketch_curve.clear();
}
void StateCylinder::display() {
    // current curve
    glBegin(GL_LINE_STRIP);
    glColor3d(state_color);
    for (auto& pn : core.common_sketch_curve)
        glVertex3dv(&pn[0]);
    glEnd();
}
void StateCylinder::mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (core.common_sketch_curve.size() < 3) {
        init();
        return;
    }
    
    core.add_cylinder_curve(core.common_sketch_curve);
    
    core.generate_patches();
    
    init();
}
