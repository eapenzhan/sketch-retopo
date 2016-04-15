#include "../SketchRetopo.hh"
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateSpine::StateSpine()
    : State("spine sketching", Vector3d(0.4, 0.5, 1))
{}

void StateSpine::init() {
    core.common_sketch_curve.clear();
}
void StateSpine::display() {
    // current curve
    glBegin(GL_LINE_STRIP);
    glColor3d(state_color);
    for (auto& pn : core.common_sketch_curve)
        glVertex3dv(&pn[0]);
    glEnd();
    
}
void StateSpine::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_spine.adjust_init(mouse_x);
        return;
    }
}
void StateSpine::mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_spine.adjust_done();
        return;
    }
    
    // in symmetric mode, curve should not pass negative x side
    bool is_symmetric_negative =
        core.configSaved.symmetric &&
        find_if(core.common_sketch_curve.begin(), core.common_sketch_curve.end(), [] (PointNormal& pn) { return pn.x() < 0; }) != core.common_sketch_curve.end();
    
    if (core.common_sketch_curve.size() < 3 || core.common_sketch_curve.is_loop || is_symmetric_negative) {
        init();
        return;
    }
    
    core.memento_store();
    
    core.add_spine_curve(core.common_sketch_curve);
    
    core.generate_patches();
    
    init();
}
