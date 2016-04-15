#include "../SketchRetopo.hh"
#include <kt84/graphics/phong_tessellation.hh>
#include <kt84/graphics/graphics_util.hh>
#include <kt84/MinSelector.hh>
#include <kt84/glut_util.hh>
using namespace std;
using namespace kt84;
using namespace kt84::graphics_util;
using namespace Eigen;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateAutocmpl::StateAutocmpl()
    : State("auto-completion", Vector3d(0.5, 1, 0.5))
{}

void StateAutocmpl::init() {
    suggested_autocmpl.clear();
    selected_autocmpl_index = -1;
    h_selected              = nullptr;
    last_time               = 0;
    favor_triangle          = false;
}
void StateAutocmpl::display() {
    if (core.common_button_left_pressed &&
        core.common_pn_mouse &&
        core.common_sketch_curve.size() > 3 &&
        (clock() - last_time) * 1000 / CLOCKS_PER_SEC > core.configTemp.autocmpl_preview_interval)
    {
        init();
        last_time = clock();
        MinSelector<curvenetwork::Halfedge*> h_min;
        for (auto& h : core.curvenetwork.halfedges) {
            if (h.is_deleted || h.patch() || !h.opposite->patch()) continue;
            
            double dist = pn_norm(*core.common_pn_mouse - h.vertex->pn);
            if (dist > core.configTemp.snapSize()) continue;
            h_min.update(dist, &h);
        }
        if (!h_min.value)
            return;
        
        h_selected = h_min.value;
        
        suggested_autocmpl = core.suggest_autocmpl(h_selected, false);
        selected_autocmpl_index = suggested_autocmpl.size() - 1;
    }
    
    if (selected_autocmpl_index == -1)
        return;
    
    auto& autocmpl = suggested_autocmpl[selected_autocmpl_index];
    
    // render corners with phong tessellation
    glDisable(GL_DEPTH_TEST);
    glColor4d(state_color, 0.2);
    phong_tessellation::push_config();
    phong_tessellation::enable();
    phong_tessellation::subdiv() = 5;
    phong_tessellation::weight() = 0.7;
    phong_tessellation::begin(phong_tessellation::Mode::TRIANGLE_FAN);
    for (auto& pn : autocmpl.corners)
        phong_tessellation::vertex(pn);
    phong_tessellation::end();
    phong_tessellation::pop_config();
}
void StateAutocmpl::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_autocmpl.adjust_init(mouse_x);
        return;
    }
}
void StateAutocmpl::mouse_up(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_autocmpl.adjust_done();
        return;
    }
    if (selected_autocmpl_index != -1 && core.common_pn_mouse && core.common_sketch_curve.size() == 1) {
        // accept currently selected suggestion
        core.memento_store();
        core.compute_autocompl(suggested_autocmpl[selected_autocmpl_index]);
        core.generate_patches();
        init();
    }
    core.common_sketch_curve.clear();
}
void StateAutocmpl::keyboard(unsigned char key, int x, int y) {
    if (selected_autocmpl_index == -1 || core.common_button_left_pressed)
        return;
    
    if (key == 'z' || key == 'x') {
        // select other alternatives
        selected_autocmpl_index = util::clamp<int>(selected_autocmpl_index + (key == 'z' ? -1 : 1), 0, suggested_autocmpl.size() - 1);
        
    } else if (key == 'c') {
        // toggle "favor triangle" mode
        util::flip_bool(favor_triangle);
        suggested_autocmpl = core.suggest_autocmpl(h_selected, favor_triangle);
        selected_autocmpl_index = suggested_autocmpl.size() - 1;
    }
}
