#include "../SketchRetopo.hh"
#include <kt84/MinSelector.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;

namespace {
    auto& core = SketchRetopo::get_instance();
}

StateMoveVertex::StateMoveVertex()
    : State("move vertex", Vector3d(0.8, 1, 0.2))
{}
void StateMoveVertex::init() {
    pn_prev = boost::none;
    is_smoothing = false;
}
void StateMoveVertex::mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_moveVertex.adjust_init(mouse_x);
        return;
    }
    
    if (!core.common_pn_mouse) return;
    Vector3d mouse_pos = core.common_pn_mouse->head(3);
    
    if (ctrl_pressed) {         // Harmonic parameterization
        MinSelector<curvenetwork::Patch*> patch_selector;           // Select closest patch
        for (auto& patch : core.curvenetwork.patches) {
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
            core.memento_store();
            auto hit = core.intersect(*core.common_pn_mouse);
            auto f = core.basemesh.face_handle(hit.id0);
            core.compute_patch_interior_pn_harmonic(patch_selector.value, f);
        }
    
    } else {
        core.memento_store();
        pn_prev = *core.common_pn_mouse;
        is_smoothing = shift_pressed;
    }
}
void StateMoveVertex::mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::RIGHT) {
        core.configTemp.brushSize_moveVertex.adjust_done();
        return;
    }
    
    init();
}
void StateMoveVertex::mouse_move(int mouse_x, int mouse_y) {
    if (!pn_prev || !core.common_pn_mouse)
        return;
    
    if (is_smoothing)
        core.vertices_smooth(*core.common_pn_mouse, core.configTemp.brushSize_moveVertex());
    else
        core.vertices_move(*core.common_pn_mouse, (*core.common_pn_mouse - *pn_prev).head(3), core.configTemp.brushSize_moveVertex());
    
    core.curvenetwork.invalidate_displist();
    
    pn_prev = *core.common_pn_mouse;
}
void StateMoveVertex::keyboard(unsigned char key, int x, int y) {
    if (key == 'g') {
        core.memento_store();
        core.vertices_smooth_global();
        core.curvenetwork.invalidate_displist();
    }
}
