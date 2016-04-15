#include "SketchRetopo.hh"
#include "curvenetwork/Circulator.hh"
#include <kt84/glut_util.hh>
#include <kt84/zenity_util.hh>
#include <kt84/container_util.hh>
#include <kt84/MinSelector.hh>
#include "helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;

bool SketchRetopo::common_mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::LEFT  ) common_button_left_pressed   = true;
    if (button == Button::RIGHT ) common_button_right_pressed  = true;
    if (button == Button::MIDDLE) common_button_middle_pressed = true;
    common_shift_pressed = shift_pressed;
    common_ctrl_pressed  = ctrl_pressed;
    common_alt_pressed   = alt_pressed;
    
    if (alt_pressed) {
        // changing view
        camera->mouse_down(mouse_x, mouse_y, 
            button == Button::LEFT   ? (ctrl_pressed ? Camera::DragMode::PAN : Camera::DragMode::ROTATE) :
            button == Button::MIDDLE ? Camera::DragMode::PAN    : Camera::DragMode::ZOOM);
        return true;
    }
    if (button == Button::RIGHT && shift_pressed && ctrl_pressed ) {
        configTemp.segmentSize.adjust_init(mouse_x);
        return true;
    }
    if (button == Button::RIGHT && shift_pressed) {
        configTemp.snapSize.adjust_init(mouse_x);
        return true;
    }
    if (button == Button::RIGHT && ctrl_pressed) {
        configTemp.quadSize.adjust_init(mouse_x);
        return true;
    }
    
    // hide mesh by drawing stroke
    if (button == Button::LEFT && common_key_pressed['h']) {
        hide_stroke.clear();
        hide_stroke.push_back({mouse_x, mouse_y});
        return true;
    }
    
    if (button == Button::LEFT && shift_pressed && ctrl_pressed && state->get_editLevel() == State::EditLevel::CURVE_NETWORK) {
        common_dragMode = CommonDragMode::Snakes;
        memento_store();
        return true;
    }
    if (button == Button::LEFT && shift_pressed                 && state->get_editLevel() == State::EditLevel::CURVE_NETWORK) {
        common_dragMode = CommonDragMode::Smooth;
        memento_store();
        return true;
    }
    
    if (button == Button::LEFT && common_pn_mouse &&
        (state == &stateSketch     ||
         state == &stateSpine      ||
         state == &stateAutocmpl   ||
         state == &stateCylinder   ||
         state == &stateDeformCurve))
    {
        common_sketch_curve.push_back(*common_pn_mouse);
    }
    
    return false;
}
bool SketchRetopo::common_mouse_up(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {
    if (button == Button::LEFT  ) common_button_left_pressed   = false;
    if (button == Button::RIGHT ) common_button_right_pressed  = false;
    if (button == Button::MIDDLE) common_button_middle_pressed = false;
    common_shift_pressed = shift_pressed;
    common_ctrl_pressed  = ctrl_pressed;
    common_alt_pressed   = alt_pressed;
    
    if (camera->drag_mode != Camera::DragMode::NONE) {
        // changing view
        
        if (camera->drag_mode == Camera::DragMode::ROTATE && shift_pressed)
            // snap to closest canonical view
            camera->snap_to_canonical();
        
        if (camera->drag_mode == Camera::DragMode::PAN && configRender.auto_camera_center) {
            // set camera center to surface point where screen center is projected onto (as in Pixologic Sculptris)
            auto hit = intersect(camera->width / 2, camera->height / 2);
            if (hit) camera->update_center(intersect_convert(hit)->head(3));
        }
        camera->mouse_up();
        return true;
    }
    
    if (configTemp.snapSize   .is_adjusting()) { configTemp.snapSize   .adjust_done(); return true; }
    if (configTemp.quadSize   .is_adjusting()) { configTemp.quadSize   .adjust_done(); curvenetwork::Patch::quadSize = configTemp.quadSize(); return true; }
    if (configTemp.segmentSize.is_adjusting()) { configTemp.segmentSize.adjust_done(); return true; }
    
    // hide mesh by drawing stroke
    if (!hide_stroke.empty()) {
        hide_mesh_by_stroke();
        hide_stroke.clear();
        return true;
    }
    
    if (common_dragMode == CommonDragMode::Smooth || common_dragMode == CommonDragMode::Snakes) {
        // collect affected patches
        vector<curvenetwork::Patch*> affected_patches;
        container_util::remove_duplicate(common_affected_edgechains);
        for (auto e : common_affected_edgechains) {
            for (int i = 0; i < 2; ++i) {
                auto patch = e->halfchain[i]->patch;
                if (!patch || patch->is_imaginary())
                    continue;
                affected_patches.push_back(patch);
            }
        }
        
        // update patch vertex positions
        container_util::remove_duplicate(affected_patches);
        for (auto patch : affected_patches)
            if (!patch->is_deleted) patch->clear();         // NOTE: patch parameter is kept intact!
        for (auto patch : affected_patches) {
            if (patch->is_deleted) continue;
            patch->generate_topology(false);
            compute_patch_interior_pn(patch);
        }
        curvenetwork.invalidate_displist();
        
        common_affected_edgechains.clear();
        common_dragMode = CommonDragMode::None;
        return true;
    }
    
    if (!common_sketch_curve.empty()) {
        // loop check, resampling, default smoothing, projection
        common_sketch_curve.is_loop = pn_norm(common_sketch_curve.back() - common_sketch_curve.front()) < configTemp.snapSize();
        common_sketch_curve.resample_by_length(configTemp.segmentSize());
        common_sketch_curve.smooth(2, 0.5, 0.5);
        project(common_sketch_curve);
    }
    
    return false;
}
bool SketchRetopo::common_mouse_move(int mouse_x, int mouse_y) {
    common_mouse_pos_prev = common_mouse_pos;
    common_mouse_pos << mouse_x, mouse_y;
    
    if (camera->drag_mode != Camera::DragMode::NONE) {
        // changing view
        camera->mouse_move(mouse_x, mouse_y);
        return true;
    }
    
    // changing real-valued parameters
    if (configTemp.snapSize.is_adjusting()) {
        configTemp.snapSize.adjust_move(mouse_x);
        // snap size shouldn't be larger than brush size.
        if (state == &stateSpine     ) configTemp.snapSize.value = min<double>(configTemp.snapSize.value, configTemp.brushSize_spine     .value * 0.9);
        if (state == &stateAutocmpl  ) configTemp.snapSize.value = min<double>(configTemp.snapSize.value, configTemp.brushSize_autocmpl  .value * 0.9);
        if (state == &stateMoveVertex) configTemp.snapSize.value = min<double>(configTemp.snapSize.value, configTemp.brushSize_moveVertex.value * 0.9);
        return true;
    }
    if (configTemp.brushSize_spine.is_adjusting()) {
        configTemp.brushSize_spine.adjust_move(mouse_x);
        // brush size shouldn't be smaller than snap size.
        configTemp.brushSize_spine.value = max<double>(configTemp.brushSize_spine.value, configTemp.snapSize.value * 0.9);
        return true;
    }
    if (configTemp.brushSize_autocmpl.is_adjusting()) {
        configTemp.brushSize_autocmpl.adjust_move(mouse_x);
        configTemp.brushSize_autocmpl.value = max<double>(configTemp.brushSize_autocmpl.value, configTemp.snapSize.value * 0.9);
        return true;
    }
    if (configTemp.brushSize_moveVertex.is_adjusting()) {
        configTemp.brushSize_moveVertex.adjust_move(mouse_x);
        configTemp.brushSize_moveVertex.value = max<double>(configTemp.brushSize_moveVertex.value, configTemp.snapSize.value * 0.9);
        return true;
    }
    if (configTemp.quadSize   .is_adjusting()) { configTemp.quadSize   .adjust_move(mouse_x); return true; }
    if (configTemp.segmentSize.is_adjusting()) { configTemp.segmentSize.adjust_move(mouse_x); return true; }
    
    // hide mesh by drawing stroke
    if (!hide_stroke.empty()) {
        hide_stroke.push_back({mouse_x, mouse_y});
        return true;
    }
    
    // compute intersection
    common_pn_mouse = intersect_convert(intersect(mouse_x, mouse_y));
    if (!common_pn_mouse)
        return true;
    
    project(*common_pn_mouse);
    
    // common curve editing tool: smoothing
    if (common_dragMode == CommonDragMode::Smooth) {
        // look for closest non-corner vertex
        MinSelector<curvenetwork::Vertex*> v_closest;
        for (auto& v : curvenetwork.vertices) {
            if (v.is_corner() || v.is_openend())
                continue;
            
            double dist = pn_norm(v.pn - *common_pn_mouse);
            v_closest.update(dist, &v);
        }
        
        if (v_closest.score < configTemp.snapSize()) {
            // perform smoothing on the corresponding edgechain
            auto c = v_closest.value->halfedge->halfchain;
            auto polyline = c->toPolyline();
            bool is_loop = polyline.is_loop;
            int  n       = polyline.size();
            
            if (!is_loop && n < 5)
                return true;
            
            polyline.smooth(25, configTemp.snakes_internal1, configTemp.snakes_internal2, configTemp.snakes_damping);
            auto h = c->halfedge_front;
            for (int i = 0; i < n; ++i) {
                auto& pn = h->from_vertex()->pn;
                pn = polyline[i];
                project(pn);
                
                h = h->next;
                if (!h)
                    break;
            }
            
            common_affected_edgechains.push_back(c->edgechain);
        }
        return true;
    }
    
    // common curve editing tool: smoothing with snakes
    if (common_dragMode == CommonDragMode::Snakes) {
        // look for closest corner/non-corner vertex
        MinSelector<curvenetwork::Vertex*> v_corner;
        MinSelector<curvenetwork::Vertex*> v_noncorner;
        for (auto& v : curvenetwork.vertices) {
            if (configSaved.symmetric && v.on_symmetry_plane())
                continue;
            
            double dist = pn_norm(v.pn - *common_pn_mouse);
            if (configTemp.snapSize() < dist)
                continue;
            (v.is_corner() || v.is_openend() ? v_corner : v_noncorner).update(dist, &v);
        }
        
        if (v_corner.value) {
            for (int i = 0; i < 5; ++i)
                snakes_move(v_corner.value);
            for (curvenetwork::VOCIter c(v_corner.value); c; ++c)
                common_affected_edgechains.push_back(c->edgechain);
            curvenetwork.invalidate_displist();
            
        } else if (v_noncorner.value) {
            auto e = v_noncorner.value->halfedge->halfchain->edgechain;
            for (int i = 0; i < 5; ++i)
                snakes_move(e);
            common_affected_edgechains.push_back(e);
        }
        return true;
    }
    
    if (common_pn_mouse && !common_sketch_curve.empty())
        common_sketch_curve.push_back(*common_pn_mouse);
    
    return false;
}
namespace {
void turntable_idle() {
    auto& core = SketchRetopo::get_instance();
    if (core.camera != &core.camera_upright) {
        auto eye    = core.camera->get_eye();
        auto center = core.camera->center;
        auto up     = core.camera->get_up();
        core.camera = &core.camera_upright;
        core.camera->init(eye, center, up);
        core.camera->update_center(center);
    }
    core.camera_upright.theta += core.configRender.turntable_speed;
    glutPostRedisplay();
}
}

bool SketchRetopo::common_keyboard(unsigned char key, int x, int y) {
    string fname;
    switch (key) {
        case 'q':
        case 'w':
        case 'e':
        case 'r':
        case 't':
        case 'y':
        case 'f':
        case 'a':
        case 's':
        case 'd':
            // mode change
            state_set(
                key == 'q' ? EnumState::Sketch       :
                key == 'w' ? EnumState::Spine        :
                key == 'e' ? EnumState::Autocmpl     :
                key == 'r' ? EnumState::Laser        :
                key == 't' ? EnumState::Cylinder     :
                key == 'y' ? EnumState::DeformCurve  :
                key == 'f' ? EnumState::EditCorner   :
                key == 'a' ? EnumState::EditTopology :
                key == 's' ? EnumState::EdgeLoop     :
                             EnumState::MoveVertex   );
            return true;
        
        case 25:        // CTRL + y
            memento_redo();
            return true;
        
        case 26:        // CTRL + z
            memento_undo();
            return true;
        
        case 6:         // CTRL + f -> toggle fullscreen
            {
                static int prev_window_x = 0;
                static int prev_window_y = 0;
                static int prev_window_width  = 0;
                static int prev_window_height = 0;
                static bool is_fullscreen = false;
                if (is_fullscreen) {
                    glutReshapeWindow (prev_window_width, prev_window_height);
                    glutPositionWindow(prev_window_x    , prev_window_y     );
                } else {
                    prev_window_x      = glutGet(GLUT_WINDOW_X     );
                    prev_window_y      = glutGet(GLUT_WINDOW_Y     );
                    prev_window_width  = glutGet(GLUT_WINDOW_WIDTH );
                    prev_window_height = glutGet(GLUT_WINDOW_HEIGHT);
                    glutFullScreen();
                }
                is_fullscreen = !is_fullscreen;
            }
            return true;
        
        case ' ':
            configRender.mode = configRender.mode == ConfigRender::Mode::QUADMESH_ONLY ? ConfigRender::Mode::DEFAULT : ConfigRender::Mode::QUADMESH_ONLY;
            return true;
        
        case 2:         // ctrl + b
            configRender.mode = configRender.mode == ConfigRender::Mode::BASEMESH_ONLY ? ConfigRender::Mode::DEFAULT : ConfigRender::Mode::BASEMESH_ONLY;
            return true;

        case 15:        // ctrl + o -> load xml
            if ((!configTemp.autoSave.unsaved || zenity_util::question("SketchRetopo", "Current data is not yet saved. Continue?")) &&   // confirm discarding unsaved data
                zenity_util::file_selection_load(fname, "SketchRetopo - load xml", "", "xml files|*.xml"))
            {
                xml_load(fname);
            }
            return true;
        
        case 19:        // ctrl + s -> save xml
            if (!configTemp.autoSave.unsaved)
                // don't save if there's no change
                return true;
            
            if (configTemp.autoSave.filename.empty()) {
                if (zenity_util::file_selection_save(fname, "SketchRetopo - save xml", "", "xml files|*.xml"))
                    xml_save(fname);
                
            } else {
                configTemp.autoSave.unsaved = false;
                configTemp.autoSave.last_time = clock();
                
                // update ausoSave.filename if it has a specific form of <fname_core>.<fname_num>.xml
                string fname = configTemp.autoSave.filename;
                string fname_core;
                int fname_num;
                if (decompose_xml_fname(fname, fname_core, fname_num)) {
                    stringstream ss;
                    ss << fname_core << '.' << (fname_num + 1) << ".xml";
                    fname = ss.str();
                }
                
                xml_save(fname);
            }
            return true;
        
        case 9:         // ctrl + i -> import obj
            if ((!configTemp.autoSave.unsaved || zenity_util::question("SketchRetopo", "Current data is not yet saved. Continue?")) &&   // confirm discarding unsaved data
                zenity_util::file_selection_load(fname, "SketchRetopo - import base mesh", "", "obj files|*.obj"))
            {
                import_basemesh(fname);
            }
            return true;
        
        case 5:         // ctrl + e -> export obj
            if (zenity_util::file_selection_save(fname, "SketchRetopo - export retopo mesh", "", "obj files|*.obj"))
                export_retopomesh(fname);
            return true;
        case '0':       // reset camera center
            camera->center = basemesh.boundingBox.center();
            return true;
        case '9':
            generate_patch_debug();
            return true;
        case '8':
            // easy camera save/load
            {
                stringstream ss;
                auto eye    = camera->get_eye();
                auto center = camera->center;
                auto up     = camera->get_up();
                int window_x = glutGet(GLUT_WINDOW_X);
                int window_y = glutGet(GLUT_WINDOW_Y);
                int window_width = glutGet(GLUT_WINDOW_WIDTH);
                int window_height = glutGet(GLUT_WINDOW_HEIGHT);
                ss <<
                    eye   [0] << " " <<
                    eye   [1] << " " <<
                    eye   [2] << " " <<
                    center[0] << " " <<
                    center[1] << " " <<
                    center[2] << " " <<
                    up    [0] << " " <<
                    up    [1] << " " <<
                    up    [2] << " " <<
                    window_x     << " " <<
                    window_y     << " " <<
                    window_width << " " <<
                    window_height;
                auto str = ss.str();
                if (!zenity_util::entry(str, "SketchRetopo", "camera parameter", str))
                    return true;
                ss.str(str);
                ss >>
                    eye   [0] >>
                    eye   [1] >>
                    eye   [2] >>
                    center[0] >>
                    center[1] >>
                    center[2] >>
                    up    [0] >>
                    up    [1] >>
                    up    [2] >>
                    window_x      >>
                    window_y      >>
                    window_width  >>
                    window_height;
                camera->init(eye, center, up);
                glutPositionWindow(window_x, window_y);
                glutReshapeWindow(window_width, window_height);
            }
            return true;
        case '7':
            // stats report
            cout << "# tri: " << basemesh.n_faces() << endl;
            {
                int n_quad = 0;
                int n_patch_tri = 0;
                int n_patch_rect = 0;
                int n_patch_pent = 0;
                for (auto& p : curvenetwork.patches) {
                    n_quad += p.n_faces();
                    int num_corners = p.num_corners();
                    ++(num_corners == 3 ? n_patch_tri : num_corners == 4 ? n_patch_rect : n_patch_pent);
                }
                cout << "# quad: " << n_quad << endl;
                cout << "# patch: " << n_patch_tri << "/" << n_patch_rect << "/" << n_patch_pent << endl;
            }
            return true;
        case '6':
            // turntable mode
            util::flip_bool(configRender.turntable_mode);
            glutIdleFunc(configRender.turntable_mode ? turntable_idle : 0);
            return true;
    }
    
    common_key_pressed[key] = true;
    
    return false;
}
bool SketchRetopo::common_keyboardup(unsigned char key, int x, int y) {
    common_key_pressed[key] = false;
    
    return false;
}

