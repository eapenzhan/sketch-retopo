#pragma once
#include "State.hh"
#include <boost/optional.hpp>
#include <kt84/geometry/PointNormal.hh>

struct StateMoveVertex : public State {
    boost::optional<kt84::PointNormal> pn_prev;
    bool is_smoothing;
    
    StateMoveVertex();
    bool show_brushSize_snap() { return false; }
    void init();
    void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void mouse_move(int mouse_x, int mouse_y);
    void keyboard(unsigned char key, int x, int y);
    EditLevel get_editLevel() const { return EditLevel::QUAD_MESH; }
};
