#pragma once
#include "State.hh"
#include <kt84/geometry/PolylineT.hh>

struct StateEdgeLoop : public State {
    kt84::Polyline_PointNormal edgeLoop_preview;
    
    StateEdgeLoop();
    void init();
    void display();
    void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void mouse_move(int mouse_x, int mouse_y);
    EditLevel get_editLevel() const { return EditLevel::QUAD_MESH; }
};
