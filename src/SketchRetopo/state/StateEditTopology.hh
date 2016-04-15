#pragma once

#include "State.hh"
#include "../curvenetwork/decl.hh"

struct StateEditTopology : public State {
    curvenetwork::Patch* current_patch;
    Eigen::Vector3d      p_prev;
    double               average_edge_length;
    int selected_variable;
    
    StateEditTopology();
    void init();
    void display();
    void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed_, bool ctrl_pressed, bool alt_pressed);
    void mouse_move(int mouse_x, int mouse_y);
    void keyboard(unsigned char key, int x, int y);
    EditLevel get_editLevel() const { return EditLevel::QUAD_MESH; }
};
