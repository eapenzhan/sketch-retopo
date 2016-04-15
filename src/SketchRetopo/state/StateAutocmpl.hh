#pragma once

#include "State.hh"
#include "../curvenetwork/decl.hh"
#include <ctime>
#include "../AutocmplInfo.hh"

struct StateAutocmpl : public State {
    std::vector<AutocmplInfo> suggested_autocmpl;
    int                       selected_autocmpl_index;
    curvenetwork::Halfedge*   h_selected;
    clock_t                   last_time;
    bool                      favor_triangle;
    
    StateAutocmpl();
    void init();
    void display();
    void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void keyboard(unsigned char key, int x, int y);
    EditLevel get_editLevel() const { return EditLevel::CURVE_NETWORK; }
};

