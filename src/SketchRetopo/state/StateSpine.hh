#pragma once

#include "State.hh"

struct StateSpine : public State {
    StateSpine();
    void init();
    void display();
    void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    void mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    EditLevel get_editLevel() const { return EditLevel::CURVE_NETWORK; }
};
