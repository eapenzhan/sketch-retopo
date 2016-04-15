#pragma once

#include "Button.hh"
#include <string>
#include <Eigen/Core>

struct State {
    std::string     state_name;
    Eigen::Vector3d state_color;
    
    State() {}
    State(const std::string& state_name_, const Eigen::Vector3d& state_color_)
        : state_name (state_name_ )
        , state_color(state_color_)
    {}
    virtual ~State() = 0;
    
    // virtual functions (with empty default)
    virtual bool show_brushSize_snap() { return true; }
    virtual void init() {}                // called when switched to this state
    virtual void display() {}
    virtual void mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {}
    virtual void mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed) {}
    virtual void mouse_move(int mouse_x, int mouse_y) {}
    virtual void keyboard(unsigned char key, int x, int y) {}
    virtual void keyboardup(unsigned char key, int x, int y) {}
    enum class EditLevel {
        CURVE_NETWORK,
        QUAD_MESH
    };
    virtual EditLevel get_editLevel() const = 0;
};
