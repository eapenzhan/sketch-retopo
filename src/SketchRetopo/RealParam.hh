#pragma once
#include <kt84/util.hh>
#include <boost/optional.hpp>
#include <utility>

struct RealParam {
    double value;
    double value_max;
    double value_min;
    bool   is_relative;     // true if the absolute value is relative to the camera distance to the scene center
    double adjust_rate;     // rate of value adjustment w.r.t. mouse move in pixels
    double camera_distance;
    boost::optional<std::pair<int, double>> adjust_init_condition;
    
    RealParam()
        : value(0)
        , value_max( kt84::util::dbl_max())
        , value_min(-kt84::util::dbl_max())
        , is_relative()
        , adjust_rate(1)
        , camera_distance(1)
    {}
    RealParam(double value_, double value_max_, double value_min_, bool is_relative_, double adjust_rate_)
        : value    (value_    )
        , value_max(value_max_)
        , value_min(value_min_)
        , is_relative(is_relative_)
        , adjust_rate(adjust_rate_)
        , camera_distance(1)
    {}
    
    double operator()() const { return (is_relative ? camera_distance : 1.0) * value; }
    
    void set_camera_distance(double camera_distance_) { camera_distance = camera_distance_; }
    
    void adjust_init(int mouse_x) { adjust_init_condition = std::make_pair(mouse_x, value); }
    
    void adjust_move(int mouse_x) {
        int dx = mouse_x - adjust_init_condition->first;
        value = adjust_init_condition->second + adjust_rate * dx;
        value = kt84::util::clamp(value, value_min, value_max);
    }
    
    void adjust_done() { adjust_init_condition = boost::none; }
    
    bool is_adjusting() const { return static_cast<bool>(adjust_init_condition); }
};
