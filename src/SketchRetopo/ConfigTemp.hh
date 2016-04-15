#pragma once

#include <Eigen/Core>
#include "RealParam.hh"
#include "AutoSave.hh"

struct ConfigTemp {
    // real-valued adjustable parameters
    RealParam snapSize;
    RealParam brushSize_spine;
    RealParam brushSize_autocmpl;
    RealParam brushSize_moveVertex;
    RealParam quadSize;
    RealParam segmentSize;
    
    static const double default_quadSize_ratio;
    static const double default_segmentSize_ratio;
    
    AutoSave  autoSave;
    int       autocmpl_preview_interval;
    int       cylinder_num_div;                 // default number of curves automatically generated when using cylinder sketching
    
    // threshold for loop orientation detection
    double loop_threshold_normal;
    double loop_threshold_area;
    
    // parameters for geometric snakes
    double snakes_internal1;
    double snakes_internal2;
    double snakes_external ;
    double snakes_damping  ;
    
    ConfigTemp();
};
