#pragma once
#include <string>

struct AutoSave {
    bool        enabled;
    std::string filename;
    int         interval;
    bool        unsaved;
    clock_t     last_time;
    
    AutoSave();
};
