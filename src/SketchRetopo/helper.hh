#pragma once

#include <string>
#include <algorithm>
#include <boost/lexical_cast.hpp>

// function to decompose fname in the format of <fname_core>.<fname_num>.xml
inline bool decompose_xml_fname(const std::string& fname, std::string& fname_core, int& fname_num) {
    auto pos1 = fname.find_last_of('.');
    if (pos1 == std::string::npos) return false;
    
    auto fname_ext = fname.substr(pos1 + 1, fname.size() - pos1 - 1);
    if (fname_ext != "xml" && fname_ext != "XML") return false;
    
    auto pos2 = fname.find_last_of('.', pos1 - 1);
    if (pos2 == std::string::npos) return false;
    
    auto num_str = fname.substr(pos2 + 1, pos1 - pos2 - 1);
    fname_num = -1;
    try { fname_num = boost::lexical_cast<unsigned int>(num_str); } catch (const boost::bad_lexical_cast&) {}
    if (fname_num == -1) return false;
    
    fname_core = fname.substr(0, pos2);
    
    return true;
};

inline std::string de_stringify (const std::string* p, size_t n) {
    return std::accumulate(p, p + n / sizeof(std::string), std::string());
}
