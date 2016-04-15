#pragma once
#include <string>
#include <patchgen/PatchParam.hh>
#include "Patch.hh"

namespace demo {
    void determine_geometry(Patch& patch, const Eigen::VectorXi& l);
    
    void generate_patch(const Eigen::VectorXi& l, patchgen::PatchParam& param, Patch& patch);
    void generate_patch(const patchgen::PatchParam& param, Patch& patch);
    
    std::string get_fname(const patchgen::PatchParam& param);
    bool save_mesh(Patch& patch, const std::string& filename);
}
