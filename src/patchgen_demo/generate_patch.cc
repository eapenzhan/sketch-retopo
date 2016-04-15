#include "decl.hh"
#include <patchgen/generate_topology.hh>
using namespace Eigen;

void demo::generate_patch(const Eigen::VectorXi& l, patchgen::PatchParam& param, Patch& patch) {
    patchgen::generate_topology(l, param, patch);
    demo::determine_geometry(patch, l);
}

void demo::generate_patch(const patchgen::PatchParam& param, Patch& patch) {
    patchgen::generate_topology(param, patch);
    demo::determine_geometry(patch, param.l);
}
