#include "decl.hh"
#include <patchgen/decl.hh>
using namespace std;

string demo::get_fname(const patchgen::PatchParam& param) {
    stringstream fname;
    int num_sides = param.get_num_sides();
    fname << "out_" << num_sides << "sided";
    fname << "_pattern=" << param.pattern_id;
    fname << "_l=(";
    for (int i = 0; i < num_sides; ++i)
        fname << (i > 0 ? "," : "") << param.l[i];
    fname << ")";
    fname << "_perm=" << param.permutation.id;
    string param_str = patchgen::get_param_str(num_sides, param.pattern_id, param);
    if (!param_str.empty()) fname << "_" << param_str;
    fname << ".obj";
    return fname.str();
}
