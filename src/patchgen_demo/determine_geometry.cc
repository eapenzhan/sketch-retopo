#include "decl.hh"
#include <patchgen/decl.hh>
using namespace Eigen;

void demo::determine_geometry(Patch& patch, const VectorXi& l) {
    // fix boundary vertices
    Patch::HHandle h;
    for (auto v : patch.vertices()) {
        if (patch.data(v).patchgen.corner_index == 0) {
            h = patch.halfedge_handle(v);
            break;
        }
    }
    int num_sides = l.size();
    for (int i = 0; i < num_sides; ++i) {
        for (int j = 0; j < l[i]; ++j) {
            auto& vdata = patch.data(patch.from_vertex_handle(h)).laplaceDirect;
            double t = i + j / static_cast<double>(l[i]);
            vdata.value << patchgen::get_boundary_geometry(num_sides, t), 0;
            vdata.is_fixed = true;
            h = patch.prev_halfedge_handle(h);
        }
    }

    // solve
    patch.laplaceDirect_factorize();
    patch.laplaceDirect_solve();
}
