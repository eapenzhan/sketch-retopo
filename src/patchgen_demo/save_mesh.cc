#include <OpenMesh/Core/IO/MeshIO.hh>
#include "decl.hh"

bool demo::save_mesh(Patch& patch, const std::string& filename) {
    // copy laplaceDirect.value to point(v)
    for (auto v : patch.vertices()) {
        auto p = patch.data(v).laplaceDirect.value;
        patch.set_point(v, Patch::Point(p.x(), p.y(), p.z()));
    }
    return OpenMesh::IO::write_mesh(patch, filename);
}
