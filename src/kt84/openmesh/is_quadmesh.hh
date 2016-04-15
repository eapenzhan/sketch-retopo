#pragma once

namespace kt84 {

template <class Mesh>
inline bool is_quadmesh(const Mesh& mesh) {
    for (auto f : mesh.faces()) {
        if (mesh.valence(f) != 4)
            return false;
    }
    return true;
}

}