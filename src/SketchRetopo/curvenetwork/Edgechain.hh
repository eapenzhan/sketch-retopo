#pragma once
#include "ElemPtrT.hh"

namespace curvenetwork {
    struct Edgechain {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;
        
        // connectivity info
        HalfchainPtr halfchain[2];
        
        // additional data
        int num_subdiv;
        
        explicit Edgechain(int id_ = -1);
        bool is_boundary() const;
        bool on_symmetry_plane() const;
    };
}
