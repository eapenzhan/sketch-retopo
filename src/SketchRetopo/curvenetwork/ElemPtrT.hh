#pragma once
#include "decl.hh"

namespace curvenetwork {
    // ElemPtrT: pointer with unique index |
    //-------------------------------------+
    template <class T>
    struct ElemPtrT {
        T* ptr;
        int id;
        
        ElemPtrT() : ptr(0), id(-1) {}
        explicit ElemPtrT(T& elem) : ptr(&elem), id(elem.id) {}              // ctor with given element, used when creating new data
        explicit ElemPtrT(int id_) : ptr(0), id(id_) {}                      // ctor with element id, used when reconstructing from saved data (e.g., loading from file, undo/redo)
        
        T& operator* () const { return *ptr; }
        T* operator->() const { return  ptr; }
        operator T*() const {  return ptr; }
        operator bool() const { return ptr != 0 && !ptr->is_deleted; }
        bool operator==(const ElemPtrT<T>& rhs) const { return id == rhs.id; }
        bool operator!=(const ElemPtrT<T>& rhs) const { return !this->operator==(rhs); }
        ElemPtrT<T>& operator=(T* ptr_) { ptr = ptr_; id = ptr ? ptr->id : -1; return *this; }
    };
    
    typedef ElemPtrT<Vertex   > VertexPtr   ;
    typedef ElemPtrT<Halfedge > HalfedgePtr ;
    typedef ElemPtrT<Halfchain> HalfchainPtr;
    typedef ElemPtrT<Edgechain> EdgechainPtr;
    typedef ElemPtrT<Patch    > PatchPtr    ;
}
