#pragma once
#include "decl.hh"
#include <iterator>
#include <type_traits>

namespace curvenetwork {
    //-- outgoing halfedges around vertex
    template <bool IsConst>
    struct VOHIterT : public std::iterator<std::bidirectional_iterator_tag, Halfedge*> {
        typedef typename std::conditional<IsConst, const Vertex  , Vertex  >::type VertexType;
        typedef typename std::conditional<IsConst, const Halfedge, Halfedge>::type HalfedgeType;
        
        VertexType* base;
        HalfedgeType* current;
        
        VOHIterT(VertexType* base_)
            : base(base_)
            , current(base->halfedge)
        {}
        
        HalfedgeType& operator* () const { return *current; }
        HalfedgeType* operator->() const { return  current; }
        operator bool() const { return current != 0; }
        VOHIterT<IsConst>& operator++() {
            if (current) {
                current = current->opposite->next;
                if (current == base->halfedge)
                    current = 0;
            }
            return *this;
        }
        VOHIterT<IsConst>& operator--() {
            if (current) {
                current = current->opposite->prev;
                if (current == base->halfedge)
                    current = 0;
            }
            return *this;
        }
        bool operator==(const VOHIterT<IsConst>& rhs) const {
            assert(base == rhs.base);
            return current == rhs.current;
        }
    };
    //-- incoming halfedges around vertex
    template <bool IsConst>
    struct VIHIterT : public VOHIterT<IsConst> {
        VIHIterT(typename VOHIterT<IsConst>::VertexType* v) : VOHIterT<IsConst>(v) {}
        typename VOHIterT<IsConst>::HalfedgeType& operator* () const { return *VOHIterT<IsConst>::current->opposite; }
        typename VOHIterT<IsConst>::HalfedgeType* operator->() const { return  VOHIterT<IsConst>::current->opposite; }
    };
    //-- incoming halfedges around vertex
    template <bool IsConst>
    struct VOCIterT : public VOHIterT<IsConst> {
        typedef typename std::conditional<IsConst, const Halfchain, Halfchain>::type HalfchainType;
        VOCIterT(typename VOHIterT<IsConst>::VertexType* v) : VOHIterT<IsConst>(v) {}
        HalfchainType& operator* () const { return *VOHIterT<IsConst>::current->halfchain; }
        HalfchainType* operator->() const { return  VOHIterT<IsConst>::current->halfchain; }
    };
    //-- halfchains around patch
    template <bool IsConst>
    struct PCIterT : public std::iterator<std::bidirectional_iterator_tag, Halfchain*> {
        typedef typename std::conditional<IsConst, const Patch    , Patch    >::type PatchType;
        typedef typename std::conditional<IsConst, const Halfchain, Halfchain>::type HalfchainType;
        
        PatchType* base;
        HalfchainType* current;
        
        PCIterT(PatchType* base_)
            : base(base_)
            , current(base->halfchain)
        {}
        
        HalfchainType& operator* () const { return *current; }
        HalfchainType* operator->() const { return  current; }
        operator bool() const { return current != 0; }
        PCIterT<IsConst>& operator++() {
            if (current) {
                current = current->next();
                if (current == base->halfchain)
                    current = 0;
            }
            return *this;
        }
        PCIterT<IsConst>& operator--() {
            if (current) {
                current = current->prev();
                if (current == base->halfchain)
                    current = 0;
            }
            return *this;
        }
        bool operator==(const PCIterT<IsConst>& rhs) const {
            assert(base == rhs.base);
            return current == rhs.current;
        }
    };
    //-- halfedges consisting halfchain
    template <bool IsConst>
    struct CHIterT : public std::iterator<std::bidirectional_iterator_tag, Halfedge*> {
        typedef typename std::conditional<IsConst, const Halfchain, Halfchain>::type HalfchainType;
        typedef typename std::conditional<IsConst, const Halfedge , Halfedge >::type HalfedgeType;
        
        HalfchainType* base;
        HalfedgeType* current;
        
        CHIterT(HalfchainType* base_)
            : base(base_)
            , current(base->halfedge_front)
        {}

        HalfedgeType& operator* () const { return *current; }
        HalfedgeType* operator->() const { return  current; }
        operator bool() const { return current != 0; }
        CHIterT<IsConst>& operator++() {
            if (current == base->halfedge_back)
                current = 0;
            else if (current)
                current = current->next;
            return *this;
        }
        bool operator==(const CHIterT<IsConst>& rhs) const {
            assert(base == rhs.base);
            return current == rhs.current;
        }
    };
    
    // typedefs
    typedef VOHIterT<false> VOHIter;
    typedef VIHIterT<false> VIHIter;
    typedef VOCIterT<false> VOCIter;
    typedef PCIterT <false> PCIter ;
    typedef CHIterT <false> CHIter ;
    typedef VOHIterT<true > ConstVOHIter;      // const version
    typedef VIHIterT<true > ConstVIHIter;
    typedef VOCIterT<true > ConstVOCIter;
    typedef PCIterT <true > ConstPCIter ;
    typedef CHIterT <true > ConstCHIter ;
}
