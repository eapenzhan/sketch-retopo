#pragma once
#include <iterator>
#include <boost/range/algorithm.hpp>
#include <boost/optional.hpp>

namespace kt84 {
namespace container_util {
    template<class Map>
    inline boost::optional<typename Map::mapped_type> at_optional(Map& map, const typename Map::key_type& key) {
        auto p = map.find(key);
        if (p == map.end())
            return boost::none;
        else
            return p->second;
    }
    template<class Container>
    inline int mod_index(const Container& c, int index) {
        while (index < 0) index += c.size();
        return index % c.size();
    }
    template<class Container>
    inline auto at_mod(Container& c, int index) -> decltype(c[index]) {
        return c[mod_index(c,index)];
    }
    template<class Container>
    inline auto at_mod(const Container& c, int index) -> decltype(c[index]) {
        return c[mod_index(c,index)];
    }
    template<class Container, class Predicate>
    inline void remove_if(Container& c, Predicate pred) {
        c.erase(boost::range::remove_if(c, pred), c.end());
    }
    template<class Container, class Type> 
    inline void bring_front(Container& c, const Type& val) {
        boost::range::rotate(c, boost::range::find(c, val));
    }
    template<class Container, class Type> 
    inline auto find(const Container& c, const Type& val) -> decltype(boost::range::find(c, val)) {
        return boost::range::find(c, val);
    }
    template<class Container> 
    inline auto erase_at(Container& c, int index) -> decltype(c.erase(c.begin())) {
        auto pos = c.begin();
        std::advance(pos, index);
        return c.erase(pos);
    }
    template<class Container>
    inline void remove_duplicate(Container& c) {
        boost::range::sort(c);
        boost::range::unique(c);
    }
}
}
