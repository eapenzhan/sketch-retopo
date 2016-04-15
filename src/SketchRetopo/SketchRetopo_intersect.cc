#include "SketchRetopo.hh"
#include <kt84/graphics/graphics_util.hh>
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

boost::optional<PointNormal> SketchRetopo::uv_to_pn(const Vector2d& uv) const {
    auto fbc = basemesh.uv_to_fbc(uv);
    embree::Hit hit;
    if (fbc.f.is_valid()) {
        hit.id0 = fbc.f.idx();
        hit.u = static_cast<float>(fbc.bc[1]);
        hit.v = static_cast<float>(fbc.bc[2]);
    }
    
    return intersect_convert(hit);
}
boost::optional<Vector2d> SketchRetopo::pn_to_uv(const PointNormal& pn) const {
    auto hit = intersect(pn);
    if (!hit) return boost::none;
    
    auto f = basemesh.face_handle(hit.id0);
    BaryCoord bc(1 - hit.u - hit.v, hit.u, hit.v);
    
    return basemesh.fbc_to_uv(f, bc);
}
boost::optional<PointNormal> SketchRetopo::intersect_convert(const embree::Hit& hit) const {
    if (!hit) return boost::none;
    
    auto f = basemesh.face_handle(hit.id0);
    auto v = basemesh.cfv_iter(f);
    
    PointNormal pn0 = basemesh.get_pn(*  v);
    PointNormal pn1 = basemesh.get_pn(*++v);
    PointNormal pn2 = basemesh.get_pn(*++v);
    
    PointNormal result =
        static_cast<double>(1 - hit.u - hit.v) * pn0 +
        static_cast<double>(hit.u)             * pn1 +
        static_cast<double>(hit.v)             * pn2;
    
    pn_normalize(result);
    return result;
}
embree::Hit SketchRetopo::intersect(const PointNormal& pn_in) const {
    if (!embree_intersector) return embree::Hit();
    
    Vector3d p = pn_in.head(3);
    Vector3d n = pn_in.tail(3);
    
    auto make_Ray = [] (const Vector3d& org, const Vector3d& dir) {
        auto org_f = org.cast<float>();
        auto dir_f = dir.cast<float>();
        return embree::Ray(
            embree::Vec3f(org_f[0], org_f[1], org_f[2]),
            embree::Vec3f(dir_f[0], dir_f[1], dir_f[2]));
    };
    
    auto delta = basemesh.boundingBox_diagonal_norm() * 0.0001 * n;
    
    embree::Hit hit[4];
    embree_intersector->intersect(make_Ray(p - delta,  n), hit[0]);
    embree_intersector->intersect(make_Ray(p + delta,  n), hit[1]);
    embree_intersector->intersect(make_Ray(p - delta, -n), hit[2]);
    embree_intersector->intersect(make_Ray(p + delta, -n), hit[3]);
    
    embree::Hit hit_min;
    for (int i = 0; i < 4; ++i) {
        if (!hit[i])
            continue;
        
        auto face_normal = o2e(basemesh.normal(basemesh.face_handle(hit[i].id0)));
        if (face_normal.dot(n) < 0)
            continue;
        
        if (abs(hit[i].t) < abs(hit_min.t))
            hit_min = hit[i];
    }
    
    return hit_min;
}
embree::Hit SketchRetopo::intersect(int mouse_x, int mouse_y) const {
    auto modelview_matrix  = get_modelview_matrix ();
    auto projection_matrix = get_projection_matrix();
    auto viewport          = get_viewport         ();
    
    auto mouse_z = read_depth(mouse_x, mouse_y);
    
    auto p = unproject(Vector3d(mouse_x, mouse_y, mouse_z), modelview_matrix, projection_matrix, viewport);
    auto q = unproject(Vector3d(mouse_x, mouse_y, 1      ), modelview_matrix, projection_matrix, viewport);
    
    Vector3d n = p - q;
    if (n.isZero())
        return embree::Hit();
    n.normalize();
    
    PointNormal pn;
    pn << p, n;
    
    return intersect(pn);
}
void SketchRetopo::project(PointNormal& pn) const {
    auto pn_intersect = intersect_convert(intersect(pn));
    if (!pn_intersect)
        return;
    pn = *pn_intersect;
    pn.head(3) += configSaved.projOffset * pn.tail(3);
}
void SketchRetopo::reproject_all() {
    for (auto& v : curvenetwork.vertices)
        project(v.pn);
    
    for (auto& p : curvenetwork.patches) {
        for (auto v : p.vertices())
            project(p.data(v).pn);
        
        p.invalidate_displist();
    }
}
