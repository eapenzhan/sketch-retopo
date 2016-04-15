#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Patch.hh"
#include "Circulator.hh"
#include <kt84/util.hh>
#include <kt84/eigen_util.hh>
#include <kt84/openmesh/append_quad_strip.hh>
#include <kt84/openmesh/flip_faces.hh>
#include <kt84/graphics/phong_tessellation.hh>
#include <kt84/graphics/graphics_util.hh>
#include <patchgen/decl.hh>
#include <patchgen/generate_topology.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

curvenetwork::Patch curvenetwork::Patch::imaginary_patch_symmetry(-2);
curvenetwork::Patch curvenetwork::Patch::imaginary_patch_boundary(-3);
double              curvenetwork::Patch::quadSize;
bool                curvenetwork::Patch::use_even_num_subdiv = true;
bool                curvenetwork::Patch::prefer_rect_proc3   = true;
const int           curvenetwork::Patch::failure_patch_num_subdiv = 4;
bool curvenetwork::Patch::is_ordinary(const Patch* patch) { return patch && !patch->is_failure && !patch->is_imaginary(); }

curvenetwork::Patch::Patch(int id_)
    : id        (id_)
    , is_deleted()
    , flag      ()
    , is_failure()
    , is_hidden()
{
    if (id < -1)
        // for imaginary patches, add a single imaginary point so that it is not seen as empty
        add_vertex(Point());
}

void curvenetwork::Patch::clear() {
    PatchBase::clear();
    invalidate_displist();
}
int curvenetwork::Patch::num_corners() const {
    int result = 0;
    for (auto c = halfchain; ; ) {
        if (!c)
            // bug!
            return 0;
        if (c->is_corner()) ++result;
        c = c->next();
        if (c == halfchain) break;
    }
    return result;
}
curvenetwork::Patch::HHandle curvenetwork::Patch::opposite_boundary_halfedge(HHandle boundary_halfedge) const {
    auto& h_data = data(boundary_halfedge);
    auto c = h_data.halfchain;
    auto c_opposite = c->opposite();
    
    auto p_opposite = opposite_patch(boundary_halfedge);
    if (!p_opposite)
        return HHandle();
    
    for (auto h : p_opposite->halfedges()) {
        auto& h_opposite_data = p_opposite->data(h);
        if (h_opposite_data.halfchain == c_opposite && h_opposite_data.index_wrt_halfchain == c->num_subdiv() - 1 - h_data.index_wrt_halfchain)
            return h;
    }
    
    return HHandle();
}
curvenetwork::Patch*         curvenetwork::Patch::opposite_patch            (HHandle boundary_halfedge) const {
    auto& h_data = data(boundary_halfedge);
    auto c = h_data.halfchain;
    auto c_opposite = c->opposite();
    
    if (!Patch::is_ordinary(c_opposite->patch) || c->patch == c_opposite->patch)
        return nullptr;
    
    return c_opposite->patch;
}
curvenetwork::Vertex*        curvenetwork::Patch::vertex_patch_to_curveNetwork(VHandle v) const {
    if (!is_boundary(v))
        return nullptr;
    
    auto h = prev_halfedge_handle(halfedge_handle(v));
    if (data(h).index_wrt_halfchain != 0)
        return nullptr;
    
    return data(h).halfchain->halfedge_front->from_vertex();
}
curvenetwork::Patch::VHandle curvenetwork::Patch::vertex_curveNetwork_to_patch(curvenetwork::Vertex* v) const {
    if (!v->is_corner())
        return VHandle();
    
    for (auto patch_v : vertices()) {
        if (!is_boundary(patch_v))
            continue;
        
        auto w = vertex_patch_to_curveNetwork(patch_v);
        if (w == v)
            return patch_v;
    }
    return VHandle();
}
curvenetwork::Patch::VHandle curvenetwork::Patch::patch_corner(int corner_index) const {
    for (auto v : vertices()) {
        if (data(v).patchgen.corner_index == corner_index) return v;
    }
    return VHandle();
}
int    curvenetwork::Patch::num_subdiv         (int side_index) const {
    int current_index = 0;
    int result        = 0;
    
    for (auto c = halfchain; ; ) {
        result += c->num_subdiv();
        
        if (c->is_corner()) {
            if (side_index == current_index)
                return result;
            
            ++current_index;
            result = 0;
        }
        
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    
    return 0;
}
int    curvenetwork::Patch::num_free_halfchains(int side_index) const {
    int current_index = 0;
    int result        = 0;
    
    for (auto c = halfchain; ; ) {
        if (c->num_subdiv() == 0)
            ++result;
        
        if (c->is_corner()) {
            if (current_index == side_index)
                return result;
            
            ++current_index;
            result = 0;
        }
        
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    return true;
}
double curvenetwork::Patch::side_length        (int side_index) const {
    int current_index = 0;
    double result     = 0;
    
    for (auto c = halfchain; ; ) {
        result += c->length();
        
        if (c->is_corner()) {
            if (side_index == current_index)
                return result;
            
            ++current_index;
            result = 0;
        }
        
        c = c->next();
        if (c == halfchain)
            break;  // this should not happen!
    }
    
    return 0;
}
void   curvenetwork::Patch::change_num_subdiv  (int side_index, int target_num_subdiv) {
    int delta_target = target_num_subdiv - num_subdiv(side_index);
    if (delta_target <= 0)
        return;
    
    vector<Halfchain*> c_free;
    int current_index = 0;
    for (auto c = halfchain; ; ) {
        if (current_index == side_index && c->num_subdiv() == 0)
            c_free.push_back(c);
        
        if (c->is_corner()) {
            if (current_index == side_index)
                break;
            ++current_index;
        }
        
        c = c->next();
        if (c == halfchain)
            break;      // this should not happen!
    }
    assert(!c_free.empty());
    
    double length_sum = accumulate(
        c_free.begin(), c_free.end(), 0.0,
        [] (double result, Halfchain* c) { return result + c->length(); });
    
    for_each(c_free.begin(), c_free.end(),
        [delta_target, length_sum] (Halfchain* c) {
            c->edgechain->num_subdiv = static_cast<int>(floor(0.5 + delta_target * c->length() / length_sum));
            
            if (use_even_num_subdiv && c->edgechain->num_subdiv % 2)
                // ensure even num_subdiv option
                ++c->edgechain->num_subdiv;
    });
    
    // not sure if this is really useful...
    //delta_target = target_num_subdiv - num_subdiv(side_index);
    //
    //if (delta_target != 0)
    //    c_free.front()->edgechain->num_subdiv += delta_target;
}
pair<curvenetwork::Patch::HHandle, curvenetwork::Patch::HHandle> curvenetwork::Patch::get_side_halfedges(int side_index) const {
    auto h_front = prev_halfedge_handle(halfedge_handle(patch_corner(side_index)));
    auto h = h_front;
    while (data(from_vertex_handle(h)).patchgen.corner_index == -1)
        h = prev_halfedge_handle(h);
    auto h_back = h;
    return make_pair(h_front, h_back);
}
Polyline_PointNormal  curvenetwork::Patch::get_side_polyline (int side_index) const {
    Polyline_PointNormal polyline;
    auto h_pair = get_side_halfedges(side_index);
    auto h = h_pair.first;
    polyline.push_back(data(to_vertex_handle(h)).pn);
    while (true) {
        polyline.push_back(data(from_vertex_handle(h)).pn);
        if (h == h_pair.second) break;
        h = prev_halfedge_handle(h);
    }
    return polyline;
}
double curvenetwork::Patch::distance(const Vector3d& p) const {
    double dist_min = util::dbl_max();
    for (auto f : faces()) {
        auto fv = cfv_iter(f);
        Vector3d p0 = data(*fv).pn.head(3);    ++fv;
        Vector3d p1 = data(*fv).pn.head(3);    ++fv;
        Vector3d p2 = data(*fv).pn.head(3);    ++fv;
        Vector3d p3 = data(*fv).pn.head(3);
        double dist0 = *eigen_util::distance_to_triangle(p0, p1, p2, p, true);
        double dist1 = *eigen_util::distance_to_triangle(p0, p2, p3, p, true);
        dist_min = util::min(dist_min, dist0, dist1);
    }
    return dist_min;
}
void curvenetwork::Patch::add_padding() {
    const int num_sides = param.get_num_sides();
    // IMPORTANT NOTE: boundary halfedges are oriented CLOCKWISE!
    Patch::HHandle h_boundary_front;
    for (auto v : vertices()) {
        if (data(v).patchgen.corner_index == 0) {
            h_boundary_front = halfedge_handle(v);
            break;
        }
    }
    
    for (int i = num_sides - 1; i >= 0; --i) {
        Patch::HHandle h_boundary_back = h_boundary_front;
        while (data(to_vertex_handle(h_boundary_back)).patchgen.corner_index == -1)
            h_boundary_back = next_halfedge_handle(h_boundary_back);
        
        for (int j = 0; j < param.p[i]; ++j) {
            // clear corner flag for current corner vertices
            data(from_vertex_handle(h_boundary_front)).patchgen.corner_index = -1;
            data(to_vertex_handle  (h_boundary_back )).patchgen.corner_index = -1;
            
            kt84::append_quad_strip(*this, h_boundary_front, h_boundary_back);
            
            // set corner flag for new corner vertices
            data(from_vertex_handle(h_boundary_front)).patchgen.corner_index = (i + 1) % num_sides;
            data(to_vertex_handle  (h_boundary_back )).patchgen.corner_index = i;
        }
        
        // go to next (clockwise) side
        h_boundary_front = next_halfedge_handle(h_boundary_back);
    }
    assert(data(from_vertex_handle(h_boundary_front)).patchgen.corner_index == 0);
}
void curvenetwork::Patch::generate_topology(bool init_param) {
    // set default num_subdiv for unspecified edgechains
    int num_sides = num_corners();
    
    if (init_param) {
        for (int i = 0; i < num_sides; ++i) {
            int j = (i + 2) % num_sides;
            bool is_fixed_i = num_free_halfchains(i) == 0;
            bool is_fixed_j = num_free_halfchains(j) == 0;
            if (!is_fixed_i && !is_fixed_j) {
                int target_num_subdiv = util::max<int>(
                    static_cast<int>(side_length(i) / quadSize),
                    static_cast<int>(side_length(j) / quadSize),
                    num_subdiv(i) + num_free_halfchains(i),
                    num_subdiv(j) + num_free_halfchains(j));
                change_num_subdiv(i, target_num_subdiv);
                change_num_subdiv(j, target_num_subdiv);
            } else if (!is_fixed_i) {
                change_num_subdiv(i, num_subdiv(j));
            } else if (!is_fixed_j) {
                change_num_subdiv(j, num_subdiv(i));
            }
        }
        
        VectorXi l(num_sides);
        for (int i = 0; i < num_sides; ++i)
            l[i] = num_subdiv(i);
        
        is_failure = l.sum() % 2 != 0;
        if (is_failure)
            l.fill(failure_patch_num_subdiv);
        
        // determine parameter
        param = patchgen::get_default_parameter(l);
    }
    
    // generate topology
    patchgen::generate_topology(param, *this);
    
    garbage_collection();
    
#ifndef NDEBUG
    debugInfo_get(false, false, false);
#endif
    
    laplaceSolver_init();
    
    set_halfedgeData();
    
    set_boundary_pn();
}
void curvenetwork::Patch::set_halfedgeData() {
    if (is_failure)
        return;
    
    auto v_corner = patch_corner(0);
    if (!v_corner.is_valid())
        // bug!
        return;
    auto h = prev_halfedge_handle(halfedge_handle(v_corner));
    
    auto c = halfchain;
    if (!c)
        // bug!
        return;
    
    // debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef NDEBUG
    cout << "hogehoge!=============\n";
    for (auto h2 = h; ; ) {
        cout << "halfedge: " << h2.idx() << ", to-vertex corner: " << data(to_vertex_handle(h2)).patchgen.corner_index << endl;
        h2 = prev_halfedge_handle(h2);
        if (h2 == h)
            break;
    }
    for (auto c2 = c; ; ) {
        cout << "halfchain: " << c2->id << ", num_subdiv: " << c2->num_subdiv() << endl;
        c2 = c2->next();
        if (c2 == c)
            break;
    }
    cout << "\n\n";
#endif
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
    
    while (true) {
        int index_wrt_halfchain = 0;
        
        for (int i = 0; i < c->num_subdiv(); ++i) {
            data(h).halfchain = c;
            data(h).index_wrt_halfchain = index_wrt_halfchain;
            
            h = prev_halfedge_handle(h);
            ++index_wrt_halfchain;
        }
        
        c = c->next();
        if (c == halfchain) {
            assert(data(to_vertex_handle(h)).patchgen.corner_index == 0);
            break;
        }
    }
}

void curvenetwork::Patch::laplaceSolver_init() {
    for (auto v : vertices()) {
        auto& vdata = data(v);
        
        vdata.laplaceDirect.is_fixed = vdata.laplaceIterative.is_fixed = is_boundary(v);
    }

    laplaceDirect_factorize();
}

void curvenetwork::Patch::set_boundary_pn() {
    auto v_corner = patch_corner(0);
    if (!v_corner.is_valid())
        // bug!
        return;
    
    auto h0 = prev_halfedge_handle(halfedge_handle(v_corner));
    
    if (is_failure) {
        // for patches which system failed to generate due to impossible num_subdiv configuration --> don't care about consistent connectivity, just uniformly sample along the boundary polyline
        HHandle    h = h0;
        Halfchain* c = halfchain;
        
        int num_corners_ = num_corners();
        for (int i = 0; i < num_corners_; ++i) {
            // construct polyline for this side
            Polyline_PointNormal polyline_all;
            while (true) {
                auto polyline_c = c->toPolyline();
                polyline_all.insert(polyline_all.end(), polyline_c.begin(), polyline_c.end());
                
                if (c->is_corner())
                    break;
                
                polyline_all.pop_back();
                c = c->next();
            }
            
            // uniform sampling
            for (int j = 0; j < failure_patch_num_subdiv; ++j, h = prev_halfedge_handle(h)) {
                double t = j / static_cast<double>(failure_patch_num_subdiv);
                data(to_vertex_handle(h)).pn = polyline_all.point_at(t);
            }
            
            c = c->next();
        }
        
    } else {
        for (auto h = h0; ; ) {
            auto& hdata = data(h);
            auto& pn = data(to_vertex_handle(h)).pn;
            auto c = hdata.halfchain;
            
            auto h_opposite = opposite_boundary_halfedge(h);
            if (h_opposite.is_valid()) {
                // copy from adjacent patch boundary
                auto opposite_patch = c->opposite()->patch;
                pn = opposite_patch->data(opposite_patch->from_vertex_handle(h_opposite)).pn;
            
            } else {
                // sample pn on Halfchain using arc-length parameterization
                double t = hdata.index_wrt_halfchain / static_cast<double>(c->num_subdiv());
                pn = c->pn_at(t);
            }
            
            h = prev_halfedge_handle(h);
            if (h == h0) break;
        }
    }
}

bool curvenetwork::Patch::is_imaginary_symmetry() const { return this == &imaginary_patch_symmetry; }
bool curvenetwork::Patch::is_imaginary_boundary() const { return this == &imaginary_patch_boundary; }
bool curvenetwork::Patch::is_imaginary() const { return is_imaginary_symmetry() || is_imaginary_boundary(); }

void curvenetwork::Patch::render_phong_fill() {
    displist_phong_fill.render([&] () {
        for (auto f : faces()) {
            phong_tessellation::begin(phong_tessellation::Mode::TRIANGLE_FAN);
            for (auto v : fv_range(f))
                phong_tessellation::vertex(data(v).pn);
            phong_tessellation::end();
        }
    });
}
void curvenetwork::Patch::render_phong_line() {
    displist_phong_line.render([&] () {
        // render all but boundary edges
        phong_tessellation::begin(phong_tessellation::Mode::LINES);
        for (auto e : edges()) {
            if (is_boundary(e))
                continue;
            
            auto v0 = to_vertex_handle(halfedge_handle(e, 0));
            auto v1 = to_vertex_handle(halfedge_handle(e, 1));
            phong_tessellation::vertex(data(v0).pn);
            phong_tessellation::vertex(data(v1).pn);
        }
        phong_tessellation::end();
    });
}
void curvenetwork::Patch::render_phong_line_boundary() {
    displist_phong_line_boundary.render([&] () {
        // boundary edges only
        phong_tessellation::begin(phong_tessellation::Mode::LINES);
        for (auto e : edges()) {
            if (!is_boundary(e))
                continue;
            
            auto v0 = to_vertex_handle(halfedge_handle(e, 0));
            auto v1 = to_vertex_handle(halfedge_handle(e, 1));
            phong_tessellation::vertex(data(v0).pn);
            phong_tessellation::vertex(data(v1).pn);
        }
        phong_tessellation::end();
    });
}
void curvenetwork::Patch::render_flat_fill() {
    displist_flat_fill.render([&] () {
        glBegin(GL_QUADS);
        for (auto f : faces()) {
            Vector3d normal = Vector3d::Zero();
            for (auto v : fv_range(f))
                normal += data(v).pn.tail(3);
            normal.normalize();
            glNormal3d(normal);
            for (auto v : fv_range(f))
                glVertex3dv(&data(v).pn[0]);
        }
        glEnd();
    });
}
void curvenetwork::Patch::render_flat_line() {
    displist_flat_line.render([&] () {
        glBegin(GL_LINES);
        for (auto e : edges()) {
            for (int i = 0; i < 2; ++i)
                glVertex3dv(&data(to_vertex_handle(halfedge_handle(e, i))).pn[0]);
        }
        glEnd();
    });
}
void curvenetwork::Patch::render_irregular_vertices() const {
    glBegin(GL_POINTS);
    for (auto v : vertices()) {
        int valence_ = valence(v);
        bool is_corner = data(v).patchgen.corner_index != -1;
        if (is_corner)
            valence_ += 2;
        else if (is_boundary(v))
            ++valence_;
        
        if (valence_ == 4)
            continue;
        glColor3d(valence_ == 3 ? Vector3d(0.2, 0.5, 1) : Vector3d(1, 0.5, 0.2));
        glVertex3dv(&data(v).pn[0]);
    }
    glEnd();
}
void curvenetwork::Patch::invalidate_displist() {
    displist_phong_fill         .invalidate();
    displist_phong_line         .invalidate();
    displist_phong_line_boundary.invalidate();
    displist_flat_fill          .invalidate();
    displist_flat_line          .invalidate();
}
