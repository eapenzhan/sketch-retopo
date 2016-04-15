#include "SketchRetopo.hh"
#include <algorithm>
#include <unordered_map>
#include <Eigen/Sparse>
#include <kt84/container_util.hh>
#include <boost/range/algorithm.hpp>
#include <kt84/eigen_util.hh>
#include <kt84/loop_util.hh>
#include <kt84/zenity_util.hh>
#include <kt84/MinSelector.hh>
#include <kt84/adjacent_pairs.hh>
#include <kt84/math/RBFKernel.hh>
#include <kt84/openmesh/edgeloop.hh>
#include <patchgen/decl.hh>
#include "curvenetwork/Circulator.hh"
#include "curvenetwork/helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;

namespace {
    inline float bc(const embree::Hit& hit, int index) {
        return
            index == 0 ? 1 - hit.u - hit.v :
            index == 1 ? hit.u :
            index == 2 ? hit.v : numeric_limits<float>::quiet_NaN();
    }
}

bool SketchRetopo::is_loop_ccw(const Polyline_PointNormal& loop) const {
    assert(loop.is_loop);
    
    int n = loop.size();
    
    // average
    PointNormal avg = accumulate(loop.begin(), loop.end(), PointNormal::Zero().eval()) / n;
    
    // reject if normal directions vary widely
    if (avg.tail(3).norm() < configTemp.loop_threshold_normal)
        return false;
    
    // rough estimate of imagined surface area surrounded by this loop
    double area_sum = 0;
    for (int i = 0; i < loop.size(); ++i) {
        Vector3d d0 = (loop[ i         ] - avg).head(3);
        Vector3d d1 = (loop[(i + 1) % n] - avg).head(3);
        Vector3d n  = loop[i].tail(3);
        area_sum += 0.5 * d0.cross(d1).dot(n);
    }
    
    double loop_length = loop.length();
    double area_max = loop_length * loop_length / (4 * util::pi());      // maximum (i.e. disc) area for the loop length
    
    return area_max * configTemp.loop_threshold_area < area_sum;
}

void SketchRetopo::generate_patches() {
    curvenetwork.invalidate_displist();
    
    // visited flag
    curvenetwork.set_flag_halfchains(curvenetwork::Core::FlagOp::SUBST, 0);
    
    for (auto& c : curvenetwork.halfchains) {
        if (c.is_deleted || c.flag || (c.patch && !c.patch->is_deleted)) continue;
        
        if (c.halfedge_front->imaginary_patch) {
            // this halfchain is adjacent to symmetry line or basemesh boundary --> set reference to imaginary patch
            c.patch = c.halfedge_front->imaginary_patch;
            continue;
        }
        
        curvenetwork::Halfchain* c_front;
        curvenetwork::Halfchain* c_back;
        curvenetwork::trace_halfchains(&c, c_front, c_back);
        
        // set visited flag, count corners
        int num_corners = 0;
        for (auto c1 = c_front; ; c1 = c1->next()) {
            c1->flag = 1;
            if (c1->is_corner()) ++num_corners;
            if (c1 == c_back)
                break;
        }
        
        if (c_front->prev() != c_back || (num_corners < 2 || 6 < num_corners))
            // skip if c doesn't form {2, 3, 4, 5, 6}-sided closed polygon
            continue;
        
        // determine if this loop is counterclockwise
        Polyline_PointNormal polyline;
        for (auto c1 = c_front; ; c1 = c1->next()) {
            auto polyline_sub = c1->toPolyline();
            polyline_sub.pop_back();
            polyline.insert(polyline.end(), polyline_sub.begin(), polyline_sub.end());
            if (c1 == c_back)
                break;
        }
        polyline.is_loop = true;
        
        double loop_length = polyline.length();
        if (loop_length > basemesh.boundingBox_diagonal_norm() * 2)
            // loop is too long
            continue;
        
        if (!is_loop_ccw(polyline))
            // doesn't form a reasonable counterclockwise loop --> skip
            continue;

        // instantiate new patch
        auto patch = curvenetwork.new_patch();
        
        // set connectivity info
        for (auto c1 = c_front; ; ) {
            c1->patch = patch;
            c1 = c1->next();
            if (c1 == c_front) break;
        }
        
        // ensure that patch->halfchain->prev is corner
        patch->halfchain = c_front;
        while(!patch->halfchain->is_corner())
            patch->halfchain = patch->halfchain->next();
        patch->halfchain = patch->halfchain->next();
        
        // generate topology and set uv in abstract parameter space
        patch->generate_topology();
        
        // compute pn info for patch interior vertices (maybe using local harmonic parameterization?)
        compute_patch_interior_pn(patch);
    }
}
void SketchRetopo::generate_patch_debug() {
    if (curvenetwork.patches.size() != 1)
        // only do testing when there is a single patch
        return;
    
    auto& patch = curvenetwork.patches.front();
    
    int sum = 0;
    for (curvenetwork::PCIter c(&patch); c; ++c) {
        int num_subdiv = rand() % 20 + 1;
        c->edgechain->num_subdiv = num_subdiv;
        sum += num_subdiv;
    }
    if (sum % 2)
        ++patch.halfchain->edgechain->num_subdiv;
    
    patch.generate_topology();
    compute_patch_interior_pn(&patch);
}

void SketchRetopo::add_curve_open_delegate(Polyline_PointNormal& curve, bool corner_at_openend) {
    // helper function to get offset to nearby curvenetwork::Vertex
    auto find_snap = [this] (const PointNormal& pn, PointNormal& result) -> bool {
        MinSelector<curvenetwork::Vertex*> v_openend;
        MinSelector<curvenetwork::Vertex*> v_corner ;
        MinSelector<curvenetwork::Vertex*> v_side   ;
        
        for (auto& v : curvenetwork.vertices) {
            if (v.is_deleted)
                continue;
            
            bool is_openend = v.is_openend();
            bool is_corner  = v.is_corner ();
            
            double dist = pn_norm(pn - v.pn);
            if (configTemp.snapSize() < dist)
                continue;
            
            if (is_openend)
                v_openend.update(dist, &v);
            else if (is_corner)
                v_corner.update(dist, &v);
            else
                v_side.update(dist, &v);
        }
        
        if (v_openend.value) { result = v_openend.value->pn; return true; }
        if (v_corner .value) { result = v_corner .value->pn; return true; }
        if (v_side   .value) { result = v_side   .value->pn; return true; }
        
        return false;
    };
    
    // offset endpoints based on snapping
    PointNormal pn_snap;
    bool is_snapped = false;
    if (find_snap(curve.front(), pn_snap)) {
        curve.offset_front(pn_snap);
        is_snapped = true;
    }
    if (find_snap(curve.back(), pn_snap)) {
        curve.offset_back(pn_snap);
        is_snapped = true;
    }
    
    if (is_snapped)
        // project again
        project(curve);
    
    // speed-dependent snap radius
    double snapSize = configTemp.snapSize();
    double delta_max = 10;
    double mouse_delta = min<double>((common_mouse_pos - common_mouse_pos_prev).cast<double>().norm(), delta_max);
    snapSize *= 1 + mouse_delta / delta_max;
    
    curvenetwork.add_curve_open(curve, configTemp.snapSize(), corner_at_openend);
}

void SketchRetopo::compute_patch_interior_pn(curvenetwork::Patch* patch) {
    // initial solution by direct Laplace solve
    int n_interior_vertices = 0;
    for (auto v : patch->vertices()) {
        if (!patch->is_boundary(v)) {
            ++n_interior_vertices;
            continue;
        }
        
        auto& vdata = patch->data(v);
        vdata.laplaceDirect.value = vdata.pn;
    }
    patch->laplaceDirect_solve();
    
    // pass result to laplace iterative solver
    for (auto v : patch->vertices()) {
        auto& vdata = patch->data(v);
        vdata.laplaceIterative.value = vdata.laplaceDirect.value;
    }
    
    // interleaved iterations of Laplacian smoothing + projection
    int n_iterations = max<int>(static_cast<int>(n_interior_vertices * 0.01), 5);
    for (int i = 0; i < n_iterations; ++i) {
        // laplace smoothing
        patch->laplaceIterative_compute(1);
        
        for (auto v : patch->vertices()) {
            if (patch->is_boundary(v))
                continue;
            
            auto& pn = patch->data(v).laplaceIterative.value;
                
            if (pn.tail(3).isZero())
                continue;
            
            // projection
            pn_normalize(pn);
            project(pn);
        }
    }
    
    // copy result back
    for (auto v : patch->vertices()) {
        if (patch->is_boundary(v))
            continue;
        
        patch->data(v).pn = patch->data(v).laplaceIterative.value;
    }
    
    patch->invalidate_displist();
}
void SketchRetopo::compute_patch_interior_pn_harmonic(curvenetwork::Patch* patch, BaseMesh::FHandle submesh_seed) {
    using Patch = curvenetwork::Patch;
    
    Patch::VHandle p_corner0;                           // Patch corner vertex with flag 0
    for (auto p_v : patch->vertices()) {                // For clarity, handles for Patch and BaseMesh are prefixed by p_ and b_, respectively.
        if (patch->data(p_v).patchgen.corner_index == 0) {
            p_corner0 = p_v;
            break;
        }
    }
    
    int num_sides = patch->param.get_num_sides();
    VectorXi l    = patch->param.l;
    
    // Step I: Compute harmonic parameterization for the patch |
    //---------------------------------------------------------+
    auto p_h = patch->halfedge_handle(p_corner0);
    for (int i = 0; i < num_sides; ++i) {
        for (int j = 0; j < l[i]; ++j, p_h = patch->prev_halfedge_handle(p_h)) {
            double t = i + j / static_cast<double>(l[i]);
            patch->data(patch->from_vertex_handle(p_h)).laplaceDirect.value << patchgen::get_boundary_geometry(num_sides, t), 0, 0, 0, 0; // Fixed boundary
        }
    }
    patch->laplaceDirect_solve();
    for (auto p_v : patch->vertices()) {
        auto& vdata = patch->data(p_v);
        vdata.harmonic_uv = vdata.laplaceDirect.value.head(2);
    }
    
    // Step II: Compute harmonic parameterization for the base mesh |
    //--------------------------------------------------------------+
    struct Constraint {                     // Tuple representing a constraint. Can handle face/edge/vertex-wise constraints.
        BaseMesh::VHandle vertex[3];
        double coord[3];
        Vector2d value;
    };
    vector<Constraint> constraints;
    
    // Clear floodfill flag
    for (auto f : basemesh.faces())
        basemesh.data(f).floodfill_flag = false;
    
    vector<BaseMesh::FHandle> submesh_faces;
    
    // For every patch boundary edge, set constraints
    p_h = patch->halfedge_handle(p_corner0);
    for (int i = 0; i < num_sides; ++i) {
        for (int j = 0; j < l[i]; ++j, p_h = patch->prev_halfedge_handle(p_h)) {
            auto p_v0v1 = patch->util_halfedge_to_vertex_pair(patch->opposite_halfedge_handle(p_h));
            auto p_v0 = p_v0v1.first;
            auto p_v1 = p_v0v1.second;
            
            // Face constraint
            auto hit = intersect(patch->data(p_v0).pn);
            auto b_f = basemesh.face_handle(hit.id0);
            auto b_fv = basemesh.util_fv_vec(b_f);
            constraints.push_back(Constraint{
                {b_fv[0]          ,     b_fv[1],    b_fv[2] },
                {1 - hit.u - hit.v,     hit.u  ,    hit.v   },
                patch->data(p_v0).harmonic_uv
            });
            
            // Traverse along the geodesics toward the next patch boundary vertex, adding edge constraints
            embree::Hit hit_next = intersect(patch->data(p_v1).pn);
    	    geodesic::SurfacePoint source(&geodesic_mesh.faces()[hit     .id0], 1 - hit     .u - hit     .v, hit     .u, hit     .v);
    	    geodesic::SurfacePoint target(&geodesic_mesh.faces()[hit_next.id0], 1 - hit_next.u - hit_next.v, hit_next.u, hit_next.v);
		    vector<geodesic::SurfacePoint> sources(1, source);
		    vector<geodesic::SurfacePoint> stop_points(1, target);
		    vector<geodesic::SurfacePoint> path;
            geodesic_algorithm.propagate(sources, geodesic::GEODESIC_INF, &stop_points);
		    geodesic_algorithm.trace_back(target, path);
            boost::reverse(path);                              // originally path is oriented from target to source
            double path_length_total = geodesic::length(path);
            double path_length_acc = 0;
            for (size_t j = 1; j < path.size() - 1; ++j) {
                geodesic::SurfacePoint& p = path[j];
                path_length_acc += p.distance(&path[j - 1]);
                double t = path_length_acc / path_length_total;
                Vector2d uv = (1 - t) * patch->data(p_v0).harmonic_uv +
                               t      * patch->data(p_v1).harmonic_uv;
                
                if (p.type() == geodesic::VERTEX) {
                    geodesic::vertex_pointer vertex = static_cast<geodesic::vertex_pointer>(p.base_element());
                    // vertex constraint
                    constraints.push_back(Constraint{
                        { basemesh.vertex_handle(vertex->id()), BaseMesh::VHandle(), BaseMesh::VHandle() },
                        { 1                                   , 0                  , 0                   },
                        uv
                    });
                
                } else if (p.type() == geodesic::EDGE) {
                    geodesic::edge_pointer edge = static_cast<geodesic::edge_pointer>(p.base_element());
                    auto b_v0 = basemesh.vertex_handle(edge->adjacent_vertices()[0]->id());
                    auto b_v1 = basemesh.vertex_handle(edge->adjacent_vertices()[1]->id());
                    auto p0 = o2e(basemesh.point(b_v0));
                    auto p1 = o2e(basemesh.point(b_v1));
                    double s = (Vector3d(p.x(), p.y(), p.z()) - p0).dot(p1 - p0) / (p1 - p0).squaredNorm();
                    // edge constraint
                    constraints.push_back(Constraint{
                        { b_v0 , b_v1, BaseMesh::VHandle() },
                        { 1 - s, s   , 0                   },
                        uv
                    });
                }
                
                // set floodfill flag for adjacent faces
                for (size_t k = 0; k < p.base_element()->adjacent_faces().size(); ++k) {
                    auto b_f = basemesh.face_handle(p.base_element()->adjacent_faces()[k]->id());
                    basemesh.data(b_f).floodfill_flag = true;
                    submesh_faces.push_back(b_f);
                }
            }
        }
    }
    
    // Collect basemesh face IDs inside the patch boundary by flood filling starting from the seed face
    const size_t n_submesh_faces_limit = submesh_faces.size() * 10000;
    vector<BaseMesh::FHandle> candidates = { submesh_seed };
    while (!candidates.empty()) {
        if (n_submesh_faces_limit < submesh_faces.size()) {
            cout << "Sub-mesh extraction process doesn't seem to terminate!\n";
            return;
        }
        auto b_f = candidates.back();
        candidates.pop_back();
        basemesh.data(b_f).floodfill_flag = true;
        submesh_faces.push_back(b_f);
        for (auto b_ff : basemesh.ff_range(b_f))
            if (!basemesh.data(b_ff).floodfill_flag)
                candidates.push_back(b_ff);
    }
    
    container_util::remove_duplicate(submesh_faces);
    
    // Mapping from vertices in the submesh to indices in the system matrix
    struct Hash {
        size_t operator()(const BaseMesh::VHandle& key) const { return key.idx(); }
    };
    unordered_map<BaseMesh::VHandle, int, Hash> vtx2var;
    for (auto b_f : submesh_faces)
        for (auto b_v : basemesh.fv_range(b_f))
            if (vtx2var.find(b_v) == vtx2var.end())
                vtx2var[b_v] = vtx2var.size();
    int n = vtx2var.size();
    
    // Construct Laplacian matrix
    vector<Triplet<double, int>> L_triplets;
    for (auto p_vtx_var : vtx2var) {
        auto b_v = p_vtx_var.first;
        int i    = p_vtx_var.second;
        
        double weight_sum = 0;
        for (auto b_ve : basemesh.ve_range(b_v)) {
            auto b_vv = basemesh.util_opposite_vertex(b_ve, b_v);
            auto j = container_util::at_optional(vtx2var, b_vv);
            if (!j) continue;
            //double weight = basemesh.cotanWeight(b_ve);           // Not sure if this is appropriate.
            double weight = 0;
            for (int k = 0; k < 2; ++k) {
                auto b_h = basemesh.halfedge_handle(b_ve, k);
                if (basemesh.is_boundary(b_h)) continue;
                if (basemesh.data(basemesh.face_handle(b_h)).floodfill_flag)
                    weight += basemesh.data(b_h).cotangent;
            }
            L_triplets.push_back( { i, *j, weight } );                            // Off-diagonal
            weight_sum += weight;
        }
        L_triplets.push_back({ i, i, -weight_sum });                              // Diagonal
    }
    SparseMatrix<double> L(n, n);
    L.setFromTriplets(L_triplets.begin(), L_triplets.end());
    
    // Weight for soft constraint
    double constraint_weight = 0;
    for (auto triplet : L_triplets)
        constraint_weight += abs(triplet.value());
    constraint_weight *= 0.1;
    
    // Specify soft constraints on the matrix and right hand side
    vector<Triplet<double, int>> C_triplets;
    MatrixX2d rhs = MatrixX2d::Zero(n, 2);
    for (auto constraint : constraints) {
        for (int k0 = 0; k0 < 3; ++k0) {
            auto i = container_util::at_optional(vtx2var, constraint.vertex[k0]);
            if (!i) continue;
            for (int k1 = 0; k1 < 3; ++k1) {
                auto j = container_util::at_optional(vtx2var, constraint.vertex[k1]);
                if (!j) continue;
                C_triplets.push_back({*i, *j, constraint_weight * constraint.coord[k0] * constraint.coord[k1]});
            }
            rhs.row(*i) += Vector2d(constraint_weight * constraint.coord[k0] * constraint.value);
        }
    }
    SparseMatrix<double> C(n, n);
    C.setFromTriplets(C_triplets.begin(), C_triplets.end());
    
    // Solve
    SparseMatrix<double> A = -L + C;
    //SparseMatrix<double> A = L * L + C;                 // This may be better? The sign could be wrong...
    SimplicialCholesky<SparseMatrix<double>> solver(A);
    MatrixX2d x = solver.solve(rhs);
    
    // Copy results
    for (auto p_vtx_var : vtx2var) {
        auto b_v = p_vtx_var.first;
        int i    = p_vtx_var.second;
        basemesh.data(b_v).harmonic_uv = x.row(i);
    }
    
    // Step III: For each patch interior vertex, find position on the basemesh with the same harmonic_uv |
    //---------------------------------------------------------------------------------------------------+
    for (auto p_v : patch->vertices()) {                            // Brute force search. Could be accelerated.
        if (patch->is_boundary(p_v)) continue;
        Vector2d p_uv = patch->data(p_v).harmonic_uv;
        
        for (auto b_f : submesh_faces) {
            auto b_fv = basemesh.util_fv_vec(b_f);
            Vector2d b_uv0 = basemesh.data(b_fv[0]).harmonic_uv;
            Vector2d b_uv1 = basemesh.data(b_fv[1]).harmonic_uv;
            Vector2d b_uv2 = basemesh.data(b_fv[2]).harmonic_uv;
            Vector3d t;
            if (!eigen_util::project_to_triangle(b_uv0, b_uv1, b_uv2, p_uv, t)) continue;
            if (t.minCoeff() < 0) continue;
            
            PointNormal b_pn0;      b_pn0 << o2e(basemesh.point(b_fv[0])), o2e(basemesh.normal(b_fv[0]));
            PointNormal b_pn1;      b_pn1 << o2e(basemesh.point(b_fv[1])), o2e(basemesh.normal(b_fv[1]));
            PointNormal b_pn2;      b_pn2 << o2e(basemesh.point(b_fv[2])), o2e(basemesh.normal(b_fv[2]));
            PointNormal& p_pn = patch->data(p_v).pn;
            p_pn = t[0] * b_pn0 +
                   t[1] * b_pn1 +
                   t[2] * b_pn2;

            pn_normalize(p_pn);
            project(p_pn);
            
            break;
        }
    }
    
    patch->invalidate_displist();
}

void SketchRetopo::add_uv_line(
    const Vector2d& uv0, const Vector2d& uv1,
    const boost::optional<Vector2d>& unnormalized_dir0, const boost::optional<Vector2d>& unnormalized_dir1)
{
    Polyline_PointNormal curve;
    int num_segments = static_cast<int>((uv1 - uv0).norm() / configTemp.segmentSize());
    boost::optional<Vector2d> a, b;
    if (unnormalized_dir0 && unnormalized_dir1) {
        double r = (uv1 - uv0).norm();
        Vector2d derivative0 = r * unnormalized_dir0->normalized();         // better not to use auto keyword with Eigen!
        Vector2d derivative1 = r * unnormalized_dir1->normalized();         // http://codepad.org/kZo6bPDY
        a =  uv0 - uv1 + derivative0;
        b = -uv0 + uv1 - derivative1;
    }
    for (int i = 0; i <= num_segments; ++i) {
        double t = i / static_cast<double>(num_segments);
        
        Vector2d uv;
        if (unnormalized_dir0 && unnormalized_dir1)
            // cubic spline interpolation
            uv = (1 - t) * uv0 + t * uv1 + t * (1 - t) * ((1 - t) * *a + t * *b);
        else
            uv = (1 - t) * uv0 + t * uv1;
        
        auto pn = uv_to_pn(uv);
        if (!pn)
            continue;
        
        project(*pn);
        curve.push_back(*pn);
    }
    
    if (curve.size() < 3)
        // this is weird, unusual situation
        return;
    
    add_curve_open_delegate(curve);
}

void SketchRetopo::add_spine_curve(const Polyline_PointNormal& spine_curve) {
    const double brush_size = configTemp.brushSize_spine();
    const double curve_length = spine_curve.length();
    const double angle_threshold = cos(3.1415 / 6);
    
    // set seed information for discrete exponential map
    vector<BaseMesh::VHandle> seeds;
    vector<Vector3d         > seeds_basis_u;
    vector<Vector2d         > seeds_uv     ;
    double seed_v = 0;
    for (int i = 0; i < spine_curve.size(); ++i) {
        auto hit = intersect(spine_curve[i]);
        if (!hit)
            // assume all curve points are really on surface
            return;
        
        // set seed id to one with highest barycentric coordinate
        auto f_hit = basemesh.face_handle(hit.id0);
        auto fv_hit = basemesh.fv_iter(f_hit);
        BaseMesh::VHandle v0_hit = *fv_hit;
        BaseMesh::VHandle v1_hit = *++fv_hit;
        BaseMesh::VHandle v2_hit = *++fv_hit;
        auto seed =
            bc(hit, 0) < bc(hit, 1) && bc(hit, 1) < bc(hit, 2) ? v2_hit :
            bc(hit, 0) < bc(hit, 1)                          ? v1_hit :
                                                             v0_hit ;
        
        // set seed basis direction
        Vector3d seed_normal  = o2e(basemesh.normal(seed));
        Vector3d seed_basis_v = Vector3d::Zero();
        if (0 < i                     ) seed_basis_v += (spine_curve[i    ] - spine_curve[i - 1]).head(3);
        if (i < spine_curve.size() - 1) seed_basis_v += (spine_curve[i + 1] - spine_curve[i    ]).head(3);
        seed_basis_v -= seed_normal.dot(seed_basis_v) * seed_normal;
        seed_basis_v.normalize();
        Vector3d seed_basis_u = seed_basis_v.cross(seed_normal);
        
        // set seed uv
        Vector2d seed_uv(0, seed_v);
        Vector3d d = o2e(basemesh.point(seed)) - spine_curve[i].head(3);
        seed_uv[0] += d.dot(seed_basis_u);
        seed_uv[1] += d.dot(seed_basis_v);
        if (i > 0) seed_v += pn_norm(spine_curve[i] - spine_curve[i - 1]);
        
        // add to array
        seeds        .push_back(seed        );
        seeds_basis_u.push_back(seed_basis_u);
        seeds_uv     .push_back(seed_uv     );
    }
    
    // parameterize the local surface area w/ discrete exponential map. this would cover parameterized area [-brush_size, brush_size] X [0, curve_length]
    Vector2d uv_min(-1.5 * brush_size, -0.5 * brush_size);
    Vector2d uv_max( 1.5 * brush_size, curve_length + 0.5 * brush_size);
    basemesh.expmap_compute(seeds, seeds_basis_u, seeds_uv, brush_size * 1.8, uv_min, uv_max);
    
    // debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    program_basemesh.enable();
    program_basemesh.set_uniform_1<float>("isoline_interval", static_cast<float>(brush_size));
    program_basemesh.disable();
    basemesh.invalidate_displist();
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug
    
    // halfedge flag indicating whether it is inside the brushed area or not.
    enum HalfedgeFlag {
        UNDETERMINED = 0,
        INVALID,                // 1: deleted or non-border
        SNAPPED_V,              // 2: snapped in v-direction
        OUTSIDE,                // 3: inside the parameterized area, but outside the brushed area
        INSIDE,                 // 4: inside the brushed area
        CONNECTED_TO_CORNER     // 5: connected to corner, could be slightly off the brushed area in v-direction, but considered as inside
    };
    // set invalid halfedge flag first
    for (auto& h : curvenetwork.halfedges) {
        h.flag = HalfedgeFlag::UNDETERMINED;
        if (h.is_deleted || h.patch() || !h.opposite->patch())
            h.flag = HalfedgeFlag::INVALID;
    }
    
    // determine corner configuration
    struct CornerInfo {
        bool is_snapped_u;
        bool is_snapped_v;
        curvenetwork::Halfedge* snapped_h;
        Vector2d uv;
        
        CornerInfo()
            : is_snapped_u(false)
            , is_snapped_v(false)
            , snapped_h(0)
        {}
    };
    CornerInfo corner_front_left ;
    CornerInfo corner_front_right;
    CornerInfo corner_back_left  ;
    CornerInfo corner_back_right ;
    
    // set default value for CornerInfo::uv
    for (int is_front = 0; is_front < 2; ++is_front)
    for (int is_left  = 0; is_left  < 2; ++is_left )
    {
        auto& corner = is_front ?
            (is_left ? corner_front_left : corner_front_right) :
            (is_left ? corner_back_left  : corner_back_right );
        corner.uv = Vector2d(is_left ? -brush_size : brush_size, is_front ? 0 : curve_length);
    }
    
    // check if corners snap in v-direction
    for (int is_front = 0; is_front < 2; ++is_front) {
        auto pn_endpoint = is_front ? spine_curve.front() : spine_curve.back();
        auto& corner_left  = is_front ? corner_front_left  : corner_back_left ;
        auto& corner_right = is_front ? corner_front_right : corner_back_right;
        
        double dist_min = DBL_MAX;
        curvenetwork::Halfedge* h_min = 0;
        for (auto& h : curvenetwork.halfedges) {
            if (h.flag != HalfedgeFlag::UNDETERMINED) continue;
            
            if (h.vertex->is_corner()) continue;        // don't snap to a corner
            
            double dist = pn_norm(h.vertex->pn - pn_endpoint);
            if (dist < dist_min) {
                dist_min = dist;
                h_min = &h;
            }
        }
        
        if (configTemp.snapSize() < dist_min) continue; // no snapping for this endpoint
        
        // walk forward on halfchain until hitting a corner or the brush region boundary
        auto h_forward = h_min;
        for (; ; h_forward = h_forward->next) {
            if (h_forward->vertex->is_corner())
                // stop when hitting a corner
                break;
            
            auto uv = pn_to_uv(h_forward->vertex->pn);
            if (!uv || abs(uv->x()) > brush_size) {
                // stop when exitting the parameterized area
                h_forward = h_forward->prev;
                break;
            }
        }
        
        // walk backward on halfchain until hitting a corner or the brush region boundary
        auto h_backward = h_min;
        for (; ; h_backward = h_backward->prev) {
            if (h_backward->vertex->is_corner())
                // stop when hitting a corner
                break;
            
            auto uv = pn_to_uv(h_backward->vertex->pn);
            if (!uv || abs(uv->x()) > brush_size) {
                // stop when exitting the parameterized area
                h_backward = h_backward->next;
                break;
            }
        }
        
        // set CornerInfo appropriately
        auto& corner_forward  = is_front ? corner_right : corner_left ;
        auto& corner_backward = is_front ? corner_left  : corner_right;
        
        corner_forward .is_snapped_v = true;
        corner_backward.is_snapped_v = true;
        corner_forward .snapped_h = h_forward ;
        corner_backward.snapped_h = h_backward;
        corner_forward .uv = *pn_to_uv(h_forward ->vertex->pn);
        corner_backward.uv = *pn_to_uv(h_backward->vertex->pn);
        
        // halfedges snapped in v-direction are excluded from testing snapping in u-direction
        for (auto h = h_forward; h != h_backward; h = h->prev)
            h->flag = HalfedgeFlag::SNAPPED_V;
    }
    
    // check if corners snap in u-directions
    for (int is_front = 0; is_front < 2; ++is_front)
    for (int is_left  = 0; is_left  < 2; ++is_left )
    {
        auto& corner = is_front ? (is_left ? corner_front_left : corner_front_right) : (is_left ? corner_back_left : corner_back_right);
        
        if (corner.is_snapped_v) {
            if (!corner.snapped_h->vertex->is_corner())
                // if snapped vertex is not a corner --> no snapping in u-direction
                continue;
            
            // check if snapped corner has an outgoing halfedge oriented similar to v-axis
            for (curvenetwork::VOHIter h(corner.snapped_h->vertex); h; ++h) {
                auto uv = pn_to_uv(h->vertex->pn);
                if (!uv) continue;
                
                *uv -= corner.uv;
                uv->normalize();
                if (is_front && uv->y() > angle_threshold || !is_front && uv->y() < -angle_threshold) {
                    corner.is_snapped_u = true;
                    break;
                }
            }
        
        } else {
            // check if there is a boundary vertex within the brush area
            double                  corner_dist_min = DBL_MAX;
            curvenetwork::Halfedge* corner_h_min = 0;
            double                  nocorner_dist_min = DBL_MAX;
            curvenetwork::Halfedge* nocorner_h_min = 0;
            for (auto& h : curvenetwork.halfedges) {
                if (h.flag != HalfedgeFlag::UNDETERMINED) continue;
                
                auto uv = pn_to_uv(h.vertex->pn);
                
                if (!uv) continue;      // outside parameterized area
                
                double dist_v = abs(is_front ?  uv->y() : uv->y() - curve_length);
                if (configTemp.snapSize() < dist_v) continue;        // too far from front/back side
                
                double dist_u = is_left ? -uv->x() : uv->x();
                if (brush_size < dist_u || dist_u <= 0) continue;           // too far from spine, or lying on the opposite side of spine
                
                auto is_corner = h.vertex->is_corner();
                if (is_corner && dist_u < corner_dist_min) {
                    corner_dist_min = dist_u;
                    corner_h_min = &h;
                } else if (!is_corner && dist_u < nocorner_dist_min) {
                    nocorner_dist_min = dist_u;
                    nocorner_h_min = &h;
                }
            }
            
            if (corner_h_min || nocorner_h_min) {
                // close vertex found --> snap to it
                corner.snapped_h = corner_h_min ? corner_h_min : nocorner_h_min;
                corner.uv = *pn_to_uv(corner.snapped_h->vertex->pn);
                corner.is_snapped_u = true;
            }
        }
    }
    
    // trace halfedges from snapped corner in v-direction
    for (int is_front = 0; is_front < 2; ++is_front)
    for (int is_left  = 0; is_left  < 2; ++is_left )
    {
        auto& corner          = is_front ? (is_left ? corner_front_left : corner_front_right) : (is_left ? corner_back_left  : corner_back_right );
        auto& corner_opposite = is_front ? (is_left ? corner_back_left  : corner_back_right ) : (is_left ? corner_front_left : corner_front_right);
        
        if (corner.is_snapped_u) {
            bool is_forward = is_front && !is_left || !is_front && is_left;
            
            // find halfedge outgoing from corner orienting similar to v-axis
            curvenetwork::Halfedge* h = 0;
            double angle_max = -DBL_MAX;
            for (curvenetwork::VOHIter voh(corner.snapped_h->vertex); voh; ++voh) {
                auto uv = pn_to_uv(voh->vertex->pn);
                if (!uv) continue;
                
                *uv -= corner.uv;
                uv->normalize();
                double angle = is_front ? uv->y() : -uv->y();
                if (angle_max < angle) {
                    angle_max = angle;
                    h = is_forward ? &*voh : voh->opposite;
                }
            }
            
            // set OUTSIDE flag to all undetermined one-ring halfedges but h
            for (curvenetwork::VOHIter voh(corner.snapped_h->vertex); voh; ++voh) {
                if (&*voh         != h && voh          ->flag == HalfedgeFlag::UNDETERMINED) voh          ->flag = HalfedgeFlag::OUTSIDE;
                if (voh->opposite != h && voh->opposite->flag == HalfedgeFlag::UNDETERMINED) voh->opposite->flag = HalfedgeFlag::OUTSIDE;
            }
            
            // trace halfedge while setting flag to CONNECTED_TO_CORNER. stop when exitting the parameterized area in u-direction
            // double traced_length = 0;
            while (true) {
                if (h->flag != HalfedgeFlag::UNDETERMINED || h == corner_opposite.snapped_h)
                    // reached the opposite side corner
                    break;
                
                h->flag = HalfedgeFlag::CONNECTED_TO_CORNER;
                if (h != corner.snapped_h)
                    h->is_corner = false;
                
                auto h_new = is_forward ? h->next : h->prev;
                auto uv     = pn_to_uv((is_forward ? h    ->vertex : h    ->from_vertex())->pn);
                auto uv_new = pn_to_uv((is_forward ? h_new->vertex : h_new->from_vertex())->pn);
                auto d = (*uv_new - *uv).normalized();
                double angle = is_front ? d.y() : -d.y();
                
                if (!uv_new ||
                    is_left && uv_new->x() < -brush_size ||
                    !is_left && brush_size < uv_new->x() ||
                    angle < angle_threshold)
                {
                    // outside the brushed area
                    break;
                }
                
                // traced_length += h->toVector3d().norm();
                h = h_new;
            }
        }
    }
    
    // halfedge coupled with UV coordinates of its two vertices
    struct HalfedgeUV {
        curvenetwork::Halfedge* h;
        boost::optional<Vector2d> uv_from;
        boost::optional<Vector2d> uv_to;
        
        bool operator<(const HalfedgeUV& rhs) const {
            if (!(*this) || !rhs) return false;
            
            return uv().y() < rhs.uv().y();
        }
        Vector2d uv() const { return uv_from ? *uv_from : *uv_to; }
        operator bool() const { return h && (uv_from || uv_to); }
    };
    
    // helper functions begin =================================================================================
    auto make_HalfedgeUV = [&] (curvenetwork::Halfedge* h) -> HalfedgeUV {
        HalfedgeUV result;
        result.h = h;
        result.uv_from = pn_to_uv(h->from_vertex()->pn);
        result.uv_to   = pn_to_uv(h->vertex       ->pn);
        return result;
    };
    
    auto is_inside = [&] (const HalfedgeUV& huv) -> bool {
        if (!huv.uv_from || !huv.uv_to) return false;
        
        for (int i = 0; i < 2; ++i) {
            auto& uv = i ? *huv.uv_from : *huv.uv_to;
            
            if (uv.x() < -brush_size || brush_size < uv.x() || uv.y() < 0 || curve_length < uv.y()) return false;
        }
        
        Vector2d d = *huv.uv_to - *huv.uv_from;
        d.normalize();
        
        return angle_threshold < abs(d.y());
    };
    // ================================================================================= helper functions end
    
    // set HalfedgeFlag::{ OUTSIDE | INSIDE }
    for (auto& h : curvenetwork.halfedges) {
        if (h.flag != HalfedgeFlag::UNDETERMINED) continue;
        
        auto huv = make_HalfedgeUV(&h);
        
        if (!huv || !is_inside(huv)) {
            h.flag = HalfedgeFlag::OUTSIDE;
            continue;
        }
        
        bool is_left = huv.uv().x() < 0;
        double v_min = (is_left ? corner_front_left : corner_front_right).uv.y();
        double v_max = (is_left ? corner_back_left  : corner_back_right ).uv.y();
        if (v_min < huv.uv().y() && huv.uv().y() < v_max)
            h.flag = HalfedgeFlag::INSIDE;
        else
            h.flag = HalfedgeFlag::OUTSIDE;
    }

    // sequence of halfedge which is inside and whose adjacent neighbor is outside.
    list<HalfedgeUV> huv_sequence_left;
    list<HalfedgeUV> huv_sequence_right;
    auto is_corner_vertex = [&] (curvenetwork::Vertex* v) {
        return
            corner_front_left .snapped_h && v == corner_front_left .snapped_h->vertex || 
            corner_front_right.snapped_h && v == corner_front_right.snapped_h->vertex || 
            corner_back_left  .snapped_h && v == corner_back_left  .snapped_h->vertex || 
            corner_back_right .snapped_h && v == corner_back_right .snapped_h->vertex;
    };
    for (auto& h : curvenetwork.halfedges) {
        if (h.flag < HalfedgeFlag::INSIDE) continue;    // h is not inside --> skip
        
        if (h.next->flag >= HalfedgeFlag::INSIDE && h.prev->flag >= HalfedgeFlag::INSIDE) {
            // both of its neighbor are also inside --> skip
            if (!is_corner_vertex(h      .vertex)) h      .is_corner = false;
            if (!is_corner_vertex(h.next->vertex)) h.next->is_corner = false;
            if (!is_corner_vertex(h.prev->vertex)) h.prev->is_corner = false;
            continue;
        }
        
        auto huv = make_HalfedgeUV(&h);
        bool is_left = huv.uv().x() < 0;
        
        if (h.next->flag < HalfedgeFlag::INSIDE && !is_corner_vertex(h.vertex)) {
            // next halfedge is outside
            auto huv_tmp = huv;
            huv_tmp.uv_from = boost::none;
            (is_left ? huv_sequence_left : huv_sequence_right).push_back(huv_tmp);
        }
        
        if (h.prev->flag < HalfedgeFlag::INSIDE && !is_corner_vertex(h.from_vertex())) {
            // prev halfedge is outside
            auto huv_tmp = huv;
            huv_tmp.uv_to = boost::none;
            (is_left ? huv_sequence_left : huv_sequence_right).push_back(huv_tmp);
        }
    }
    
    // sort huv_sequence based on its v-value
    huv_sequence_left .sort();
    huv_sequence_right.sort();
    
    // add curves in v-direction
    for (int is_left = 0; is_left < 2; ++is_left) {
        auto& corner_front = is_left ? corner_front_left : corner_front_right;
        auto& corner_back  = is_left ? corner_back_left  : corner_back_right ;
        auto& huv_sequence = is_left ? huv_sequence_left : huv_sequence_right;
        
        if (!corner_front.is_snapped_u && !corner_back.is_snapped_u && huv_sequence.empty()) {
            // no snapping in u-direction --> simply connect corner_front & corner_back
            // TODO: feature-aware curve modification
            add_uv_line(corner_front.uv, corner_back.uv);
            continue;
        }
        
        for (int is_front = 0; is_front < 2; ++is_front) {
            auto& corner = is_front ? corner_front : corner_back;
            
            if (!corner.is_snapped_u && !huv_sequence.empty()) {
                // corner_{front|back} is not snapped in u-direction --> connect huv_sequence.{front|back} and corner_{front|back}
                auto& huv = is_front ? huv_sequence.front() : huv_sequence.back();
                add_uv_line(corner.uv, huv.uv());
                
                // set t-junction flag
                (is_left && is_front || !is_left && !is_front ? huv.h : huv.h->prev)->is_corner = false;
                
                // for convenience in the later processing
                if (is_front) huv_sequence.pop_front();
                else          huv_sequence.pop_back ();
            }
        }
        
        // connect pairs in huv_sequence
        for (auto huv = huv_sequence.begin(); huv != huv_sequence.end(); ++huv) {
            auto uv0 = huv->uv();
            auto h0 = huv->h;
            ++huv;
            auto uv1 = huv->uv();
            auto h1 = huv->h;
            add_uv_line(uv0, uv1);
            (is_left ? h0->prev : h0)->is_corner = false;
            (is_left ? h1 : h1->prev)->is_corner = false;
        }
    }
    
    // add curves (u direction)
    for (int is_front = 0; is_front < 2; ++is_front) {
        auto& corner_left  = is_front ? corner_front_left  : corner_back_left ;
        auto& corner_right = is_front ? corner_front_right : corner_back_right;
        
        if (!corner_left.is_snapped_v) {
            assert(!corner_right.is_snapped_v);
            // TODO: feature-aware curve modification
            add_uv_line(corner_left.uv, corner_right.uv);
        }
    }
}

void SketchRetopo::compute_autocompl(const AutocmplInfo& autocmpl_info) {
    for (auto& uv_line : autocmpl_info.uv_lines) {
        add_uv_line(uv_line.uv0,
                    uv_line.uv1,
                    uv_line.unnormalized_dir0,
                    uv_line.unnormalized_dir1);
    }
}
vector<AutocmplInfo> SketchRetopo::suggest_autocmpl(curvenetwork::Halfedge* h, bool favor_triangle) {
    const double brush_size = configTemp.brushSize_autocmpl();
    
    // local parameterization using exponential map
    BaseMesh::VHandle seed;
    Vector3d seed_basis = h->toVector3d();
    Vector2d seed_uv = Vector2d::Zero();
    {
        auto hit = intersect(h->vertex->pn);
        if (!hit) return vector<AutocmplInfo>();
        // set seed id to one with highest barycentric coordinate
        auto f_hit  = basemesh.face_handle(hit.id0);
        auto fv_hit = basemesh.fv_iter(f_hit);
        BaseMesh::VHandle v0_hit = *fv_hit;
        BaseMesh::VHandle v1_hit = *++fv_hit;
        BaseMesh::VHandle v2_hit = *++fv_hit;
        seed =
            bc(hit, 0) < bc(hit, 1) && bc(hit, 1) < bc(hit, 2) ? v2_hit :
            bc(hit, 0) < bc(hit, 1)                          ? v1_hit :
                                                             v0_hit ;
        
        // set seed basis direction
        Vector3d seed_normal  = o2e(basemesh.normal(seed));
        eigen_util::orthonormalize(seed_normal, seed_basis);
        
        // set seed uv
        Vector3d d = o2e(basemesh.point(seed)) - h->vertex->pn.head(3);
        seed_uv[0] += d.dot(seed_basis);
        seed_uv[1] += d.dot(seed_normal.cross(seed_basis));
    }
    
    Vector2d uv_min = brush_size * Vector2d(-1.2, -0.3);
    Vector2d uv_max = brush_size * Vector2d( 1.2,  2.2);
    basemesh.expmap_compute(seed, seed_basis, seed_uv, brush_size * 3.8, uv_min, uv_max);
    
    // debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    program_basemesh.enable();
    program_basemesh.set_uniform_1<float>("isoline_interval", static_cast<float>(brush_size));
    program_basemesh.disable();
    basemesh.invalidate_displist();
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug
    
    vector<curvenetwork::Halfedge*> relevant_halfedges;
    for (auto& h : curvenetwork.halfedges) {
        if (h.is_deleted            ||
        
            h.patch()               ||
            !h.opposite->patch()    ||
            !pn_to_uv(h.vertex->pn) ||
            !pn_to_uv(h.from_vertex()->pn))
            continue;
        
        relevant_halfedges.push_back(&h);
    }
    
    vector<AutocmplInfo> result;
    
    const int num_steps = 10;
    for (int step = 1; step <= num_steps; ++step) {
        const double brush_size_current = step * brush_size / num_steps;
        AutocmplInfo autocmpl_info_current;
        
        enum HalfedgeFlag {
            OUTSIDE = 0,
            INSIDE,
            INSIDE_GROUPED,
        };
        
        auto is_inside = [&] (const Vector2d& uv) {
            return AlignedBox2d(Vector2d(-1.0, -0.25), Vector2d(1.0, 2.0)).contains(uv / brush_size_current);
        };
        
        // check all relevant halfedges inside the parameterized domain R
        for (auto h : relevant_halfedges) {
            h->flag = HalfedgeFlag::OUTSIDE;
            
            curvenetwork::Vertex* v[2] = {
                h->vertex,
                h->from_vertex()
            };
            
            for (int i = 0; i < 2; ++i) {
                auto uv = *pn_to_uv(v[i]->pn);
                if (is_inside(uv)) {
                    h->flag = HalfedgeFlag::INSIDE;
                    break;
                }
            }
        }
        
        struct HalfedgeSequence {
            curvenetwork::Halfedge* h_front;
            curvenetwork::Halfedge* h_back;
            
            HalfedgeSequence()
                : h_front()
                , h_back()
            {}
            int num_corners() const {
                int result = 0;
                for (auto h = h_front; h != h_back; h = h->next)
                    if (h->vertex->is_corner()) ++result;
                return result;
            }
            curvenetwork::Halfedge* get_corner(int i) const {
                int i_current = 0;
                for (auto h = h_front; h != h_back; h = h->next) {
                    if (h->vertex->is_corner()) {
                        if (i_current == i) return h;
                        ++i_current;
                    }
                }
                return 0;
            }
            PointNormal pn_front() const { return h_front->from_vertex()->pn; }
            PointNormal pn_back () const { return h_back ->vertex       ->pn; }
            PointNormal pn_corner(int i) const { return get_corner(i)->vertex->pn; }
            PointNormal sub_pn_front() const {
                int nc = num_corners();
                double dist_target = nc < 2
                    ? pn_norm(pn_back  ()       - pn_corner(0)     )
                    : pn_norm(pn_corner(nc - 1) - pn_corner(nc - 2));
                auto h_start = get_corner(0);
                auto pn_start = h_start->vertex->pn;
                for (auto h = h_start; h != h_front; h = h->prev) {
                    auto pn = h->from_vertex()->pn;
                    double dist = pn_norm(pn - pn_start);
                    if (dist_target < dist) return pn;
                }
                return pn_front();
            }
            PointNormal sub_pn_back() const {
                int nc = num_corners();
                double dist_target = nc < 2
                    ? pn_norm(pn_front ()  - pn_corner(0))
                    : pn_norm(pn_corner(1) - pn_corner(0));
                auto h_start = get_corner(nc - 1);
                auto pn_start = h_start->vertex->pn;
                for (auto h = h_start->next; h != h_back; h = h->next) {
                    auto pn = h->vertex->pn;
                    double dist = pn_norm(pn - pn_start);
                    if (dist_target < dist) return pn;
                }
                return pn_back();
            }
        };
        
        vector<HalfedgeSequence> halfedge_sequences;
        for (auto h : relevant_halfedges) {
            if (h->flag != HalfedgeFlag::INSIDE) continue;
            
            HalfedgeSequence s;
            // trace forward;
            s.h_back = h;
            while (true) {
                s.h_back->flag = HalfedgeFlag::INSIDE_GROUPED;
                if (s.h_back->next->flag != HalfedgeFlag::INSIDE) break;
                s.h_back = s.h_back->next;
            }
            // trace backward
            s.h_front = h;
            while (true) {
                s.h_front->flag = HalfedgeFlag::INSIDE_GROUPED;
                if (s.h_front->prev->flag != HalfedgeFlag::INSIDE) break;
                s.h_front = s.h_front->prev;
            }
            
            halfedge_sequences.push_back(s);
        }
        
        if (halfedge_sequences.empty() || halfedge_sequences.size() > 2)
            continue;
        
        const auto is_sharp      = [] (double angle) -> bool { return angle > 0 && angle < util::pi() * 0.66; };
        const auto is_very_sharp = [] (double angle) -> bool { return angle > 0 && angle < util::pi() * 0.33; };
        
        if (halfedge_sequences.size() == 1) {
            auto s = halfedge_sequences[0];
            int nc = s.num_corners();
            
            auto uv_front = pn_to_uv(s.pn_front());
            auto uv_back  = pn_to_uv(s.pn_back ());
            
            if (!uv_front || !uv_back) continue;
            
            auto uv0 = 0 < nc ? pn_to_uv(s.get_corner(0)->vertex->pn) : boost::none;
            auto uv1 = 1 < nc ? pn_to_uv(s.get_corner(1)->vertex->pn) : boost::none;
            auto uv2 = 2 < nc ? pn_to_uv(s.get_corner(2)->vertex->pn) : boost::none;
            auto uv3 = 3 < nc ? pn_to_uv(s.get_corner(3)->vertex->pn) : boost::none;
            
            auto sub_uv_front = nc != 0 ? pn_to_uv(s.sub_pn_front()) : boost::none;
            auto sub_uv_back  = nc != 0 ? pn_to_uv(s.sub_pn_back ()) : boost::none;
            
            if (nc == 0) {
                Vector2d d = eigen_util::rotate90(*uv_back - *uv_front).normalized();
                
                Vector2d uv_corner0 = *uv_front + (2 * brush_size_current) * d;
                Vector2d uv_corner1 = *uv_back  + (2 * brush_size_current) * d;
                
                auto pn_corner0 = uv_to_pn(uv_corner0);
                auto pn_corner1 = uv_to_pn(uv_corner1);
                if (!pn_corner0 || !pn_corner1) continue;

                // TODO: feature-aware curve modification
                autocmpl_info_current.add_uv_line(*uv_front , uv_corner0);
                autocmpl_info_current.add_uv_line(*uv_back  , uv_corner1);
                autocmpl_info_current.add_uv_line(uv_corner0, uv_corner1);
                autocmpl_info_current.set_corners(s.pn_front(),
                                                  s.pn_back (),
                                                  *pn_corner1,
                                                  *pn_corner0);
            } else if (nc == 1) {
                if (!uv0) continue;
                
                Vector2d d0 = *uv_front - *uv0;
                Vector2d d1 = *uv_back  - *uv0;
                double angle = eigen_util::angle(d1, d0);
                
                if (is_very_sharp(angle) || favor_triangle) {
                    if (d0.norm() < d1.norm()) {
                        autocmpl_info_current.add_uv_line(*uv_front, *sub_uv_back, eigen_util::rotate90(-d0), eigen_util::rotate90(-d1));
                        autocmpl_info_current.set_corners(s.pn_front(),
                                                          s.pn_corner(0),
                                                          s.sub_pn_back());
                    } else {
                        autocmpl_info_current.add_uv_line(*sub_uv_front, *uv_back, eigen_util::rotate90(-d0), eigen_util::rotate90(-d1));
                        autocmpl_info_current.set_corners(s.sub_pn_front(),
                                                          s.pn_corner(0),
                                                          s.pn_back());
                    }
                } else if (is_sharp(angle)) {
                    Vector2d uv_diag = *uv_front + *uv_back - *uv0;
                    
                    auto pn_diag = uv_to_pn(uv_diag);
                    if (!pn_diag) continue;
                    
                    autocmpl_info_current.add_uv_line(*uv_front, uv_diag);
                    autocmpl_info_current.add_uv_line(*uv_back , uv_diag);
                    autocmpl_info_current.set_corners(s.pn_front(),
                                                      s.pn_corner(0),
                                                      s.pn_back(),
                                                      *pn_diag);
                }
                
            } else if (nc == 2) {
                if (!uv0 || !uv1) continue;
                
                Vector2d d0 = (*uv0     - *uv_front).normalized();
                Vector2d d1 = (*uv1     - *uv0     ).normalized();
                Vector2d d2 = (*uv_back - *uv1     ).normalized();
                
                double angle0 = eigen_util::angle(d1, -d0);
                double angle1 = eigen_util::angle(d2, -d1);
                
                bool is_sharp0 = is_sharp(angle0);
                bool is_sharp1 = is_sharp(angle1);
                
                double r = (*uv1 - *uv0).norm();
                Vector2d uv_corner0 = *uv0 + r * eigen_util::rotate90(d1);
                Vector2d uv_corner1 = *uv1 + r * eigen_util::rotate90(d1);
                if (is_sharp1) uv_corner0 = *uv0 + r * d2;
                if (is_sharp0) uv_corner1 = *uv1 - r * d0;
                
                auto pn_corner0 = uv_to_pn(uv_corner0);
                auto pn_corner1 = uv_to_pn(uv_corner1);
                if (!pn_corner0 || !pn_corner1) continue;
                
                if (is_sharp0 && is_sharp1) {
                    autocmpl_info_current.add_uv_line(*uv_front, *uv_back);
                    autocmpl_info_current.set_corners(s.pn_front(),
                                                      s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_back());
                } else if (!is_sharp0 && is_sharp1) {
                    if (is_very_sharp(angle1) || favor_triangle) {
                        autocmpl_info_current.add_uv_line(*uv0, *sub_uv_back, eigen_util::rotate90(d1), eigen_util::rotate90(-d2));
                        autocmpl_info_current.set_corners(s.get_corner(0)->vertex->pn,
                                                          s.get_corner(1)->vertex->pn,
                                                          s.sub_pn_back());
                    } else {
                        autocmpl_info_current.add_uv_line(*uv0        , uv_corner0);
                        autocmpl_info_current.add_uv_line(*sub_uv_back, uv_corner0);
                        autocmpl_info_current.set_corners(s.pn_corner(0),
                                                          s.pn_corner(1),
                                                          s.sub_pn_back(),
                                                          *pn_corner0);
                    }
                    
                } else if (is_sharp0 && !is_sharp1) {
                    if (is_very_sharp(angle0) || favor_triangle) {
                        autocmpl_info_current.add_uv_line(*uv1, *sub_uv_front, eigen_util::rotate90(d1), eigen_util::rotate90(-d0));
                        autocmpl_info_current.set_corners(s.pn_corner(0),
                                                          s.pn_corner(1),
                                                          s.sub_pn_front());
                    } else {
                        autocmpl_info_current.add_uv_line(*uv1         , uv_corner1);
                        autocmpl_info_current.add_uv_line(*sub_uv_front, uv_corner1);
                        autocmpl_info_current.set_corners(s.pn_corner(0),
                                                          s.pn_corner(1),
                                                          *pn_corner1,
                                                          s.sub_pn_front());
                    }
                    
                } else if (!is_sharp0 && !is_sharp1) {
                    // TODO: feature-aware curve modification
                    autocmpl_info_current.add_uv_line(*uv0, uv_corner0);
                    autocmpl_info_current.add_uv_line(*uv1, uv_corner1);
                    autocmpl_info_current.add_uv_line(uv_corner0, uv_corner1);
                    autocmpl_info_current.set_corners(s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      *pn_corner1,
                                                      *pn_corner0);
                }
                
            } else if (nc == 3) {
                if (!uv0 || !uv1 || !uv2) continue;
                
                Vector2d d0 = *uv0     - *uv_front;
                Vector2d d1 = *uv1     - *uv0     ;
                Vector2d d2 = *uv2     - *uv1     ;
                Vector2d d3 = *uv_back - *uv2     ;
                
                double angle0 = eigen_util::angle(d1, -d0);
                double angle1 = eigen_util::angle(d2, -d1);
                double angle2 = eigen_util::angle(d3, -d2);
                
                bool is_sharp0 = is_sharp(angle0);
                bool is_sharp1 = is_sharp(angle1);
                bool is_sharp2 = is_sharp(angle2);
                
                if (!is_sharp1 || is_sharp0 && is_sharp2)
                    continue;
                
                if (is_very_sharp(angle1) || favor_triangle) {
                    autocmpl_info_current.add_uv_line(*uv0, *uv2, eigen_util::rotate90(d1), eigen_util::rotate90(-d2));
                    autocmpl_info_current.set_corners(s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_corner(2));
                } else if (is_sharp0) {
                    autocmpl_info_current.add_uv_line(*sub_uv_front, *uv2);
                    autocmpl_info_current.set_corners(s.sub_pn_front(),
                                                      s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_corner(2));
                } else if (is_sharp2) {
                    autocmpl_info_current.add_uv_line(*sub_uv_back, *uv0);
                    autocmpl_info_current.set_corners(s.sub_pn_back(),
                                                      s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_corner(2));
                } else {
                    auto uv_corner = *uv0 + *uv2 - *uv1;
                    
                    auto pn_corner = uv_to_pn(uv_corner);
                    if (!pn_corner) continue;
                    
                    autocmpl_info_current.add_uv_line(*uv0, uv_corner);
                    autocmpl_info_current.add_uv_line(*uv2, uv_corner);
                    autocmpl_info_current.set_corners(s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_corner(2),
                                                      *pn_corner);
                }
                
            } else if (nc == 4) {
                if (!uv0 || !uv1 || !uv2 || !uv3) continue;
                
                Vector2d d0 = *uv1 - *uv0;
                Vector2d d1 = *uv2 - *uv1;
                Vector2d d2 = *uv3 - *uv2;
                
                double angle0 = eigen_util::angle(d1, -d0);
                double angle1 = eigen_util::angle(d2, -d1);
                
                bool is_sharp0 = is_sharp(angle0);
                bool is_sharp1 = is_sharp(angle1);
                
                if (is_sharp0 && is_sharp1) {
                    autocmpl_info_current.add_uv_line(*uv0, *uv3, eigen_util::rotate90(d0), eigen_util::rotate90(-d2));
                    autocmpl_info_current.set_corners(s.pn_corner(0),
                                                      s.pn_corner(1),
                                                      s.pn_corner(2),
                                                      s.pn_corner(3));
                } else if (is_sharp0) {
                    if (is_very_sharp(angle0) || favor_triangle) {
                        autocmpl_info_current.add_uv_line(*uv0, *uv2, eigen_util::rotate90(d0), eigen_util::rotate90(-d1));
                        autocmpl_info_current.set_corners(s.pn_corner(0),
                                                          s.pn_corner(1),
                                                          s.pn_corner(2));
                    } else {
                        auto uv_corner = *uv0 + *uv2 - *uv1;
                        
                        auto pn_corner = uv_to_pn(uv_corner);
                        if (!pn_corner) continue;
                        
                        autocmpl_info_current.add_uv_line(*uv0, uv_corner);
                        autocmpl_info_current.add_uv_line(*uv2, uv_corner);
                        autocmpl_info_current.set_corners(s.pn_corner(0),
                                                          s.pn_corner(1),
                                                          s.pn_corner(2),
                                                          *pn_corner);
                    }
                } else if (is_sharp1) {
                    if (is_very_sharp(angle1) || favor_triangle) {
                        autocmpl_info_current.add_uv_line(*uv1, *uv3, eigen_util::rotate90(d1), eigen_util::rotate90(-d2));
                        autocmpl_info_current.set_corners(s.pn_corner(1),
                                                          s.pn_corner(2),
                                                          s.pn_corner(3));
                    } else {
                        auto uv_corner = *uv1 + *uv3 - *uv2;
                        
                        auto pn_corner = uv_to_pn(uv_corner);
                        if (!pn_corner) continue;
                        
                        autocmpl_info_current.add_uv_line(*uv1, uv_corner);
                        autocmpl_info_current.add_uv_line(*uv3, uv_corner);
                        autocmpl_info_current.set_corners(s.pn_corner(1),
                                                          s.pn_corner(2),
                                                          s.pn_corner(3),
                                                          *pn_corner);
                    }
                }
            }
        
        } else {    // halfedge_sequences.size() == 2
            auto& s0 = halfedge_sequences[0];
            auto& s1 = halfedge_sequences[1];
            
            if (s0.num_corners() < s1.num_corners()) swap(s0, s1);
            
            int nc0 = s0.num_corners();
            int nc1 = s1.num_corners();
            
            auto s0_front = pn_to_uv(s0.pn_front());        auto s0_back = pn_to_uv(s0.pn_back());  if (!s0_front || !s0_back) continue;
            auto s1_front = pn_to_uv(s1.pn_front());        auto s1_back = pn_to_uv(s1.pn_back());  if (!s1_front || !s1_back) continue;
            
            auto s0_uv0 = nc0 == 2 ? pn_to_uv(s0.get_corner(0)->vertex->pn) : boost::none;
            auto s0_uv1 = nc0 == 2 ? pn_to_uv(s0.get_corner(1)->vertex->pn) : boost::none;
            auto s1_uv0 = nc1 == 2 ? pn_to_uv(s1.get_corner(0)->vertex->pn) : boost::none;
            auto s1_uv1 = nc1 == 2 ? pn_to_uv(s1.get_corner(1)->vertex->pn) : boost::none;
            
            if (nc0 == 0 && nc1 == 0) {
                Vector2d d0 = *s0_back - *s0_front;
                Vector2d d1 = *s1_back - *s1_front;
                autocmpl_info_current.add_uv_line(*s0_front, *s1_back , eigen_util::rotate90(d0), eigen_util::rotate90(-d1));
                autocmpl_info_current.add_uv_line(*s0_back , *s1_front, eigen_util::rotate90(d0), eigen_util::rotate90(-d1));
                autocmpl_info_current.set_corners(s0.pn_front(),
                                                  s0.pn_back (),
                                                  s1.pn_front(),
                                                  s1.pn_back ());
            } else if (nc0 == 2 && nc1 == 0) {
                if (!s0_uv0 || !s0_uv1) continue;
                
                Vector2d d0 = (*s0_uv0  - *s0_front).normalized();
                Vector2d d1 = (*s0_uv1  - *s0_uv0  ).normalized();
                Vector2d d2 = (*s0_back - *s0_uv1  ).normalized();
                
                double angle0 = eigen_util::angle(d1, -d0);
                double angle1 = eigen_util::angle(d2, -d1);
                
                bool is_sharp0 = is_sharp(angle0);
                bool is_sharp1 = is_sharp(angle1);
                
                if (is_sharp0 || is_sharp1)
                    continue;
                
                Vector2d e0 = *s0_uv1  - *s0_uv0  ;
                Vector2d e1 = *s1_back - *s1_front;
                autocmpl_info_current.add_uv_line(*s0_uv0, *s1_back , eigen_util::rotate90(e0), eigen_util::rotate90(-e1));
                autocmpl_info_current.add_uv_line(*s0_uv1, *s1_front, eigen_util::rotate90(e0), eigen_util::rotate90(-e1));
                autocmpl_info_current.set_corners(s0.pn_corner(0),
                                                  s0.pn_corner(1),
                                                  s1.pn_front(),
                                                  s1.pn_back ());
            } else if (nc0 == 2 && nc1 == 2) {
                if (!s0_uv0 || !s0_uv1 || !s1_uv0 || !s1_uv1) continue;
                
                Vector2d s0_d0 = (*s0_uv0  - *s0_front).normalized();
                Vector2d s0_d1 = (*s0_uv1  - *s0_uv0  ).normalized();
                Vector2d s0_d2 = (*s0_back - *s0_uv1  ).normalized();
                Vector2d s1_d0 = (*s1_uv0  - *s1_front).normalized();
                Vector2d s1_d1 = (*s1_uv1  - *s1_uv0  ).normalized();
                Vector2d s1_d2 = (*s1_back - *s1_uv1  ).normalized();
                
                double s0_angle0 = eigen_util::angle(s0_d1, -s0_d0);
                double s0_angle1 = eigen_util::angle(s0_d2, -s0_d1);
                double s1_angle0 = eigen_util::angle(s1_d1, -s1_d0);
                double s1_angle1 = eigen_util::angle(s1_d2, -s1_d1);
                
                bool is_s0_sharp0 = is_sharp(s0_angle0);
                bool is_s0_sharp1 = is_sharp(s0_angle1);
                bool is_s1_sharp0 = is_sharp(s1_angle0);
                bool is_s1_sharp1 = is_sharp(s1_angle1);
                
                if (is_s0_sharp0 || is_s0_sharp1 || is_s1_sharp0 || is_s1_sharp1)
                    continue;
                
                Vector2d e0 = *s0_uv1 - *s0_uv0;
                Vector2d e1 = *s1_uv1 - *s1_uv0;
                autocmpl_info_current.add_uv_line(*s0_uv0, *s1_uv1, eigen_util::rotate90(e0), eigen_util::rotate90(-e1));
                autocmpl_info_current.add_uv_line(*s0_uv1, *s1_uv0, eigen_util::rotate90(e0), eigen_util::rotate90(-e1));
                autocmpl_info_current.set_corners(s0.pn_corner(0),
                                                  s0.pn_corner(1),
                                                  s1.pn_corner(0),
                                                  s1.pn_corner(1));
            }
        }
        
        if (!autocmpl_info_current.uv_lines.empty())
            result.push_back(autocmpl_info_current);
    }
    
    return result;
}

void SketchRetopo::add_cylinder_curve(Polyline_PointNormal cylinder_curve) {
    // set flags for all halfedges that are on patch boundaries and form loops
    curvenetwork.set_flag_halfedges(curvenetwork::Core::FlagOp::SUBST, 0);
    for (auto& c : curvenetwork.halfchains) {
        if (c.patch)
            // this halfchain is not on a patch boundary
            continue;
        
        if (c.halfedge_front->flag)
            // we already know it's on a loop
            continue;
        
        // see if this halfchain is on a loop
        curvenetwork::Halfchain* c_front(nullptr);
        curvenetwork::Halfchain* c_back (nullptr);
        curvenetwork::trace_halfchains(&c, c_front, c_back);
        if (c_front->prev() == c_back) {
            // set flags for all halfedges on the same loop
            int num_corner = 0;
            for (auto h = c.halfedge_front; ; ) {
                h->flag = 1;
                if (h->vertex->is_corner())
                    ++num_corner;
                h = h->next;
                if (h == c.halfedge_front)
                    break;
            }
            if (num_corner == 0) {
                // if there's no corner on this loop, set a special flag for each halfedge
                for (auto h = c.halfedge_front; ; ) {
                    h->flag = 2;
                    h = h->next;
                    if (h == c.halfedge_front)
                        break;
                }
            }
        }
    }
    
    // find a halfedge each endpoint snaps to
    curvenetwork::Halfedge* h_snapped[2] = { 0, 0 };
    for (int i = 0; i < 2; ++i) {
        auto& pn_endpoint = i == 0 ? cylinder_curve.front() : cylinder_curve.back();
        Vector3d d;
        if (i == 0) d = (*++cylinder_curve. begin() - *cylinder_curve. begin()).head(3);
        else        d = (*++cylinder_curve.rbegin() - *cylinder_curve.rbegin()).head(3);
        d.normalize();
        
        double dist_min = util::dbl_max();
        for (auto& h : curvenetwork.halfedges) {
            if (h.is_deleted || !h.flag) continue;
            
            if (h.flag == 1 && !h.vertex->is_corner())
                // if this loop has any corners, endpoint should always snap to corners
                continue;
            
            double dist = pn_norm(h.vertex->pn - pn_endpoint);
            
            double dp = h.toVector3d().cross(d).dot(pn_endpoint.tail(3));       // to avoid choosing halfedge on the wrong side
            
            if (dist < dist_min && dist < configTemp.snapSize() && dp > 0) {
                dist_min = dist;
                h_snapped[i] = &h;
            }
        }
    }
    
    if (!h_snapped[0] || !h_snapped[1])
        // both endpoints have to be snapped!
        return;
    
    // make sure these loops are different ones, check number of corners
    int num_corners[2] = { 0, 0 };
    for (int i = 0; i < 2; ++i) {
        for (auto h = h_snapped[i]; ; ) {
            if (h == h_snapped[(i + 1) % 2]) {
                // two halfedges are on the same edge flow!
                h_snapped[0] = 0;
                break;
            }
            
            if (h->vertex->is_corner())
                ++num_corners[i];
            
            h = h->next;
            if (h == h_snapped[i])
                break;
        }
        
        if (!h_snapped[0])
            return;
    }
    
    // check if corners on the loops match
    if (num_corners[0] && num_corners[1] && num_corners[0] != num_corners[1])
        return;
    
    // obtain loop curve geometry
    Polyline_PointNormal loop_curve[2];
    for (int i = 0; i < 2; ++i) {
        for (auto h = h_snapped[i]; ; ) {
            loop_curve[i].push_back(h->vertex->pn);
            
            h = i == 0 ? h->next : h->prev;
            if (h == h_snapped[i])
                break;
        }
        loop_curve[i].is_loop = true;
    }
    
    // reorder such that num_corners[1] is equal or greater than num_corner[0]
    if (num_corners[0] > num_corners[1]) {
        reverse(cylinder_curve.begin(), cylinder_curve.end());
        swap(h_snapped  [0], h_snapped  [1]);
        swap(num_corners[0], num_corners[1]);
        swap(loop_curve [0], loop_curve [1]);
        loop_curve[0].reverse();
        loop_curve[1].reverse();
    }
    
    vector<PointNormal> geodesic_endpoint[2];
    
    if (num_corners[0] > 0) {
        assert(num_corners[0] == num_corners[1]);
        
        // just connect corresponding corners
        for (int i = 0; i < 2; ++i) {
            for (auto h = i == 0 ? h_snapped[i]->next : h_snapped[i]->prev; ; ) {
                if (h->vertex->is_corner())
                    geodesic_endpoint[i].push_back(h->vertex->pn);
                
                if (h == h_snapped[i])
                    break;
                h = i == 0 ? h->next : h->prev;
            }
        }
        
    } else if (num_corners[1] > 0) {
        assert(!num_corners[0]);
        
        // get arc length parameter from num_corners[1] side
        vector<double> arc_length_parameter(num_corners[1], 0);
        auto c = h_snapped[1]->halfchain;
        for (int i = 0; i < num_corners[1]; ++i, c = c->prev())
            arc_length_parameter[i] = (i == 0 ? 0 : arc_length_parameter[i - 1]) + c->length();
        assert(c == h_snapped[1]->halfchain);
        
        for (int i = 0; i < num_corners[1]; ++i)
            arc_length_parameter[i] /= arc_length_parameter.back();
        
        // sample point on loop_curve[0] according to the above arc length parameter
        for (int i = 0; i < num_corners[1]; ++i)
            geodesic_endpoint[0].push_back(loop_curve[0].point_at(arc_length_parameter[i]));
        
        // sample corners on loop_curve[1]
        for (auto h = h_snapped[1]->prev; ; ) {
            if (h->vertex->is_corner())
                geodesic_endpoint[1].push_back(h->vertex->pn);
            
            if (h == h_snapped[1])
                break;
            h = h->prev;
        }
    
    } else {
        assert(!num_corners[0] && !num_corners[1]);
        
        for (int i = 0; i < 2; ++i) {
            for (int j = 1; j <= configTemp.cylinder_num_div; ++j)
                geodesic_endpoint[i].push_back(loop_curve[i].point_at(j / static_cast<double>(configTemp.cylinder_num_div)));
        }

        assert(!num_corners[0] && !num_corners[1]);
    }
    assert(geodesic_endpoint[0].size() == geodesic_endpoint[1].size());
    
    // last points correspond to cylinder_curve --> omit
    geodesic_endpoint[0].pop_back();
    geodesic_endpoint[1].pop_back();
    
    memento_store();
    
    // add geodesic traced curve
    loop_util::for_each(
        geodesic_endpoint[0].begin(), geodesic_endpoint[0].end(), geodesic_endpoint[1].begin(), 
        [&] (const PointNormal& pn_from, const PointNormal& pn_to) {
            auto geodesic_curve = geodesic_compute(pn_from, pn_to);
            if (geodesic_curve.empty())
                return;
            
            // reampling, projection
            geodesic_curve.resample_by_length(configTemp.segmentSize());
            for_each(
                geodesic_curve.begin(), geodesic_curve.end(),
                [&] (PointNormal& pn) { project(pn); });
            
            add_curve_open_delegate(geodesic_curve);
    });
    
    // add sketched curve
    add_curve_open_delegate(cylinder_curve);
    
    // check if the opposite side of each snapped halfedge forms a reasonable loop. if yes, activate the corner flag
    for (int i = 0; i < 2; ++i) {
        if (i == 0)
            loop_curve[i].reverse();

        if (!is_loop_ccw(loop_curve[i]))
            continue;
        
        curvenetwork::Halfedge* h_start = h_snapped[i]->opposite;
        for (auto h = h_start; ; ) {
            if (h->vertex->is_corner())
                h->is_corner = true;
            
            h = h->next;
            if (h == h_start)
                break;
        }
    }
}

Polyline_PointNormal SketchRetopo::trace_on_plane(const PointNormal& pn_start, const Vector3d& plane_normal, double dist_max) {
    const Polyline_PointNormal empty_result;
    
    // intersected basemesh face
    auto hit = intersect(pn_start);
    if (!hit)
        return empty_result;
    
    // helper functions
    auto plane_func = [pn_start, plane_normal] (const Vector3d& p) { return plane_normal.dot(p - pn_start.head(3)); };
    
    auto is_crossing = [plane_func] (const BaseMesh& mesh, const BaseMesh::HHandle h) {
        auto p01 = mesh.util_halfedge_to_point_pair(h);
        return plane_func(o2e(p01.first)) * plane_func(o2e(p01.second)) <= 0;
    };
    
    // find a starting halfedge
    BaseMesh::HHandle h_start;
    for (auto h : basemesh.fh_range(basemesh.face_handle(hit.id0))) {
        if (is_crossing(basemesh, h)) {
            h_start = h;
            break;
        }
    }
    if (!h_start.is_valid())
        return empty_result;
    
    // check if the chain of halfedges is open or closed
    for (auto h = h_start; ; ) {
        auto h1 = basemesh.next_halfedge_handle(h);
        auto h2 = basemesh.next_halfedge_handle(h1);
        
        auto h_next = is_crossing(basemesh, h1) ? h1 : h2;
        h_next = basemesh.opposite_halfedge_handle(h_next);
        
        if (basemesh.is_boundary(h_next)) {
            // reached a boundary
            h_start = basemesh.opposite_halfedge_handle(h_next);
            break;
        
        } else if (h_next == h_start) {
            // loop detected!
            break;
        }
        
        h = h_next;
    }
    
    // start walking on basemesh
    Polyline_PointNormal result;
    double dist_sum = 0;
    for (auto h = h_start; ; ) {
        // get edge segment
        auto v0v1 = basemesh.util_halfedge_to_vertex_pair(h);
        PointNormal pn0 = basemesh.get_pn(v0v1.first );
        PointNormal pn1 = basemesh.get_pn(v0v1.second);
        
        // get intersected point
        double t0 = plane_func(pn0.head(3));
        double t1 = plane_func(pn1.head(3));
        assert(t0 * t1 <= 0);
        double u = - t0 / (t1 - t0);
        PointNormal pn = (1 - u) * pn0 + u * pn1;
        
        pn_normalize(pn);
        project(pn);
        
        if (!result.empty()) {
            dist_sum += pn_norm(result.back() - pn);
            if (dist_max < dist_sum)
                // have traveled too far
                return empty_result;
        }
        
        result.push_back(pn);
        
        // look for next halfedge
        auto h1 = basemesh.next_halfedge_handle(h);
        auto h2 = basemesh.next_halfedge_handle(h1);
        
        auto h_next = is_crossing(basemesh, h1) ? h1 : h2;
        h_next = basemesh.opposite_halfedge_handle(h_next);
        
        if (basemesh.is_boundary(h_next) || h_next == h_start) {        // reached a boundary, or loop detected
            if (h_next == h_start)
                result.is_loop = true;
            
            // resampling
            //result.resample_by_length(configTemp.segmentSize());
            result.insert_points_per_segment(9);
            project(result);
            
            return result;
        }
        
        h = h_next;
    }
    
    // should not reach here
    return empty_result;
}

void                 SketchRetopo::edgeLoop_insert (curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start, double t) {
    if (!patch_start || !h_start.is_valid() || !patch_start->is_boundary(h_start) || t <= 0 || 1 <= t)
        return;
    
    memento_store();
    
    // check if this global (i.e., inter-patch) edge loop is open or closed
    auto patch = patch_start;
    auto h     = h_start;
    while (true) {
        auto edgeLoop = find_edgeloop(*patch, h);
        assert(edgeLoop.front() != edgeLoop.back());
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        assert(patch->is_boundary(h_opposite));
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        if (!patch_next || patch_next->is_imaginary()) {
            // reached a boundary
            patch_start = patch;
            h_start     = h_opposite;
            t = 1 - t;
            break;
        } else if (patch_next == patch_start && h_next == h_start) {
            // loop detected
            break;
        }
        
        patch = patch_next;
        h     = h_next;
    }
    
    vector<curvenetwork::Edgechain*> e_affected;
    vector<curvenetwork::Patch    *> patch_affected;
    
    patch = patch_start;
    h     = h_start;
    while (true) {
        auto edgeLoop = find_edgeloop(*patch, h);
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        
        e_affected    .push_back(patch->data(h         ).halfchain->edgechain);
        e_affected    .push_back(patch->data(h_opposite).halfchain->edgechain);
        patch_affected.push_back(patch);
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        insert_edgeloop<curvenetwork::Patch>(*patch, h, t,
            [&] (curvenetwork::Patch& patch_, OpenMesh::HalfedgeHandle h_old, OpenMesh::VertexHandle v_new) {
                auto v0v1 = patch_.util_halfedge_to_vertex_pair(h_old);
                auto& vdata0    = patch_.data(v0v1.first );
                auto& vdata1    = patch_.data(v0v1.second);
                auto& vdata_new = patch_.data(v_new      );
                
                vdata_new.pn = (1 - t) * vdata0.pn + t * vdata1.pn;
                pn_normalize(vdata_new.pn);
                project(vdata_new.pn);
        });
        
        if (!patch_next || patch_next->is_imaginary() || patch_next == patch_start && h_next == h_start)
            break;
        patch = patch_next;
        h     = h_next;
    }
    
    container_util::remove_duplicate(e_affected);
    for_each(e_affected.begin(), e_affected.end(),
        [] (curvenetwork::Edgechain* e) { ++e->num_subdiv; });
    
    container_util::remove_duplicate(patch_affected);
    for_each(patch_affected.begin(), patch_affected.end(),
        [] (curvenetwork::Patch* patch_) {
            patch_->set_halfedgeData();
            patch_->invalidate_displist();
    });
}
Polyline_PointNormal SketchRetopo::edgeLoop_preview(curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start, double t) const {
    auto patch_start_orig = patch_start;
    auto h_start_orig = h_start;
    double t_orig = t;
    
    if (!patch_start || !h_start.is_valid() || !patch_start->is_boundary(h_start) || t <= 0 || 1 <= t)
        return Polyline_PointNormal();
    
    Polyline_PointNormal result;
    
    // check if this global (i.e., inter-patch) edge loop is open or closed
    auto patch = patch_start;
    auto h     = h_start;
    while (true) {
        if (!h.is_valid())
            // bug!
            return result;
        
        auto edgeLoop = find_edgeloop(*patch, h);
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        if (!patch_next || patch_next->is_imaginary()) {
            // reached a boundary
            patch_start = patch;
            h_start     = h_opposite;
            t = 1 - t;
            break;
        } else if (patch_next == patch_start && h_next == h_start) {
            // loop detected
            break;
        }
        
        patch = patch_next;
        h     = h_next;
    }
    
    patch = patch_start;
    h     = h_start;
    while (true) {
        auto edgeLoop = find_edgeloop(*patch, h);
        
        for_each(
            edgeLoop.rbegin(), edgeLoop.rend(),
            [&] (curvenetwork::Patch::HHandle h) {
                auto v0v1 = patch->util_halfedge_to_vertex_pair(h);
                auto& pn0 = patch->data(v0v1.first ).pn;
                auto& pn1 = patch->data(v0v1.second).pn;
                PointNormal pn = (1 - t) * pn0 + t * pn1;
                pn_normalize(pn);
                result.push_back(pn);
        });
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        if (!patch_next || patch_next->is_imaginary() || patch_next == patch_start && h_next == h_start)
            break;
        patch = patch_next;
        h     = h_next;
    }
    
    return result;
}

vector<pair<curvenetwork::Patch*, curvenetwork::Patch::HHandle>> SketchRetopo::edgeLoop_walk(curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start) const {
    typedef vector<pair<curvenetwork::Patch*, curvenetwork::Patch::HHandle>> Result;
    
    auto patch_start_orig = patch_start;
    auto h_start_orig = h_start;
    
    if (!patch_start || !h_start.is_valid() || !patch_start->is_boundary(h_start))
        return Result();
    
    Result result;
    
    // check if this global (i.e., inter-patch) edge loop is open or closed
    auto patch = patch_start;
    auto h     = h_start;
    while (true) {
        auto edgeLoop = find_edgeloop(*patch, h);
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        if (!patch_next || patch_next->is_imaginary()) {
            // reached a boundary
            patch_start = patch;
            h_start     = h_opposite;
            break;
        } else if (patch_next == patch_start && h_next == h_start) {
            // loop detected
            break;
        }
        
        patch = patch_next;
        h     = h_next;
    }
    
    patch = patch_start;
    h     = h_start;
    while (true) {
        result.push_back(make_pair(patch, h));
        
        auto edgeLoop = find_edgeloop(*patch, h);
        
        auto h_opposite = patch->opposite_halfedge_handle(edgeLoop.front());
        
        auto patch_next = patch->opposite_patch            (h_opposite);
        auto h_next     = patch->opposite_boundary_halfedge(h_opposite);
        
        if (!patch_next || patch_next->is_imaginary()) {
            result.push_back(make_pair(patch, h_opposite));
            break;
        }
        
        if (patch_next == patch_start && h_next == h_start)
            break;
        
        patch = patch_next;
        h     = h_next;
    }
    
    return result;
}

void SketchRetopo::vertices_move  (const PointNormal& pn_center, const Vector3d& center_offset, double radius) {
    for (auto& patch : curvenetwork.patches) {
        for (auto v : patch.vertices()) {
            auto& pn = patch.data(v).pn;
            
            double dist = pn_norm(pn_center - pn);
            double dot  = pn_center.tail(3).dot(pn.tail(3));
            if (dist > radius || dot < 0)
                continue;
            
            pn.head(3) += RBFKernel_Wendland()(dist / radius) * center_offset;
            project(pn);
            
            // check if v is on the symmetry plane
            bool is_symmetric = false;
            if (configSaved.symmetric && patch.is_boundary(v)) {
                auto h = patch.halfedge_handle(v);
                curvenetwork::Halfchain* c[2] = {
                    patch.data(h                            ).halfchain,
                    patch.data(patch.prev_halfedge_handle(h)).halfchain
                };
                for (int i = 0; i < 2; ++i) {
                    if (c[i]->opposite()->patch && c[i]->opposite()->patch->is_imaginary_symmetry()) {
                        is_symmetric = true;
                        break;
                    }
                }
            }
            
            if (is_symmetric)
                pn[0] = pn[3] = 0;
            
            patch.invalidate_displist();
        }
    }
    
    vector<curvenetwork::Edgechain*> affected_edgechains;
    for (auto& v : curvenetwork.vertices) {
        double dist = pn_norm(pn_center - v.pn);
        double dot = pn_center.tail(3).dot(v.pn.tail(3));
        if (dist > radius || dot < 0)
            continue;
        
        for (curvenetwork::VOCIter c(&v); c; ++c)
            affected_edgechains.push_back(c->edgechain);
    }
    
    container_util::remove_duplicate(affected_edgechains);
    for (auto e : affected_edgechains) {
        for (int i = 0; i < 2; ++i) {
            if (curvenetwork::Patch::is_ordinary(e->halfchain[i]->patch)) {
                set_halfchain_pn_from_patch(e->halfchain[i]);
                break;
            }
        }
    }
}
void SketchRetopo::vertices_smooth(const PointNormal& pn_center, double radius) {
    vector<curvenetwork::Edgechain*> affected_edgechains;
    
    for (auto& patch : curvenetwork.patches) {
        for (auto patch_v : patch.vertices()) {
            auto& vdata = patch.data(patch_v);
            vdata.pn_temp = vdata.pn;
            
            // too far or too differently oriented vertices are skipped
            double dist = pn_norm(vdata.pn - pn_center);
            double dot = vdata.pn.tail(3).dot(pn_center.tail(3));
            if (dist > radius || dot < 0)
                continue;
            
            auto average = vertices_smooth_sub(&patch, patch_v, affected_edgechains);
            if (!average)
                continue;
            
            // store falloff-weighted result into pn_temp
            double falloff = RBFKernel_Wendland()(dist / radius);
            vdata.pn_temp = (1 - falloff) * vdata.pn + falloff * (*average);
            pn_normalize(vdata.pn_temp);
            
            // projection
            project(vdata.pn_temp);
            
            // update diplay list
            patch.invalidate_displist();
        }
    }
    
    // copy result from pn_temp
    for (auto& patch : curvenetwork.patches) {
        for (auto v : patch.vertices()) {
            auto& vdata = patch.data(v);
            vdata.pn = vdata.pn_temp;
        }
    }
    
    // propagate position update to curve network vertices
    container_util::remove_duplicate(affected_edgechains);
    for (auto e : affected_edgechains) {
        for (int i = 0; i < 2; ++i) {
            if (curvenetwork::Patch::is_ordinary(e->halfchain[i]->patch)) {
                set_halfchain_pn_from_patch(e->halfchain[i]);
                break;
            }
        }
    }
}
auto SketchRetopo::vertices_smooth_sub(curvenetwork::Patch* patch, curvenetwork::Patch::VHandle patch_v, std::vector<curvenetwork::Edgechain*>& affected_edgechains) -> boost::optional<PointNormal> const {
    PointNormal sum_pn     = PointNormal::Zero();
    double      sum_weight = 0;
    
    if (patch->is_boundary(patch_v)) {
        if (patch->is_failure)
            return boost::none;
        
        auto v = patch->vertex_patch_to_curveNetwork(patch_v);
        if (v) {
            // patch corner vertex
            if (v->is_boundary())
                // if this corner is on the boundary of curve network --> skip it
                return boost::none;
            
            // walk over the connected patches
            for (curvenetwork::VOCIter c(v); c; ++c) {
                auto c_patch = c->patch;
                if (c_patch->is_failure)
                    continue;
                
                auto patch_v2 = c_patch->vertex_curveNetwork_to_patch(v);
                if (!patch_v2.is_valid())
                    // bug!
                    continue;
                
                // sum up adjacent patch vertices
                for (auto w : c_patch->vv_range(patch_v2)) {
                    double weight = c_patch->is_boundary(w) ? 0.5 : 1.0;
                    sum_pn     += weight * c_patch->data(w).pn;
                    sum_weight += weight;
                }
                affected_edgechains.push_back(c->edgechain);
            }
            
        } else {
            // patch side vertex
            auto patch_h = patch->prev_halfedge_handle(patch->halfedge_handle(patch_v));
            auto c = patch->data(patch_h).halfchain;
            
            affected_edgechains.push_back(c->edgechain);
            
            auto& pn_prev = patch->data(patch->from_vertex_handle(patch_h                             )).pn;
            auto& pn_next = patch->data(patch->to_vertex_handle  (patch->next_halfedge_handle(patch_h))).pn;
            
            sum_pn     += pn_prev + pn_next;
            sum_weight += 2;
            
            auto patch_opposite = patch->opposite_patch(patch_h);
            if (curvenetwork::Patch::is_ordinary(patch_opposite)) {
                // if vertex is not on the boundary of curve network, also add up patch interior vertices
                for (auto w : patch->vv_range(patch_v)) {
                    if (patch->is_boundary(w))
                        continue;
                    sum_pn     += patch->data(w).pn;
                    sum_weight += 1;
                }
                
                // from opposite patch
                auto patch_opposite_h = patch->opposite_boundary_halfedge(patch_h);
                auto patch_opposite_v = patch_opposite->from_vertex_handle(patch_opposite_h);
                for (auto w : patch_opposite->vv_range(patch_opposite_v)) {
                    if (patch_opposite->is_boundary(w))
                        continue;
                    sum_pn     += patch_opposite->data(w).pn;
                    sum_weight += 1;
                }
            }
        }
        
    } else {
        // patch interior vertex
        for (auto w : patch->vv_range(patch_v)) {
            sum_pn     += patch->data(w).pn;
            sum_weight += 1;
        }
        
    }
    
    return sum_pn /= sum_weight;
}
void SketchRetopo::vertices_smooth_global() {
    for (auto& patch : curvenetwork.patches) {
        for (auto patch_v : patch.vertices()) {
            auto& vdata = patch.data(patch_v);
            vdata.pn_temp = vdata.pn;
            
            vector<curvenetwork::Edgechain*> affected_edgechains;
            auto average = vertices_smooth_sub(&patch, patch_v, affected_edgechains);
            if (!average)
                continue;
            
            vdata.pn_temp = *average;
            pn_normalize(vdata.pn_temp);
            project(vdata.pn_temp);
        }
        patch.invalidate_displist();
    }
    // copy result from pn_temp
    for (auto& patch : curvenetwork.patches) {
        for (auto patch_v : patch.vertices()) {
            auto& vdata = patch.data(patch_v);
            vdata.pn = vdata.pn_temp;
        }
    }
    
    // propagate position update to curve network vertices
    for (auto& e : curvenetwork.edgechains) {
        for (int i = 0; i < 2; ++i) {
            if (curvenetwork::Patch::is_ordinary(e.halfchain[i]->patch)) {
                set_halfchain_pn_from_patch(e.halfchain[i]);
                break;
            }
        }
    }
}
void SketchRetopo::set_halfchain_pn_from_patch(curvenetwork::Halfchain* c) {
    auto patch = c->patch;
    if (!curvenetwork::Patch::is_ordinary(patch))
        return;
    
    // find patch boundary halfedge corresponding to this halfchain
    auto patch_h_start = patch->prev_halfedge_handle(patch->halfedge_handle(patch->patch_corner(0)));
    auto patch_h = patch_h_start;
    while (true) {
        if (patch->data(patch_h).halfchain == c)
            break;
        
        patch_h = patch->prev_halfedge_handle(patch_h);
        if (patch_h == patch_h_start)
            // this should never happen!
            return;
    }
    
    // sample patch boundary points
    Polyline_PointNormal patch_boundary_curve;
    patch_boundary_curve.push_back(patch->data(patch->to_vertex_handle(patch_h)).pn);
    for (int i = 0; i < c->num_subdiv(); ++i) {
        patch_boundary_curve.push_back(patch->data(patch->from_vertex_handle(patch_h)).pn);
        patch_h = patch->prev_halfedge_handle(patch_h);
    }
    
    // update curve network vertex positions
    int n = c->toPolyline().size() - 1;
    auto h = c->halfedge_front;
    h->from_vertex()->pn = patch_boundary_curve.front();
    for (int i = 1; i <= n; ++i) {
        double t = i / static_cast<double>(n);
        h->vertex->pn = patch_boundary_curve.point_at(t);
        project(h->vertex->pn);
        
        h = h->next;
    }
}

void SketchRetopo::snakes_move(curvenetwork::Edgechain* e) {
    auto c = e->halfchain[0];
    
    vector<PointNormal> result;
    
    for (auto h = c->halfedge_front->next; h->next != c->halfedge_back; h = h->next) {
        auto pn_m2 = h->prev->from_vertex()->pn;
        auto pn_m1 = h->prev->vertex->pn;
        auto pn_0  = h->vertex->pn;
        auto pn_p1 = h->next->vertex->pn;
        auto pn_p2 = h->next->next->vertex->pn;
        
        // internal energy term
        Vector3d gradient_i1 = ((pn_m1 + pn_p1) / 2.0                          - pn_0).head(3);
        Vector3d gradient_i2 = ((-pn_m2 + 4 * pn_m1 + 4 * pn_p1 - pn_p2) / 6.0 - pn_0).head(3);
        
        // external energy term
        Vector3d gradient_e = Vector3d::Zero();
        auto hit = intersect(pn_0);
        if (hit) {
            auto f = basemesh.face_handle(hit.id0);
            auto fv = basemesh.fv_iter(f);
            for (int i = 0; i < 3; ++i, ++fv)
                gradient_e += static_cast<double>(bc(hit, i)) * basemesh.data(*fv).snake_energy_gradient;
        }
        
        Vector3d gradient =
            ( configTemp.snakes_internal1 * gradient_i1
            + configTemp.snakes_internal2 * gradient_i2
            + configTemp.snakes_external  * gradient_e)
            / (configTemp.snakes_internal1 + configTemp.snakes_internal2 + configTemp.snakes_external);
        
        pn_0.head(3) += configTemp.snakes_damping * gradient;
        project(pn_0);
        result.push_back(pn_0);
    }

    auto i = result.begin();
    for (auto h = c->halfedge_front->next; h->next != c->halfedge_back; h = h->next, ++i) {
        h->vertex->pn = *i;
    }
}
void SketchRetopo::snakes_move(curvenetwork::Vertex   * v) {
    auto hit = intersect(v->pn);
    if (!hit)
        return;
    
    Vector3d gradient_e = Vector3d::Zero();
    auto f = basemesh.face_handle(hit.id0);
    auto fv = basemesh.fv_iter(f);
    for (int i = 0; i < 3; ++i, ++fv)
        gradient_e += static_cast<double>(bc(hit, i)) * basemesh.data(*fv).snake_energy_gradient;
    
    Vector3d offset_value =
        configTemp.snakes_damping *
        configTemp.snakes_external / (configTemp.snakes_internal1 + configTemp.snakes_internal2 + configTemp.snakes_external) *
        gradient_e;
    
    auto offset_func = [&] (PointNormal& pn, double w) {
        pn.head(3) += w * offset_value;
        project(pn);
    };
    
    // apply offset to selected vertex and its neighbors
    offset_func(v->pn, 1);
    for (curvenetwork::VOCIter c(v); c; ++c) {
        int n = c->toPolyline().size() - 1;
        auto h = c->halfedge_front;
        for (int i = 1; i < n; ++i, h = h->next) {
            double w = 1.0 - i / static_cast<double>(n);
            offset_func(h->vertex->pn, w);
        }
    }
}

void SketchRetopo::toggle_symmetric() {
    if (!curvenetwork.empty() && !zenity_util::question("SketchRetopo", "All current curves will be lost. Continue?"))
        return;
    
    curvenetwork.clear();
    memento.init();
    configTemp.autoSave.unsaved = true;
    
    configSaved.symmetric = !configSaved.symmetric;
    
    if (configSaved.symmetric) {
        vector<int> edge_visited(basemesh.n_edges(), 0);
        
        auto is_crossing = [this] (const BaseMesh::HHandle h) -> bool { auto p0p1 = basemesh.util_halfedge_to_point_pair(h); return p0p1.first[0] * p0p1.second[0] <  0;                                                };
        auto is_touching = [this] (const BaseMesh::HHandle h) -> bool { auto p0p1 = basemesh.util_halfedge_to_point_pair(h); return p0p1.first[0] * p0p1.second[0] == 0 && (p0p1.first[0] != 0 || p0p1.second[0] != 0); };
        auto is_coincide = [this] (const BaseMesh::HHandle h) -> bool { auto p0p1 = basemesh.util_halfedge_to_point_pair(h); return                                         p0p1.first[0] == 0 && p0p1.second[0] == 0 ; };
        
        for (auto h_iter : basemesh.halfedges()) {
            int e_idx = basemesh.edge_handle(h_iter).idx();
            if (edge_visited[e_idx]) continue;
            edge_visited[e_idx] = 1;
            
            if (!is_crossing(h_iter) && !is_coincide(h_iter))
                continue;
            
            BaseMesh::HHandle h_start = h_iter;
            bool is_loop = false;
            
            for (auto h_current = h_start; ; ) {
                BaseMesh::HHandle h_next;
                
                if (is_coincide(h_current)) {     // h_current is coinciding with x=0
                    auto v = basemesh.to_vertex_handle(h_current);
                    for (auto voh : basemesh.voh_range(v)) {
                        if (basemesh.opposite_halfedge_handle(voh) == h_current)
                            continue;
                        
                        if (is_coincide(voh)) {
                            h_next = voh;
                            break;
                        }
                        
                        auto h2 = basemesh.next_halfedge_handle(voh);
                        if (is_crossing(h2)) {
                            h_next = h2;
                            break;
                        }
                    }
                    
                } else {                                    // h_current is crossing x=0
                    auto h_opposite = basemesh.opposite_halfedge_handle(h_current);
                    if (basemesh.is_boundary(h_opposite)) {
                        h_start = h_opposite;
                        break;
                    }
                    
                    auto h0 = basemesh.next_halfedge_handle(basemesh.opposite_halfedge_handle(h_current));
                    auto h1 = basemesh.next_halfedge_handle(h0);
                    if (is_touching(h0)) {
                        auto v = basemesh.to_vertex_handle(h0);
                        for (auto voh : basemesh.voh_range(v)) {
                            if (is_coincide(voh)) {
                                h_next = voh;
                                break;
                            }
                            
                            auto h2 = basemesh.next_halfedge_handle(voh);
                            if (basemesh.opposite_halfedge_handle(h2) == h_current)
                                continue;
                            
                            if (is_crossing(h2)) {
                                h_next = h2;
                                break;
                            }
                        }
                        
                    } else {                // either h0 or h1 is crossing x=0
                        h_next = is_crossing(h0) ? h0 : h1;
                    }
                }
                
                if (h_next == h_start) {
                    // loop detected
                    is_loop = true;
                    break;
                }
                
                if (!h_next.is_valid()) {
                    // open end detected
                    h_start = basemesh.opposite_halfedge_handle(h_current);
                    break;
                }
                
                h_current = h_next;
            }
            
            // helper functors for adding new curve network vertex
            auto vertex_by_v = [this] (BaseMesh::VHandle basemesh_v) {
                auto pn = basemesh.get_pn(basemesh_v);
                pn[0] = pn[3] = 0;
                pn_normalize(pn);
                pn.head(3) += configSaved.projOffset * pn.tail(3);
                auto cn_v = curvenetwork.new_vertex();
                cn_v->pn = pn;
                return cn_v;
            };
            auto vertex_by_h = [this] (BaseMesh::HHandle basemesh_h, double t) {
                auto v0v1 = basemesh.util_halfedge_to_vertex_pair(basemesh_h);
                auto pn0 = basemesh.get_pn(v0v1.first );
                auto pn1 = basemesh.get_pn(v0v1.second);
                PointNormal pn = (1 - t) * pn0 + t * pn1;
                pn[0] = pn[3] = 0;
                pn_normalize(pn);
                pn.head(3) += configSaved.projOffset * pn.tail(3);
                auto cn_v = curvenetwork.new_vertex();
                cn_v->pn = pn;
                return cn_v;
            };
            
            // do the exactly same walking again, this time adding curve network vertices
            vector<curvenetwork::Vertex*> vertices;
            
            if (!is_loop) {
                // preprocess for boundary case
                if (is_coincide(h_start)) {
                    vertices.push_back(vertex_by_v(basemesh.from_vertex_handle(h_start)));
                    
                }else if (!basemesh.is_boundary(h_start)) {
                    auto h = basemesh.next_halfedge_handle(h_start);
                    assert(is_touching(h));
                    vertices.push_back(vertex_by_v(basemesh.to_vertex_handle(h)));
                }
            }
            
            for (auto h_current = h_start; ; ) {
                edge_visited[basemesh.edge_handle(h_current).idx()] = 1;
                
                BaseMesh::HHandle h_next;
                
                if (is_coincide(h_current)) {     // h_current is coinciding with x=0
                    auto v = basemesh.to_vertex_handle(h_current);
                    
                    vertices.push_back(vertex_by_v(v));
                    
                    for (auto voh : basemesh.voh_range(v)) {
                        if (basemesh.opposite_halfedge_handle(voh) == h_current)
                            continue;
                        
                        if (is_coincide(voh)) {
                            h_next = voh;
                            break;
                        }
                        
                        auto h2 = basemesh.next_halfedge_handle(voh);
                        if (is_crossing(h2)) {
                            h_next = h2;
                            break;
                        }
                    }
                    
                } else {                                    // h_current is crossing x=0
                    auto p0p1 = basemesh.util_halfedge_to_point_pair(h_current);
                    double t = -p0p1.first[0] / (p0p1.second[0] - p0p1.first[0]);
                    
                    vertices.push_back(vertex_by_h(h_current, t));
                    
                    auto h_opposite = basemesh.opposite_halfedge_handle(h_current);
                    if (basemesh.is_boundary(h_opposite))
                        // reached other end of open cross section
                        break;
                    
                    auto h0 = basemesh.next_halfedge_handle(basemesh.opposite_halfedge_handle(h_current));
                    auto h1 = basemesh.next_halfedge_handle(h0);
                    if (is_touching(h0)) {
                        auto v = basemesh.to_vertex_handle(h0);
                        for (auto voh : basemesh.voh_range(v)) {
                            if (is_coincide(voh)) {
                                h_next = voh;
                                break;
                            }
                            
                            auto h2 = basemesh.next_halfedge_handle(voh);
                            if (h2 == h_opposite)
                                continue;
                            
                            if (is_crossing(h2)) {
                                h_next = h2;
                                break;
                            }
                        }
                        
                    } else {                // either h0 or h1 is crossing x=0
                        h_next = is_crossing(h0) ? h0 : h1;
                    }
                }
                
                if (!h_next.is_valid()) {
                    // postprocess for boundary case
                    if (!is_coincide(h_current) && !basemesh.is_boundary(h_current)) {
                        auto h = basemesh.next_halfedge_handle(h_current);
                        assert(is_touching(h));
                        vertices.push_back(vertex_by_v(basemesh.to_vertex_handle(h)));
                    }
                    break;
                }
                
                if (h_next == h_start)
                    // end of loop
                    break;
                
                h_current = h_next;
            }
            
            // insert more between each pair of vertices
            auto vertices_temp = vertices;
            auto insert_pos = vertices_temp.begin();
            for (auto i : adjacent_pairs(vertices, is_loop)) {
                int num_insertions = 9;
                auto d = (i.second->pn - i.first->pn) / (num_insertions + 1.0);
                PointNormal pn = i.first->pn + d;
                ++insert_pos;
                for (int j = 0; j < num_insertions; ++j, ++insert_pos, pn += d) {
                    auto v = curvenetwork.new_vertex();
                    v->pn = pn_normalized(pn);
                    insert_pos = vertices_temp.insert(insert_pos, &*v);
                }
            }
            vertices = vertices_temp;
            
            // reverse vertices if h_start's right side corresponds to negative x
            bool do_reverse = false;
            {
                auto pn0 = vertices[0]->pn;
                auto pn1 = vertices[1]->pn;
                Vector3d d = (pn1 - pn0).head(3);
                Vector3d n = pn0.tail(3);
                do_reverse = d.cross(n).x() < 0;
            }
            if (do_reverse)
                reverse(vertices.begin(), vertices.end());
            
            // add halfedge without next/prev info
            vector<curvenetwork::Halfedge*> halfedges_forward;
            vector<curvenetwork::Halfedge*> halfedges_backward;
            for (auto i : adjacent_pairs(vertices, is_loop)) {
                auto v0 = i.first;
                auto v1 = i.second;
                
                auto h01 = curvenetwork.new_halfedge();
                auto h10 = curvenetwork.new_halfedge();
                curvenetwork::set_vertex_halfedge(v0, v1, h01, h10);
                halfedges_forward .push_back(h01);
                halfedges_backward.push_back(h10);
                    
                h01->imaginary_patch = &curvenetwork::Patch::imaginary_patch_symmetry;
            }
            
            // set halfedge prev/next info
            int n_halfedges = halfedges_forward.size();
            for (int i = 0; i < n_halfedges; ++i) {
                auto h_forward  = halfedges_forward [i];
                auto h_backward = halfedges_backward[i];
                
                h_forward ->next = !is_loop && i == n_halfedges - 1 ? nullptr : halfedges_forward [(i + 1) % n_halfedges];
                h_forward ->prev = !is_loop && i == 0 ? nullptr : halfedges_forward [(i + n_halfedges - 1) % n_halfedges];
                h_backward->prev = !is_loop && i == n_halfedges - 1 ? nullptr : halfedges_backward[(i + 1) % n_halfedges];
                h_backward->next = !is_loop && i == 0 ? nullptr : halfedges_backward[(i + n_halfedges - 1) % n_halfedges];
            }
        }
    }
    
    trace_basemesh_boundary();
}

void SketchRetopo::trace_basemesh_boundary() {
    vector<int> halfedge_visited(basemesh.n_halfedges(), 0);
    
    for (auto h : basemesh.halfedges()) {
        if (!basemesh.is_boundary(h) || halfedge_visited[h.idx()])
            continue;
        
        auto h_start = h;
        Polyline_PointNormal curve;
        for (auto i = h_start; ; ) {
            halfedge_visited[i.idx()] = 1;
            
            auto pn = basemesh.get_pn(basemesh.from_vertex_handle(i));
            
            curve.push_back(pn);
            
            i = basemesh.next_halfedge_handle(i);
            if (i == h_start)
                break;
        }
        
        bool is_symmetric_open = false;
        if (configSaved.symmetric) {
            // don't add boundary loop that is entirely on the negative side
            bool is_all_negative = boost::count_if(curve, [] (PointNormal& pn) { return pn.x() >= 0; }) == 0;
            if (is_all_negative)
                continue;
            
            // if there exists a segment having both positive and negative x coordinates, move one endpoint to x=0
            for (int i = 0; i < curve.size(); ++i) {
                int i_next = (i + 1) % curve.size();
                auto& pn0 = curve[i];
                auto& pn1 = curve[i_next];
                if (pn0.x() * pn1.x() >= 0)
                    continue;
                
                double t = -pn0.x() / (pn1.x() - pn0.x());
                PointNormal pn = (1 - t) * pn0 + t * pn1;
                pn[0] = pn[3] = 0;
                pn_normalize(pn);
                pn0 = pn;
            }
            
            // look for curve point with negative x adjacent to x=0
            int rotate_pos = -1;
            for (int i = 0; i < curve.size(); ++i) {
                int i_prev = i == 0 ? curve.size() - 1 : i - 1;
                if (curve[i_prev].x() < 0 && curve[i].x() >= 0) {
                    assert(curve[i].x() == 0);
                    rotate_pos = i;
                    break;
                }
            }
            
            // erase all points with negative x
            if (rotate_pos != -1) {
                rotate(curve.begin(), curve.begin() + rotate_pos, curve.end());
                boost::remove_if(curve, [] (PointNormal& pn) { return pn.x() < 0; });
                is_symmetric_open = true;
            }
        }
        
        // projOffset
        for (auto& pn : curve)
            pn.head(3) += configSaved.projOffset * pn.tail(3);
        
        // add vertex
        vector<curvenetwork::Vertex*> vertices;
        for (int i = 0; i < curve.size(); ++i) {
            auto& pn = curve[i];
            
            if (is_symmetric_open && (i == 0 || i == curve.size() - 1)) {
                MinSelector<curvenetwork::Vertex*> v_min;
                for (auto& v : curvenetwork.vertices) {
                    double dist = (pn - v.pn).norm();
                    v_min.update(dist, &v);
                }
                assert(v_min.score < 0.00001);         // should be zero
                vertices.push_back(v_min.value);
                
            } else {
                auto v = curvenetwork.new_vertex();
                v->pn = pn;
                vertices.push_back(v);
            }
        }
        
        // add halfedge without next/prev info
        vector<curvenetwork::Halfedge*> halfedges_forward;
        vector<curvenetwork::Halfedge*> halfedges_backward;
        for (int i = 0; i < curve.size(); ++i) {
            if (is_symmetric_open && i == curve.size() - 1)
                continue;
            
            curvenetwork::Vertex* v0 = vertices[i];
            curvenetwork::Vertex* v1 = vertices[i < curve.size() - 1 ? i + 1 : 0];
            
            auto h01 = curvenetwork.new_halfedge();
            auto h10 = curvenetwork.new_halfedge();
            
            if (is_symmetric_open && i == 0) {
                curvenetwork::set_next_prev(v0->halfedge->opposite, h01);
                curvenetwork::set_next_prev(h10, v0->halfedge);
                h01->prev->is_corner = true;
                h10      ->is_corner = true;
                
            } else if (is_symmetric_open && i == curve.size() - 2) {
                curvenetwork::set_next_prev(h01, v1->halfedge);
                curvenetwork::set_next_prev(v1->halfedge->opposite, h10);
                h01      ->is_corner = true;
                h10->prev->is_corner = true;
            }
            
            curvenetwork::set_vertex_halfedge(v0, v1, h01, h10);
            
            halfedges_forward .push_back(h01);
            halfedges_backward.push_back(h10);
            
            h01->imaginary_patch = &curvenetwork::Patch::imaginary_patch_boundary;
        }
        
        // set halfedge prev/next info
        int n_halfedges = halfedges_forward.size();
        for (int i = 0; i < n_halfedges; ++i) {
            auto h_forward  = halfedges_forward [i];
            auto h_backward = halfedges_backward[i];
            
            if (!(is_symmetric_open && i == n_halfedges - 1)) {
                h_forward ->next = halfedges_forward [(i + 1) % n_halfedges];
                h_backward->prev = halfedges_backward[(i + 1) % n_halfedges];
            }
            
            if (!(is_symmetric_open && i == 0)) {
                h_forward ->prev = halfedges_forward [(i + n_halfedges - 1) % n_halfedges];
                h_backward->next = halfedges_backward[(i + n_halfedges - 1) % n_halfedges];
            }
        }
    }
    
    curvenetwork.generate_halfchains();
    generate_patches();
}

Polyline_PointNormal SketchRetopo::geodesic_compute(const PointNormal& pn0, const PointNormal& pn1) {
    auto hit0 = intersect(pn0);
    auto hit1 = intersect(pn1);
    if (!hit0 || !hit1)
        return Polyline_PointNormal();
    
    geodesic::SurfacePoint source(&geodesic_mesh.faces()[hit0.id0], 1 - hit0.u - hit0.v, hit0.u, hit0.v);
    geodesic::SurfacePoint target(&geodesic_mesh.faces()[hit1.id0], 1 - hit1.u - hit1.v, hit1.u, hit1.v);
    
    std::vector<geodesic::SurfacePoint> path;
    geodesic_algorithm.geodesic(source, target, path);
    
    Polyline_PointNormal result;
    for (auto& sp : path) {
        auto sp_type = sp.type();
        if (sp_type == geodesic::UNDEFINED_POINT)
            continue;
        
        PointNormal pn = PointNormal::Zero();
        pn.head(3) << sp.x(), sp.y(), sp.z();
        if (sp_type == geodesic::FACE) {
            for (int i = 0; i < 3; ++i) {
                auto v = basemesh.vertex_handle(sp.base_element()->adjacent_vertices()[i]->id());
                pn.tail(3) += o2e(basemesh.normal(v));
            }
        
        } else if (sp_type == geodesic::EDGE) {
            for (int i = 0; i < 2; ++i) {
                auto v = basemesh.vertex_handle(sp.base_element()->adjacent_vertices()[i]->id());
                pn.tail(3) += o2e(basemesh.normal(v));
            }
        
        } else if (sp.type() == geodesic::VERTEX) {
            auto v = basemesh.vertex_handle(sp.base_element()->id());
            pn.tail(3) = o2e(basemesh.normal(v));
        }
        pn_normalize(pn);
        
        result.push_back(pn);
    }
    return result;
}
