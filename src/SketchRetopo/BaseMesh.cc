#include "SketchRetopo.hh"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <GL/glew.h>
#include <kt84/Clock.hh>
#include <kt84/util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;

namespace {
    auto& core = SketchRetopo::get_instance();

//#ifdef _WIN32
//    inline FILE* popen(const char* command, const char* mode) { return _popen(command, mode); }
//    inline int pclose(FILE* file) { return _pclose(file); }
//#endif
    
    Vector3d scalar_to_color(double t) {
        double r = std::max<double>(0, 2 * t - 1);
        double g = 1 - 2 * std::abs<double>(t - 0.5);
        double b = std::max<double>(0, 1 - 2 * t);
        return Vector3d(r, g, b);
    }
}

void BaseMesh::init() {
    update_normals();
    
#ifndef NDEBUG
    debugInfo_get(true, true, true);
#endif
    
    cotanWeight_compute();
    centerOfMass_compute();
    boundingBox_compute();
    
    compute_k1k2();
    normalCurvature_compute();
    normalVariation_compute();
    compute_snake_energy();
    
    invalidate_displist();
}

void BaseMesh::smooth_normal() {
    update_normals();
    
    for (auto v : vertices()) {
        auto& vdata = data(v).laplaceIterative;
        vdata.value    = o2e(normal(v));
        vdata.is_fixed = false;
    }
    
    for (int i = 0; i < normalSmoothIter; ++i) {
        laplaceIterative_compute(1);
        
        for (auto v : vertices())
            data(v).laplaceIterative.value.normalize();
    }
}

void BaseMesh::smooth_expmap() {
    // fix unprocessed and their 1-ring vertices
    for (auto v : vertices()) {
        auto& vdata = data(v);
        
        vdata.laplaceIterative.value.head(2) = vdata.expmap.uv;
        
        if (!vdata.expmap.paramed) {
            vdata.laplaceIterative.is_fixed = true;
            continue;
        }
        
        vdata.laplaceIterative.is_fixed = false;
        for (auto w : vv_range(v)) {
            if (data(w).expmap.paramed) continue;
            vdata.laplaceIterative.is_fixed = true;
            break;
        }
    }
    
    // smooth expmap.uv using Laplacian
    laplaceIterative_compute(1);
    
    // copy result
    for (auto v : vertices()) {
        auto& vdata = data(v);
        if (vdata.laplaceIterative.is_fixed) continue;
        
        vdata.expmap.uv = vdata.laplaceIterative.value.head(2);
    }
    
    invalidate_displist();
}
void BaseMesh::compute_k1k2() {
    // TODO: compute principal curvature
    for (auto v : vertices()) {
        auto& k1 = data(v).k1k2[0];
        auto& k2 = data(v).k1k2[1];
        
        k1 = k2 = 0;
    }
    
    //// output current mesh to temporary file in OFF format
    //{
    //    ClkSimple clk("OpenMesh::IO::write_mesh");
    //    OpenMesh::IO::write_mesh(*this, "temp_mesh.off");
    //}
    //
    //// run CurvatureX3. http://people.cs.umass.edu/~kalo/papers/curvature
    //const char* command = "CurvatureX3 -n temp_mesh.off temp_curvature.txt";
    //cout << command << endl;
    //FILE* p = popen(command, "r");
    //char buf[4096];     // get result
    //while (fgets(buf, 4096, p) != NULL) {}
    //string result = buf;
    //cout << buf << endl;   // print output
    //pclose(p);
    //
    //// read from output file
    //ifstream fin("temp_curvature.txt");
    //
    //string dummy;
    //for (int i = 0; i < 8; ++i) fin >> dummy;
    //
    //for (auto v : vertices()) {
    //    auto& k1 = data(v).k1k2[0];
    //    auto& k2 = data(v).k1k2[1];
    //    
    //    double k1v_x, k1v_y, k1v_z;
    //    double k2v_x, k2v_y, k2v_z;
    //    
    //    fin >> k1 >> k2 >> k1v_x >> k1v_y >> k1v_z >> k2v_x >> k2v_y >> k2v_z;
    //}
}
void BaseMesh::compute_snake_energy() {
    //---------------------------------+
    // define per-vertex feature value |
    //---------------------------------+
    for (auto v : vertices()) {
        auto& vdata = data(v);
        vdata.snake_feature =
            snake_feature_type == SnakeFeatureType::MEAN_CURVATURE   ? 0.5 * (vdata.k1k2[0] + vdata.k1k2[1]) :    // mean curvature
            snake_feature_type == SnakeFeatureType::NORMAL_VARIATION ? vdata.normalVariation                 :    // simple normal variation (as in the original geometric snakes paper)
            snake_feature_type == SnakeFeatureType::NORMAL_CURVATURE ? vdata.normalCurvature                 :    // simple normal curvature
            0;
    }
    
    //---------------------------------+
    // smooth feature value (optional) |
    //---------------------------------+
    for (auto v : vertices()) {
        auto& vdata = data(v);
        vdata.laplaceIterative.is_fixed = false;
        vdata.laplaceIterative.value[0] = vdata.snake_feature;
    }
    laplaceIterative_compute(featureSmoothIter);
    for (auto v : vertices()) {
        auto& vdata = data(v);
        vdata.snake_feature = vdata.laplaceIterative.value[0];
    }
    
    //--------------------------+
    // define per-vertex energy |
    //--------------------------+
    if (use_feature_as_energy) {
        // simply copied from feature
        snake_energy_max = -util::dbl_max();
        snake_energy_min =  util::dbl_max();
        for (auto v : vertices()) {
            auto& vdata = data(v);
            vdata.snake_energy = vdata.snake_feature;
            
            if (negate_energy)
                vdata.snake_energy *= -1;
            
            snake_energy_max = max<double>(snake_energy_max, vdata.snake_energy);
            snake_energy_min = min<double>(snake_energy_min, vdata.snake_energy);
        }
        
    } else {
        // compute per-face feature gradient
        for (auto v : vertices())
            data(v).gradient_input = data(v).snake_feature;
        gradient_compute();
        
        // compute average of magniutude of per-face feature gradient of one-ring faces
        snake_energy_max = 0;
        snake_energy_min = util::dbl_max();
        for (auto v : vertices()) {
            auto& snake_energy = data(v).snake_energy;
            snake_energy = 0;
            
            for (auto f : vf_range(v)) {
                auto f_feature_gradient = data(f).gradient_output;
                snake_energy += f_feature_gradient.norm();
            }
            snake_energy /= valence(v);
            
            if (negate_energy)
                snake_energy *= -1;
            
            snake_energy_max = max<double>(snake_energy_max, snake_energy);
            snake_energy_min = min<double>(snake_energy_min, snake_energy);
        }
    }
    
    //----------------------------------+
    // compute per-face energy gradient |
    //----------------------------------+
    for (auto v : vertices())
        data(v).gradient_input = data(v).snake_energy;
    gradient_compute();
    
    //---------------------------------------------------------------------------+
    // compute per-vertex energy gradient as average of per-face energy gradient |
    //---------------------------------------------------------------------------+
    for (auto v : vertices()) {
        auto& snake_energy_gradient = data(v).snake_energy_gradient;
        snake_energy_gradient.setZero();
        
        for (auto f : vf_range(v)) {
            auto f_energy_gradient = data(f).gradient_output;
            snake_energy_gradient += f_energy_gradient;
        }
        snake_energy_gradient /= valence(v);
    }
    
    invalidate_displist();
}
void BaseMesh::render_fill() {
    displist_fill.render([&] () {
        glBegin(GL_TRIANGLES);
        for (auto f : faces()) {
            if (data(f).is_hidden) continue;
            bool is_expmap_paramed = true;
            for (auto v : fv_range(f)) {
                if (data(v).expmap.paramed) continue;
                is_expmap_paramed = false;
                break;
            }
            for (auto h : fh_range(f)) {
                auto v = to_vertex_handle(h);
                auto p  = point (v);
                auto n  = normal(v);
                auto expmap_uv  = data(v).expmap.uv;
                auto overlay_uv = o2e(texcoord2D(h));
                double snake_energy = data(v).snake_energy;
                double scalar = (snake_energy - snake_energy_min) / (snake_energy_max - snake_energy_min);
                scalar = atan(core.configRender.snake_energy_scale * scalar) / (util::pi() * 0.5);
                glColor3d(scalar, scalar, scalar);
                glTexCoord2dv(expmap_uv.data());
                glMultiTexCoord2dv(GL_TEXTURE1, overlay_uv.data());
                glNormal3dv  (n.data());
                glVertex3dv  (p.data());
            }
        }
        glEnd();
    });
}
void BaseMesh::render_line() {
    displist_line.render([&] () {
        glBegin(GL_LINES);
        for (auto e : edges()) {
            for (int i = 0; i < 2; ++i) {
                auto v = to_vertex_handle(halfedge_handle(e, i));
                auto p  = point (v);
                auto n  = normal(v);
                glNormal3dv(&n[0]);
                glVertex3dv(&p[0]);
            }
        }
        glEnd();
    });
}
PointNormal BaseMesh::get_pn(VHandle v) const {
    PointNormal pn;
    pn << o2e(point(v)), o2e(normal(v));
    return pn;
}
