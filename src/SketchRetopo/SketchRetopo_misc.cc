#include "SketchRetopo.hh"
#include "curvenetwork/Circulator.hh"
#include <sstream>
#include <base64/decode.h>
#include <lodepng/lodepng.h>
#include <kt84/math/RBF.hh>
#include <kt84/adjacent_pairs.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;

void SketchRetopo::state_set(EnumState enumState) {
    State* state_new = nullptr;
    if      (enumState == EnumState::Sketch      ) state_new = &stateSketch      ;
    else if (enumState == EnumState::Spine       ) state_new = &stateSpine       ;
    else if (enumState == EnumState::Autocmpl    ) state_new = &stateAutocmpl    ;
    else if (enumState == EnumState::Laser       ) state_new = &stateLaser       ;
    else if (enumState == EnumState::Cylinder    ) state_new = &stateCylinder    ;
    else if (enumState == EnumState::DeformCurve ) state_new = &stateDeformCurve ;
    else if (enumState == EnumState::EditCorner  ) state_new = &stateEditCorner  ;
    else if (enumState == EnumState::EditTopology) state_new = &stateEditTopology;
    else if (enumState == EnumState::EdgeLoop    ) state_new = &stateEdgeLoop    ;
    else if (enumState == EnumState::MoveVertex  ) state_new = &stateMoveVertex  ;
    
    if (!state_new || state == state_new) return;
    
    state = state_new;
    state->init();
}
SketchRetopo::EnumState SketchRetopo::state_get() const {
    if      (state == &stateSketch      ) return EnumState::Sketch      ;
    else if (state == &stateSpine       ) return EnumState::Spine       ;
    else if (state == &stateAutocmpl    ) return EnumState::Autocmpl    ;
    else if (state == &stateLaser       ) return EnumState::Laser       ;
    else if (state == &stateCylinder    ) return EnumState::Cylinder    ;
    else if (state == &stateDeformCurve ) return EnumState::DeformCurve ;
    else if (state == &stateEditCorner  ) return EnumState::EditCorner  ;
    else if (state == &stateEditTopology) return EnumState::EditTopology;
    else if (state == &stateEdgeLoop    ) return EnumState::EdgeLoop    ;
    else if (state == &stateMoveVertex  ) return EnumState::MoveVertex  ;
    assert(false);
    return static_cast<EnumState>(-1);
}
void SketchRetopo::memento_store() {
    memento.store(curvenetwork);
    configTemp.autoSave.unsaved = true;
}
void SketchRetopo::memento_undo() {
    if (memento.undo(curvenetwork)) {
        curvenetwork.validate_all_ptr();
        for (auto& p : curvenetwork.patches)
            p.invalidate_displist();
        curvenetwork.invalidate_displist();
        state->init();
        configTemp.autoSave.unsaved = true;
    }
}
void SketchRetopo::memento_redo() {
    if (memento.redo(curvenetwork)) {
        curvenetwork.validate_all_ptr();
        for (auto& p : curvenetwork.patches)
            p.invalidate_displist();
        curvenetwork.invalidate_displist();
        state->init();
        configTemp.autoSave.unsaved = true;
    }
}
void SketchRetopo::material_texture_update() {
    // decode PNG from base64 string
    stringstream ss_encoded, ss_decoded;
    ss_encoded.str(material_texture_png_base64);
    base64::decoder D;
	D.decode(ss_encoded, ss_decoded);
    
    // read PNG using awesome lodepng library!
    vector<unsigned char> pixels;   // pixels in RGB
    unsigned int width, height;
    auto str_decoded = ss_decoded.str();
    lodepng::decode(pixels, width, height, reinterpret_cast<unsigned char*>(&str_decoded[0]), str_decoded.size(), LCT_RGB);
    
    // send pixel data to GPU
    material_texture.bind();
    material_texture.allocate(width, height, GL_RGB, GL_RGB);
    material_texture.copy_cpu2gpu(GL_UNSIGNED_BYTE, &pixels[0]);
    material_texture.unbind();
}
void SketchRetopo::hide_mesh_by_stroke() {
    if (hide_stroke.size() < 3) {
        // reset hiding
        for (auto  f : basemesh.faces())      basemesh.data(f).is_hidden = false;
        for (auto& v : curvenetwork.vertices) v               .is_hidden = false;
        for (auto& p : curvenetwork.patches)  p               .is_hidden = false;
    
    } else {
        // resample with certain length
        double s = (camera->width + camera->height) * 0.001;
        hide_stroke.resample_by_length(s);
        
        // construct screen-space scalar field using Radial Basis Function
        typedef RBF<2, 1, RBFKernel_SquaredLog, 1> RBF;
        RBF rbf;
        for (auto i : adjacent_pairs(hide_stroke, false)) {
            Vector2d p = (i.first + i.second) * 0.5;
            Vector2d d = eigen_util::rotate90(i.second - i.first).normalized() * s;
            RBF::Value value;
            value <<  s;    rbf.add_constraint(p + d, value);
            value << -s;    rbf.add_constraint(p - d, value);
        }
        rbf.factorize_and_solve();
        
        // sample rbf at various places
#pragma omp parallel for
        for (int i = 0; i < basemesh.n_faces(); ++i) {
            auto f = basemesh.face_handle(i);
            Vector3d world_pos  = o2e(basemesh.util_face_center(f));
            Vector2d screen_pos = graphics_util::project(world_pos).head(2);
            basemesh.data(f).is_hidden = rbf(screen_pos)[0] < 0;
        }
        for (auto& v : curvenetwork.vertices) {
            Vector3d world_pos  = v.pn.head(3);
            Vector2d screen_pos = graphics_util::project(world_pos).head(2);
            v.is_hidden = rbf(screen_pos)[0] < -1.5 * s;
        }
        for (auto& p : curvenetwork.patches) {
            p.is_hidden = true;
            for (curvenetwork::PCIter c(&p); c; ++c) {
                for (curvenetwork::CHIter h(&*c); h; ++h) {
                    if (!h->vertex->is_hidden) {
                        p.is_hidden = false;
                        break;
                    }
                    if (!p.is_hidden) break;
                }
                if (!p.is_hidden) break;
            }
        }
    }
    
    embree_init();
    curvenetwork.invalidate_displist();
    basemesh.invalidate_displist();
}
