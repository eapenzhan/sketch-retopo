#include "decl.hh"
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <patchgen/decl.hh>
#include <kt84/graphics/graphics_util.hh>
#include <kt84/tw_util.hh>
#include <kt84/glut_util.hh>
#include <kt84/util.hh>
#include <kt84/eigen_util.hh>
#include <kt84/geometry/CameraFree.hh>
#include <kt84/MinSelector.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

struct Globals {
    TwBar* bar = nullptr;
    CameraFree camera;
    
    demo::Patch patch;
    patchgen::PatchParam param;
    VectorXi l;
    
    enum struct Mode {
        ChangeSubdiv,
        AdjustParam,
        MoveVertex,
    } mode = Mode::ChangeSubdiv;
    
    int selected_side = 0;
    int selected_variable = 0;
    demo::Patch::VHandle selected_vertex;
    
    int max_subdiv = 20;
    bool hide_auxiliary = false;
    
    void init() {
        camera.eye << 0, 0, 5;
        l = Vector3i(2, 2, 2);
        demo::generate_patch(l, param, patch);
    }
} g;

void init_gl() {
    glClearColor(1, 1, 1, 1);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);    
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    g.camera.auto_flip_y = false;
}

// >>AntTweakBar setup-------------------------------------------------
void init_bar() {
    g.bar = TwNewBar("patchgen");
    TwDefine("patchgen fontsize=3 size='300 320' valueswidth=fit color='153 204 51' alpha=200 text=dark");
    
    tw_util::AddButton(g.bar, "save",
        [&](){
            demo::save_mesh(g.patch, demo::get_fname(g.param));
        }, "key=s");
    tw_util::AddVarCB<int>(g.bar, "num_sides", TW_TYPE_INT32,
        [&](const int& value){
            if (value == g.param.get_num_sides()) return;
            g.l.setConstant(value, 2);
            demo::generate_patch(g.l, g.param, g.patch);
            g.selected_side = 0;
            g.selected_variable = 0;
        },
        [&](int& value){
            value = g.param.get_num_sides();
        }, "min=2 max=6 keyincr=UP keydecr=DOWN");
    TwAddVarRW(g.bar, "mode", TwDefineEnum("mode_type", 0, 0), &g.mode, "enum='0{change subdiv}, 1{adjust param}, 2{move vertex}' key=SPACE");
    tw_util::AddButton(g.bar, "increment",
        [&](){
            if (g.mode == Globals::Mode::ChangeSubdiv) {
                ++g.l[g.selected_side];
                if (g.l.sum() % 2 != 0)
                    g.patch.clear();
                else
                    demo::generate_patch(g.l, g.param, g.patch);
            } else if (adjust_parameter(g.param, g.selected_variable, true)) {
                demo::generate_patch(g.param, g.patch);
            }
        }, "key=v");
    tw_util::AddButton(g.bar, "decrement",
        [&](){
            if (g.mode == Globals::Mode::ChangeSubdiv) {
                if (g.l[g.selected_side] < 2) return;
                --g.l[g.selected_side];
                if (g.l.sum() % 2 != 0)
                    g.patch.clear();
                else
                    demo::generate_patch(g.l, g.param, g.patch);
            } else if (patchgen::adjust_parameter(g.param, g.selected_variable, false)) {
                demo::generate_patch(g.param, g.patch);
            }
        }, "key=z");
    tw_util::AddButton(g.bar, "increment_by_2",
        [&](){
            if (g.mode == Globals::Mode::ChangeSubdiv) {
                g.l[g.selected_side] += 2;
                if (g.l.sum() % 2 != 0)
                    g.patch.clear();
                else
                    demo::generate_patch(g.l, g.param, g.patch);
            } else if (adjust_parameter(g.param, g.selected_variable, true)) {
                demo::generate_patch(g.param, g.patch);
            }
        }, "key=c");
    tw_util::AddButton(g.bar, "decrement_by_2",
        [&](){
            if (g.mode == Globals::Mode::ChangeSubdiv) {
                if (g.l[g.selected_side] < 3) return;
                g.l[g.selected_side] -= 2;
                if (g.l.sum() % 2 != 0)
                    g.patch.clear();
                else
                    demo::generate_patch(g.l, g.param, g.patch);
            } else if (adjust_parameter(g.param, g.selected_variable, false)) {
                demo::generate_patch(g.param, g.patch);
            }
        }, "key=x");
    tw_util::AddButton(g.bar, "switch_pattern_forward",
        [&](){
            if (patchgen::switch_pattern(g.param, true))
                demo::generate_patch(g.param, g.patch);
        }, "key=w");
    tw_util::AddButton(g.bar, "switch_pattern_backward",
        [&](){
            if (patchgen::switch_pattern(g.param, false))
                demo::generate_patch(g.param, g.patch);
        }, "key=W");
    
    // randomize
    TwAddVarRW(g.bar, "max_subdiv", TW_TYPE_INT32, &g.max_subdiv, "group=random min=1");
    tw_util::AddButton(g.bar, "randomize",
        [&]() {
            int num_sides = g.param.get_num_sides();
            g.l.setZero(num_sides);
            for (int i = 0; i < num_sides; ++i)
                g.l[i] = util::random_int(1, g.max_subdiv);
            if (g.l.sum() % 2 != 0) ++g.l[0];
            if (g.l.sum() < 4) g.l[0] += 2;
            
            demo::generate_patch(g.l, g.param, g.patch);
            
            g.selected_side = 0;
            g.selected_variable = 0;
        }, "group=random key=r");
    TwAddVarRW(g.bar, "hide_auxiliary", TW_TYPE_BOOLCPP, &g.hide_auxiliary, "key=h");
}
// <<AntTweakBar setup-------------------------------------------------

// >>GLUT callback functions---------------------------------------
void display_pre() {
    glut_util::defaultcb::display_pre();
    
    // set projection matrix
    double zNear = g.camera.eye.z() * 0.1;
    double zFar  = zNear * 100;
    double aspect_ratio = g.camera.width / static_cast<double>(g.camera.height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, aspect_ratio, zNear, zFar);
    
    // set modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(g.camera.eye, g.camera.center, g.camera.up);
}
void display_main() {
    auto& param = g.param;
    auto& patch = g.patch;
    
    int num_sides = param.get_num_sides();
    
    // draw patch
    ::glColor3d(0.7, 0.7, 0.7);
    patch.draw();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(2);
    ::glColor3d(0, 0, 0);
    patch.draw();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    if (g.hide_auxiliary) return;
    
    // draw boundary
    auto draw_boundary = [&] (GLenum mode) {
        for (int i = 0; i < num_sides; ++i) {
            glBegin(mode);
            for (int j = 0; j <= g.l[i]; ++j) {
                double t = i + j / static_cast<double>(g.l[i]);
                glVertex2d(patchgen::get_boundary_geometry(num_sides, t));
            }
            glEnd();
        }
    };
    glColor3d(0, 0, 0);
    glLineWidth(5);
    draw_boundary(GL_LINE_STRIP);
    glPointSize(12);
    draw_boundary(GL_POINTS);
    
    // draw selected boundary
    if (g.mode == Globals::Mode::ChangeSubdiv) {
        glColor3d(1, 0, 0);
        auto draw_selected_boundary = [&] (GLenum mode) {
            glBegin(mode);
            for (int j = 0; j <= g.l[g.selected_side]; ++j) {
                double t = g.selected_side + j / static_cast<double>(g.l[g.selected_side]);
                glVertex2d(patchgen::get_boundary_geometry(num_sides, t));
            }
            glEnd();
        };
        draw_selected_boundary(GL_LINE_STRIP);
        draw_selected_boundary(GL_POINTS);
    
    }
    
    // draw singularities
    patch.draw_singularities();
    
    // draw arrows indicating adjustable variables
    if (g.mode == Globals::Mode::AdjustParam) {
        auto draw_doublesided_arrow = [&] (Vector2d p0, Vector2d p1) {
            double arrow_size = 0.03 * g.camera.eye.z();
            Vector2d d = p1 - p0;
            p0 += 0.1 * d;
            p1 -= 0.1 * d;
            double r = d.norm();
            arrow_size = min<double>(arrow_size, r * 0.3);
            d *= arrow_size / r;
            Vector2d e1 = 0.2 * eigen_util::rotate90(d);
            glBegin(GL_QUADS);
            glVertex2d(p0 + d + e1);
            glVertex2d(p0 + d - e1);
            glVertex2d(p1 - d - e1);
            glVertex2d(p1 - d + e1);
            glEnd();
            Vector2d e2 = 0.5 * eigen_util::rotate90(d);
            glBegin(GL_TRIANGLES);
            glVertex2d(p0);    glVertex2d(p0 + d - e2);    glVertex2d(p0 + d + e2);
            glVertex2d(p1);    glVertex2d(p1 - d + e2);    glVertex2d(p1 - d - e2);
            glEnd();
        };
        
        if (!(num_sides == 3 && param.pattern_id == 0) &&
            !(num_sides == 4 && param.pattern_id == 0) &&
            !(num_sides == 5 && param.pattern_id == 0))
        {
            for (int i = 0; i < num_sides; ++i) {
                int j = param.permutation[i];
                Vector2d p0 = patchgen::get_boundary_geometry(num_sides, j + 0.49);
                Vector2d p1 = patchgen::get_boundary_geometry(num_sides, j + 0.51);
                Vector2d d = eigen_util::rotate90(p1 - p0).normalized() * 0.2;
                Vector2d q0 = p0 - d;
                Vector2d q1 = p0 + d;
                glColor3d((i == g.selected_variable ? 2 : 1) * Vector3d(0.6, 0.2, 0.0));
                draw_doublesided_arrow(q0, q1);
            }
        }
        const Vector3d color_table[4] = {
            Vector3d(0.0, 0.5, 0.0),
            Vector3d(0.5, 0.0, 0.5),
            Vector3d(0.0, 0.5, 0.5),
            Vector3d(0.5, 0.0, 0.0),
        };
        int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
        for (int i = num_sides; i < num_variables; ++i) {
            auto& variable_indicators = patchgen::get_variable_indicators(num_sides, param.pattern_id)[i - num_sides];
            for (auto& tag_pair : variable_indicators) {
                auto v0 = patchgen::find_tagged_vertex(patch, tag_pair.first);
                auto v1 = patchgen::find_tagged_vertex(patch, tag_pair.second);
                assert(v0.is_valid() && v1.is_valid());
                Vector2d p0 = patch.data(v0).laplaceDirect.value.head(2);
                Vector2d p1 = patch.data(v1).laplaceDirect.value.head(2);
                glColor3d((i == g.selected_variable ? 2 : 1) * color_table[i - num_sides]);
                draw_doublesided_arrow(p0, p1);
            }
        }
    }
    
    glColor3d(0, 0, 0);
    
    // draw number of edge subdivisions
    glPointSize(3);
    glLineWidth(3);
    for (int i = 0; i < num_sides; ++i) {
        Vector2d p0 = patchgen::get_boundary_geometry(num_sides, i + 0.49);
        Vector2d p1 = patchgen::get_boundary_geometry(num_sides, i + 0.51);
        Vector2d d = -eigen_util::rotate90(p1 - p0).normalized() * 0.3;
        Vector2d q = p0 + d;
        glPushMatrix();
        glTranslated(q.x() - 0.1, q.y() - 0.1, 0);
        glScaled(0.001, 0.001, 1);
        string str = boost::lexical_cast<string>(g.l[i]);
        for (char c : str)
            glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, c);
        glPopMatrix();
    }
    
    // draw pattern type in 2D
    glMatrixMode(GL_PROJECTION);    glPushMatrix();    glLoadIdentity();   gluOrtho2D(0, g.camera.width, 0, g.camera.height);
    glMatrixMode(GL_MODELVIEW );    glPushMatrix();    glLoadIdentity();
    string fname = demo::get_fname(param);
    glRasterPos2i(100, 10);
    for (char c : fname)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c);
    glMatrixMode(GL_PROJECTION);    glPopMatrix();
    glMatrixMode(GL_MODELVIEW );    glPopMatrix();
}
void reshape(int width, int height) {
    glut_util::defaultcb::reshape(width, height);
    g.camera.reshape(width, height);
}
Vector2d get_mouse_pos(int x, int y) {
        Vector3d world_xyz = unproject(Vector3d(x, y, 0));
        Vector3d& camera_center = g.camera.center;
        Vector3d& camera_eye    = g.camera.eye;
        return camera_center.head(2) + (world_xyz - camera_center).head(2) * camera_eye.z() / (camera_eye.z() - world_xyz.z());
}
void mouse(int glut_button, int state, int x, int y) {
    if (TwEventMouseButtonGLUT(glut_button, state, x, y)) return glutPostRedisplay();
    y = g.camera.height - y;
    
    bool shift_pressed = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0;
    bool ctrl_pressed  = (glutGetModifiers() & GLUT_ACTIVE_CTRL ) != 0;
    bool alt_pressed   = (glutGetModifiers() & GLUT_ACTIVE_ALT  ) != 0;
    
    if (state == GLUT_UP) {
        if (g.camera.drag_mode != Camera::DragMode::NONE)
            return g.camera.mouse_up();
        else if (g.selected_vertex.is_valid())
            return g.selected_vertex.invalidate();
    }
    
    if (alt_pressed) {
        g.camera.mouse_down(x, y, ctrl_pressed ? Camera::DragMode::ZOOM : Camera::DragMode::PAN);
        return;
    }
    
    Vector2d mouse_pos = get_mouse_pos(x, y);
    int num_sides = g.param.get_num_sides();
    
    if (g.mode == Globals::Mode::ChangeSubdiv) {
        MinSelector<int> selected_side(-1);
        for (int i = 0; i < num_sides; ++i) {
            for (int j = 0; j < g.l[i]; ++j) {
                double t0 = i +  j      / static_cast<double>(g.l[i]);
                double t1 = i + (j + 1) / static_cast<double>(g.l[i]);
                Vector2d p0 = patchgen::get_boundary_geometry(num_sides, t0);
                Vector2d p1 = patchgen::get_boundary_geometry(num_sides, t1);
                auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
                if (!dist) continue;
                selected_side.update(*dist, i);
            }
        }
        g.selected_side = selected_side.value;
    
    } else if (g.mode == Globals::Mode::AdjustParam) {
        patchgen::PatchParam& param = g.param;
        demo    ::Patch     & patch = g.patch;
        
        if (num_sides == 3 && param.pattern_id == 0) return;
        if (num_sides == 4 && param.pattern_id == 0) return;
        if (num_sides == 5 && param.pattern_id == 0) return;
        
        MinSelector<int> selected_variable(-1);
        int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
        for (int i = 0; i < num_variables; ++i) {
            if (i < num_sides) {
                int j = param.permutation[i];
                Vector2d p = patchgen::get_boundary_geometry(num_sides, j + 0.49);
                Vector2d d = eigen_util::rotate90(patchgen::get_boundary_geometry(num_sides, j + 0.51) - p).normalized() * 0.1;
                Vector2d p0 = p - d;
                Vector2d p1 = p + d;
                auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
                if (!dist) continue;
                selected_variable.update(*dist, i);
            } else {
                auto& variable_indicators = patchgen::get_variable_indicators(num_sides, param.pattern_id)[i - num_sides];
                for (auto& tag_pair : variable_indicators) {
                    auto v0 = patchgen::find_tagged_vertex(patch, tag_pair.first);
                    auto v1 = patchgen::find_tagged_vertex(patch, tag_pair.second);
                    assert(v0.is_valid() && v1.is_valid());
                    Vector2d p0 = patch.data(v0).laplaceDirect.value.head(2);
                    Vector2d p1 = patch.data(v1).laplaceDirect.value.head(2);
                    auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
                    if (!dist) continue;
                    selected_variable.update(*dist, i);
                }
            }
        }
        g.selected_variable = selected_variable.value;
    
    } else if (g.mode == Globals::Mode::MoveVertex) {
        MinSelector<demo::Patch::VHandle> selected_vertex;
        for (auto v : g.patch.vertices()) {
            if (g.patch.is_boundary(v)) continue;       // don't allow moving boundary vertices
            double dist = (mouse_pos - g.patch.data(v).laplaceDirect.value.head(2)).norm();
            selected_vertex.update(dist, v);
        }
        g.selected_vertex = selected_vertex.value;
    }
    
    glutPostRedisplay();
}
void motion(int x, int y) {
    if (TwEventMouseMotionGLUT(x, y)) return glutPostRedisplay();
    y = g.camera.height - y;
    
    if (g.camera.drag_mode != Camera::DragMode::NONE) {
        g.camera.mouse_move(x, y);
        return glutPostRedisplay();
    }
    if (g.selected_vertex.is_valid()) {
        Vector2d mouse_pos = get_mouse_pos(x, y);
        g.patch.data(g.selected_vertex).laplaceDirect.value << mouse_pos, 0;
        return glutPostRedisplay();
    }
}
// <<GLUT callback functions---------------------------------------

int main_glut(int argc, char* argv[]) {
    glut_util::init(argc, argv, GLUT_DOUBLE | GLUT_RGBA, 0.8, true, "patchgen_demo",
        display_pre,
        display_main,
        glut_util::defaultcb::display_post,
        reshape,
        glut_util::defaultcb::keyboard,
        nullptr,
        glut_util::defaultcb::special,
        mouse,
        motion);
    
    // GLEW
    glewInit();
    
    // AntTweakBar
    TwInit(TW_OPENGL, NULL);
    TwGLUTModifiersFunc([](){return glutGetModifiers();});
    
    init_gl();
    init_bar();
    g.init();
    
    glutMainLoop();
    return 0;
}
