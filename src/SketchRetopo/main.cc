#include <cstdlib>
#include <kt84/graphics/graphics_util.hh>
#include <AntTweakBar.h>
#include <kt84/glut_util.hh>
#include "SketchRetopo.hh"
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace kt84;
using namespace kt84::graphics_util;

// google-breakpad ==================================================================
//#define USE_GOOGLE_BREAKPAD

#ifdef USE_GOOGLE_BREAKPAD
#   ifdef WIN32

#include "client/windows/handler/exception_handler.h"
#include <csignal>

bool google_breakpad_callback(
    const wchar_t *dump_path,
    const wchar_t *id,
    void *context,
    EXCEPTION_POINTERS *exinfo,
    MDRawAssertionInfo *assertion,
    bool succeeded)
{
    if (succeeded) {
        printf("dump guid is %ws\n", id);
    } else {
        printf("dump failed\n");
    }
    fflush(stdout);
    
    return succeeded;
}

void __cdecl sigabrt_handler(int) { *static_cast<char*>(0) = 1; }
void __cdecl sigint_handler (int) { *static_cast<char*>(0) = 1; }

#   endif
#endif
// ================================================================== google-breakpad

namespace {
    auto& core = SketchRetopo::get_instance();
}

void display_pre() {
    stringstream window_title;
    window_title << "SketchRetopo";
    if (!core.configTemp.autoSave.filename.empty())
        window_title << " - " << core.configTemp.autoSave.filename;
    if (core.configTemp.autoSave.unsaved)
        window_title << "*";
    glutSetWindowTitle(window_title.str().c_str());
    
    glut_util::defaultcb::display_pre();
    
    // background color ramp ====================================================================
    glMatrixMode(GL_PROJECTION);    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW );    glLoadIdentity();
    glDepthMask(GL_FALSE);
    glBegin(GL_QUADS);
    glColor3f(core.configRender.bgcolor_bottom);    glVertex2d(-1, -1);    glVertex2d( 1, -1);    // bottom
    glColor3f(core.configRender.bgcolor_top   );    glVertex2d( 1,  1);    glVertex2d(-1,  1);    // top
    glEnd();
    glDepthMask(GL_TRUE);
    // ==================================================================== background color ramp
    
    // set projection matrix
    double zNear = core.camera->center_to_eye().norm() * 0.1;
    double zFar  = zNear * 10 + core.basemesh.boundingBox_diagonal_norm() * 10;
    double aspect_ratio = core.camera->width / static_cast<double>(core.camera->height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (core.configRender.use_ortho) {
        double ortho_y = zNear * 2.5;
        double ortho_x  = ortho_y * aspect_ratio;
        glOrtho(-ortho_x, ortho_x, -ortho_y, ortho_y, zNear, zFar);
    } else {
        gluPerspective(40, aspect_ratio, zNear, zFar);
    }
    
    // set modelview matrix
    auto eye    = core.camera->get_eye();
    auto center = core.camera->center;
    auto up     = core.camera->get_up();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x(), eye.y(), eye.z(), center.x(), center.y(), center.z(), up.x(), up.y(), up.z());
    
    core.display_pre();
}
void display_main() {
    core.state->display();
}
void display_post() {
    core.display_post();
    glut_util::defaultcb::display_post();
}
void reshape(int width, int height) {
    core.camera_free   .reshape(width, height);
    core.camera_upright.reshape(width, height);
    glut_util::defaultcb::reshape(width, height);
}
void keyboard(unsigned char key, int x, int y) {
    if (TwEventKeyboardGLUT(key, x, y)) return glutPostRedisplay();
    
    if (!core.common_keyboard(key, x, y))
        core.state-> keyboard(key, x, y);
    
    glutPostRedisplay();
}
void keyboardup(unsigned char key, int x, int y) {
    if (!core.common_keyboardup(key, x, y))
        core.state-> keyboardup(key, x, y);
    
    glutPostRedisplay();
}
void mouse(int glut_button, int state, int x, int y) {
    if (TwEventMouseButtonGLUT(glut_button, state, x, y)) {
        glutPostRedisplay();
        return;
    }
    y = core.camera->height - y;
    
    Button button =
        glut_button == GLUT_LEFT_BUTTON   ? Button::LEFT   :
        glut_button == GLUT_MIDDLE_BUTTON ? Button::MIDDLE : Button::RIGHT;
	
    bool shift_pressed = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0;
    bool ctrl_pressed  = (glutGetModifiers() & GLUT_ACTIVE_CTRL ) != 0;
    bool alt_pressed   = (glutGetModifiers() & GLUT_ACTIVE_ALT  ) != 0;
    
    if (state == GLUT_DOWN) {
        if (!core.common_mouse_down(x, y, button, shift_pressed, ctrl_pressed, alt_pressed))
            core.state-> mouse_down(x, y, button, shift_pressed, ctrl_pressed, alt_pressed);
    } else {
        if (!core.common_mouse_up(x, y, button, shift_pressed, ctrl_pressed, alt_pressed))
            core.state-> mouse_up(x, y, button, shift_pressed, ctrl_pressed, alt_pressed);
    }
    
    glutPostRedisplay();
}
void motion(int x, int y) {
    if (TwEventMouseMotionGLUT(x, y)) {
        glutPostRedisplay();
        return;
    }
    y = core.camera->height - y;
    
    if (!core.common_mouse_move(x, y))
        core.state-> mouse_move(x, y);
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    // google-breakpad
#ifdef USE_GOOGLE_BREAKPAD
#   ifdef WIN32
    new google_breakpad::ExceptionHandler(
        L".", NULL, google_breakpad_callback, NULL, google_breakpad::ExceptionHandler::HANDLER_ALL);
    signal(SIGABRT, sigabrt_handler);
    signal(SIGINT , sigint_handler );
#   endif
#endif
    
    //------+
    // GLUT |
    //------+
    glut_util::init(argc, argv, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL, 0.8, true, "SketchRetopo",
        display_pre,
        display_main,
        display_post,
        reshape,
        keyboard,
        keyboardup,
        nullptr,
        mouse,
        motion,
        motion);
    
    // GLEW
    glewInit();
    
    // AntTweakBar
    TwInit(TW_OPENGL, NULL);
    TwGLUTModifiersFunc([]() { return glutGetModifiers(); });
    
    // init core data
    core.init();
    
    // process command line argument
    if (argc == 2) {
        string fname = argv[1];
        auto fname_ext = fname.substr(fname.size() - 4, 4);
        boost::algorithm::to_lower(fname_ext);
        if (fname_ext == ".obj" || fname_ext == ".off" || fname_ext == ".ply" || fname_ext == ".stl")
            core.import_basemesh(fname);
        else if (fname_ext == ".xml")
            core.xml_load(fname);
    }
    
    glutMainLoop();
}
