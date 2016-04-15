#include "Patch.hh"
#include <kt84/util.hh>
#include <kt84/graphics/graphics_util.hh>
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

void demo::Patch::draw() const {
    glBegin(GL_QUADS);
    for (auto f : faces())
        for (auto v = cfv_iter(f); v.is_valid(); ++v)
            glVertex3d(data(*v).laplaceDirect.value);
    glEnd();
}
void demo::Patch::draw_singularities() const {
    for (auto v : vertices()) {
        int valence = 0;
        for (auto vv = cvv_iter(v); vv.is_valid(); ++vv, ++valence);
        if (is_boundary(v)) ++valence;
        if (data(v).patchgen.corner_index != -1) ++valence;
        if (valence == 4) continue;
        assert(valence == 3 || valence == 5);
        Vector3d p = data(v).laplaceDirect.value;
        Vector3d color = valence == 3 ? Vector3d(0, 0.6, 0.9) : Vector3d(1, 0.7, 0);
        glColor3d(0, 0, 0);    glPointSize(12);   glBegin(GL_POINTS);    glVertex3d(p);    glEnd();
        glColor3d(color);      glPointSize(10);   glBegin(GL_POINTS);    glVertex3d(p);    glEnd();
    }
}
