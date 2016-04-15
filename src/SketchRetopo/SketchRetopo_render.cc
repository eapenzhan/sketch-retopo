#include "SketchRetopo.hh"
#include <kt84/glut_util.hh>
#include <kt84/graphics/graphics_util.hh>
#include "helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

void SketchRetopo::display_pre() {
    // auto save |
    //-----------+
    if ( configTemp.autoSave.enabled          &&
        !configTemp.autoSave.filename.empty() &&
         configTemp.autoSave.unsaved          &&
        (clock() - configTemp.autoSave.last_time) / CLOCKS_PER_SEC > configTemp.autoSave.interval)
    {
        configTemp.autoSave.last_time = clock();
        
        // update ausoSave.filename if it has a specific form of <fname_core>.<fname_num>.xml
        string fname = configTemp.autoSave.filename;
        string fname_core;
        int fname_num;
        if (decompose_xml_fname(fname, fname_core, fname_num)) {
            stringstream ss;
            ss << fname_core << '.' << (fname_num + 1) << ".xml";
            fname = ss.str();
        }
        
        xml_save(fname);
    }
    
    // update brush size |
    //-------------------+
    double dist_center_to_eye = camera->center_to_eye().norm();
    configTemp.snapSize            .set_camera_distance(dist_center_to_eye);
    configTemp.brushSize_spine     .set_camera_distance(dist_center_to_eye);
    configTemp.brushSize_autocmpl  .set_camera_distance(dist_center_to_eye);
    configTemp.brushSize_moveVertex.set_camera_distance(dist_center_to_eye);
    configTemp.quadSize   .adjust_rate = dist_center_to_eye * 0.0002;
    configTemp.segmentSize.adjust_rate = dist_center_to_eye * 0.0002;
    
    // configure basic shader to default value |
    //-----------------------------------------+
    program_basic.enable();
    program_basic.set_uniform_3("uni_light_pos", configRender.light_pos);
    program_basic.set_uniform_1("uni_use_shading", 0);
    program_basic.set_uniform_1("uni_symmetric_mirror", 0);
    program_basic.disable();
        
    stencil_index = 0;
    
    if (configRender.mode == ConfigRender::Mode::DEFAULT      ) render_default      ();
    if (configRender.mode == ConfigRender::Mode::QUADMESH_ONLY) render_quadmesh_only();
    if (configRender.mode == ConfigRender::Mode::BASEMESH_ONLY) render_basemesh_only();
}
void SketchRetopo::display_post() {
    static const double brush_alpha = 0.1;
    static const int slices = 36;
    static const int stacks = 18;
    
    if (common_pn_mouse) {
        // configure basic shader to default value |
        //-----------------------------------------+
        program_basic.enable();
        program_basic.set_uniform_1("uni_use_shading", 1);
        program_basic.set_uniform_1("uni_symmetric_mirror", 0);
        
        // draw current mouse position
        glDepthMask(GL_FALSE);
        glPushMatrix();
        glTranslated(common_pn_mouse->x(), common_pn_mouse->y(), common_pn_mouse->z());
        
        Vector4d brush_color;
        brush_color << state->state_color, brush_alpha;
        
        if (configTemp.quadSize.is_adjusting()) {
            // quad size
            glColor4d(0, 0.2, 0.5, brush_alpha);
            glutSolidSphere(configTemp.quadSize() * 0.7, slices, stacks);
            
        } else if (configTemp.segmentSize.is_adjusting()) {
            // segment sieze
            glColor4d(0.5, 0.2, 0, brush_alpha);
            glutSolidSphere(configTemp.segmentSize() * 3, slices, stacks);
            
        } else {
            // snap brush size
            if (state->show_brushSize_snap()) {
                glColor4d((configTemp.snapSize.is_adjusting() ? 1.0 : 0.8) * Vector4d(1, 1, 1, brush_alpha));
                double snapSize = configTemp.snapSize();
                //if (!configTemp.snapSize.is_adjusting()) {        // just to see how big the snapSize gets depending on the mouse speed
                //    double delta_max = 20;
                //    double mouse_delta = min<double>((common_mouse_pos - common_mouse_pos_prev).cast<double>().norm(), delta_max);
                //    snapSize *= 1 + mouse_delta / delta_max;
                //}
                glutSolidSphere(snapSize, slices, stacks);
            }
            
            // spine brush size
            if (state == &stateSpine) {
                glColor4d((configTemp.brushSize_spine.is_adjusting() ? 1.0 : 0.8) * brush_color);
                glutSolidSphere(configTemp.brushSize_spine(), slices, stacks);
            }
            
            // autocmpl brush size
            if (state == &stateAutocmpl && configTemp.brushSize_autocmpl.is_adjusting()) {
                glColor4d(brush_color);
                glutSolidSphere(configTemp.brushSize_autocmpl(), slices, stacks);
            }
            
            // moveVertex brush size
            if (state == &stateMoveVertex) {
                glColor4d((configTemp.brushSize_moveVertex.is_adjusting() ? 1.0 : 0.8) * brush_color);
                glutSolidSphere(configTemp.brushSize_moveVertex(), slices, stacks);
            }
        }
        
        program_basic.disable();
        
        // single point under mouse cursor
        glPointSize(3);
        glDisable(GL_DEPTH_TEST);
        glBegin(GL_POINTS);
        glColor3d(1, 0, 0);
        glVertex3d(0, 0, 0);
        glEnd();
        
        glPopMatrix();
    }
    
    // >>>> 2D draw mode
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);    glPushMatrix();     glLoadIdentity();       gluOrtho2D(0, camera->width, 0, camera->height);
    glMatrixMode(GL_MODELVIEW );    glPushMatrix();     glLoadIdentity();
    
    // hiding stroke |
    //---------------+
    if (!hide_stroke.empty()) {
        glLineWidth(5);
        glColor3d(1, 0, 0);
        glBegin(GL_LINE_STRIP);
        for (auto& p : hide_stroke) glVertex2d(p);
        glEnd();
    }
    
    // display current state name |
    //----------------------------+
    glPointSize(3);
    glLineWidth(3);
    glTranslated(30, 30, 0);
    glScaled(0.25, 0.25, 0.25);
    glColor3dv(&state->state_color[0]);
    for_each(state->state_name.begin(), state->state_name.end(),
        [] (char c) { glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, c); });
    
    glMatrixMode(GL_PROJECTION);    glPopMatrix();
    glMatrixMode(GL_MODELVIEW );    glPopMatrix();
    // <<<< 2D draw mode
}

/*
    | rendering styles |
    +------------------+
    data to be rendered:
        basemesh face
        basemesh edge
        patches mesh face
        patches mesh edge
        patches irregular vertices
        curve network segments
        curve network corners
    
    (1)default:
        basemesh face
        basemesh edge
        (no depth mask, with stencil test)
        patches mesh face (alpha)
        patches mesh edge          (quad mesh level only)
        patches irregular vertices (quad mesh level only)
        curve network segments
        curve network corners
    
    (2)quad mesh only:
        patches mesh face (opaque)
        patches mesh edge (polygon offset, different color for boundary edges)
        patches irregular vertices (polygon offset)
        curve network corners (slight offset?)
        base mesh face (depth only)
    
    (3)basemesh only:
        basemesh face
        basemesh edge
*/
    
void SketchRetopo::render_default() {
    // basemesh |
    //----------+
    program_basemesh.enable();
    program_basemesh.set_uniform_1<int>("show_expmap"  , configRender.show_expmap  );
    program_basemesh.set_uniform_1<float>("overlay_alpha" , configRender.overlay_alpha );
    program_basemesh.set_uniform_1<int>("overlay_v_flipped" , configRender.overlay_v_flipped ? 1 : 0);
    program_basemesh.set_uniform_1<int>("symmetric", configSaved .symmetric);
    if (configRender.show_snake_energy) {
        program_basemesh.disable();
    }
    glActiveTexture(GL_TEXTURE1);
    texture_overlay.bind();
    glActiveTexture(GL_TEXTURE0);
    material_texture.bind();
    basemesh.render_fill();
    if (configRender.basemesh_render_line) {
        program_basemesh.enable();
        program_basemesh.set_uniform_1<int>("line_mode", 1);
        glLineWidth(1);
        basemesh.render_line();
        program_basemesh.set_uniform_1<int>("line_mode", 0);
    }
    program_basemesh.disable();
    
    // enable stencil buffer (for z-fight-free rendering) |
    //----------------------------------------------------+
    glEnable(GL_STENCIL_TEST);
    glDepthMask(GL_FALSE);

    // patches |
    //---------+
    // fill mode
    glStencilFunc(GL_GREATER, ++stencil_index, ~0);
    for (auto& p : curvenetwork.patches) {
        if (p.is_hidden)
            continue;
        int num_corners = p.num_corners();
        Vector3d color_rgb = configRender.patch_color_rgb[p.num_corners()];
        if (p.is_failure) color_rgb << 0, 0, 0;
        
        glColor4d(color_rgb, configRender.patch_color_alpha);
        p.render_phong_fill();
    }
    
    if (configSaved.symmetric) {
        // symmetric mirror
        program_basic.enable();
        glColor4d(0.5, 0.5, 0.5, 0.5);
        program_basic.set_uniform_1("uni_symmetric_mirror", 1);
        for (auto& p : curvenetwork.patches) {
            if (p.is_hidden)
                continue;
            p.render_phong_fill();
        }
        program_basic.disable();
    }
    
    if (configRender.always_show_quads || state->get_editLevel() == State::EditLevel::QUAD_MESH) {
        // line mode
        glStencilFunc(GL_GREATER, ++stencil_index, ~0);
        glColor3d(0, 0, 0);
        glLineWidth(1);
        for (auto& p : curvenetwork.patches) {
            if (p.is_failure || p.is_hidden)
                continue;
            p.render_phong_line();
        }
        
        if (configSaved.symmetric) {
            // symmetric mirror
            program_basic.enable();
            glColor4d(0.f, 0.f, 0.f, 0.5f);
            for (auto& p : curvenetwork.patches) {
            if (p.is_failure || p.is_hidden)
                    continue;
                p.render_phong_line();
                p.render_phong_line_boundary();
            }
            program_basic.disable();
        }
    }
    
    // curve network elements (edgechains, corners/t-junctions, endpoints) |
    //---------------------------------------------------------------------+
    // edgechains
    glStencilFunc(GL_GREATER, ++stencil_index, ~0);
    curvenetwork.render_edgechains();
    
    // corners/t-junctions
    glStencilFunc(GL_GREATER, ++stencil_index, ~0);
    program_corner.enable();
    program_corner.set_uniform_1<int>("uni_symmetric_mirror", 0);
    program_corner.set_uniform_1<float>("fan_radius", camera->center_to_eye().norm() * configRender.corner_radius);
    curvenetwork.render_corners();
    program_corner.disable();
    
    // endpoints
    glStencilFunc(GL_GREATER, ++stencil_index, ~0);
    curvenetwork.render_endpoints();
    
    // irregular vertices
    if (configRender.always_show_quads || state->get_editLevel() == State::EditLevel::QUAD_MESH) {
        glStencilFunc(GL_GREATER, ++stencil_index, ~0);
        for (auto& p : curvenetwork.patches) {
            if (p.is_failure || p.is_hidden)
                continue;
            p.render_irregular_vertices();
        }
    }
    
    glDepthMask(GL_TRUE);
    glDisable(GL_STENCIL_TEST);
    
    // coordinate axis |
    //-----------------+
    if (configRender.show_axis) {
        glLineWidth(2);
        glBegin(GL_LINES);
        glColor3d(1, 0, 0);     glVertex3d(0, 0, 0);    glVertex3d(1, 0, 0);
        glColor3d(0, 1, 0);     glVertex3d(0, 0, 0);    glVertex3d(0, 1, 0);
        glColor3d(0, 0, 1);     glVertex3d(0, 0, 0);    glVertex3d(0, 0, 1);
        glEnd();
    }
    
    // bounding box |
    //--------------+
    if (configRender.show_bbox) {
        auto& bbox = basemesh.boundingBox;
        glPushMatrix();
        glTranslated(bbox.center().x(), bbox.center().y(), bbox.center().z());
        glScaled(bbox.sizes().x(), bbox.sizes().y(), bbox.sizes().z());
        
        glLineWidth(2);
        glColor3d(1, 1, 1);
        glutWireCube(1);
        
        glPopMatrix();
    }
}
void SketchRetopo::render_quadmesh_only() {
    // patches |
    //---------+
    program_basic.enable();
    
    /*
    1. render patches in fill mdoe
    // right side
    set mirror off
    for each patch {
        if (show_auxiliary)
            set color based on num_corners
        else
            set color white
        render_fill
    }
    // left side
    set mirror on
    for each patch {
        if (!show_left_side)
            set color gray
        else if (show_auxiliary)
            set color based on num_corners
        else
            set color white
        render_fill
    }
    
    2. render patches in line mode
    set line with to 1
    set color black
    for each patch {
        // right side
        set mirror off;
        render_line;
        if (!show_auxiliary) render_line_border;
        
        // left side
        set mirror on;
        render_line;
        if (!show_auxiliary || !show_left_side) render_line_border;
    }
    
    3. render auxiliary data
    if (show_auxiliary) {
        // patch boundary
        set line width to 2
        set point size to 10
        for each patch {
            set mirror off;
            set color yellow;
            render_line_border;
            render_irregular_vertices;
            
            if (show_left_side) {
                set mirror on;
                set color yellow;
                render_line_border;
                render_irregular_vertices
            }
        }
        // corner
        set mirror off;
        curvenetwork.render_corner
        if (symmetric && show_left_side) {
            set mirror on;
            curvenetwork.render_corner

        }
    }
            render_corner
    */
    
    // render in fill mode
    program_basic.set_uniform_1("uni_use_shading", 1);
    // right side
    for (auto& p : curvenetwork.patches) {
        if (p.is_hidden)
            continue;
        int num_corners = p.num_corners();
        Vector4d color;
        if (configRender.quadmesh_only_show_auxiliary)
            color << configRender.patch_color_rgb[num_corners], 1;
        else
            color << configRender.patch_color_rgb[4], 1;
        color.head(3) *= 1.25;
        
        if (p.is_failure) color.head(3).setZero();
        
        glColor4d(color);
        p.render_phong_fill();
    }
    
    if (configSaved.symmetric) {
        // left side
        program_basic.set_uniform_1("uni_symmetric_mirror", 1);
        for (auto& p : curvenetwork.patches) {
            if (p.is_hidden)
                continue;
            int num_corners = p.num_corners();
            Vector4d color;
            if (!configRender.quadmesh_only_show_left_side)
                color << 0.5, 0.5, 0.5, 1;
            else if (configRender.quadmesh_only_show_auxiliary)
                color << configRender.patch_color_rgb[num_corners], 1;
            else
                color << configRender.patch_color_rgb[4], 1;
            color.head(3) *= 1.25;
            
            if (p.is_failure)  color.head(3).setZero();
            
            glColor4d(color);
            p.render_phong_fill();
        }
    }
    
    // render in line mode
    glLineWidth(configRender.edge_width_interior);
    program_basic.set_uniform_1("uni_use_shading", 0);
    glColor4d(0.f, 0.f, 0.f, 1.f);
    for (auto& p : curvenetwork.patches) {
        if (p.is_failure || p.is_hidden)
            continue;
        
        // right side
        program_basic.set_uniform_1("uni_symmetric_mirror", 0);
        p.render_phong_line();
        if (!configRender.quadmesh_only_show_auxiliary)
            p.render_phong_line_boundary();
        
        if (configSaved.symmetric) {
            // left side
            program_basic.set_uniform_1("uni_symmetric_mirror", 1);
            p.render_phong_line();
            if (!configRender.quadmesh_only_show_auxiliary || !configRender.quadmesh_only_show_left_side)
                p.render_phong_line_boundary();
        }
    }
    
    // render auxiliary data
    if (configRender.quadmesh_only_show_auxiliary) {
        // patch boundary edges
        glLineWidth(configRender.edge_width_boundary);
        glPointSize(configRender.singularity_radius);
        for (auto& p : curvenetwork.patches) {
            if (p.is_hidden)
                continue;
            // right side
            program_basic.set_uniform_1("uni_symmetric_mirror", 0);
            glColor3d(1, 0.8, 0);
            p.render_phong_line_boundary();
            p.render_irregular_vertices();
            
            if (configSaved.symmetric && configRender.quadmesh_only_show_left_side) {
                // left side
                program_basic.set_uniform_1("uni_symmetric_mirror", 1);
                glColor3d(1, 0.8, 0);
                p.render_phong_line_boundary();
                p.render_irregular_vertices();
            }
        }
        program_basic.disable();
        
        // curve network corner/t-junction vertices
        program_corner.enable();
        program_corner.set_uniform_1<float>("fan_radius", camera->center_to_eye().norm() * configRender.corner_radius);
        program_corner.set_uniform_1("uni_symmetric_mirror", 0);
        curvenetwork.render_corners();
        if (configSaved.symmetric && configRender.quadmesh_only_show_left_side) {
            program_corner.set_uniform_1("uni_symmetric_mirror", 1);
            curvenetwork.render_corners();
        }
    }
    
    program_basic.disable();
    
    // coordinate axis |
    //-----------------+
    if (configRender.show_axis) {
        glLineWidth(2);
        glBegin(GL_LINES);
        glColor3d(1, 0, 0);     glVertex3d(0, 0, 0);    glVertex3d(1, 0, 0);
        glColor3d(0, 1, 0);     glVertex3d(0, 0, 0);    glVertex3d(0, 1, 0);
        glColor3d(0, 0, 1);     glVertex3d(0, 0, 0);    glVertex3d(0, 0, 1);
        glEnd();
    }
    
    // bounding box |
    //--------------+
    if (configRender.show_bbox) {
        auto& bbox = basemesh.boundingBox;
        glPushMatrix();
        glTranslated(bbox.center().x(), bbox.center().y(), bbox.center().z());
        glScaled(bbox.sizes().x(), bbox.sizes().y(), bbox.sizes().z());
        
        glLineWidth(2);
        glColor3d(1, 1, 1);
        glutWireCube(1);
        
        glPopMatrix();
    }
    
    // render depth at least to get mouse position projected onto basemesh
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    basemesh.render_fill();
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}
void SketchRetopo::render_basemesh_only() {
    // basemesh |
    //----------+
    program_basemesh.enable();
    program_basemesh.set_uniform_1<int>("show_expmap"  , configRender.show_expmap );
    program_basemesh.set_uniform_1<float>("overlay_alpha" , configRender.overlay_alpha);
    program_basemesh.set_uniform_1<int>("overlay_v_flipped" , configRender.overlay_v_flipped ? 1 : 0);
    program_basemesh.set_uniform_1<int>("symmetric", configSaved .symmetric && !configRender.quadmesh_only_show_left_side);
    if (configRender.show_snake_energy) {
        program_basemesh.disable();
    }
    glActiveTexture(GL_TEXTURE1);
    texture_overlay.bind();
    glActiveTexture(GL_TEXTURE0);
    material_texture.bind();
    basemesh.render_fill();
    if (configRender.basemesh_render_line) {
        program_basemesh.enable();
        program_basemesh.set_uniform_1<int>("line_mode", 1);
        glLineWidth(1);
        basemesh.render_line();
        program_basemesh.set_uniform_1<int>("line_mode", 0);
    }
    program_basemesh.disable();
}
