#include "SketchRetopo.hh"
#include <algorithm>
#include <kt84/zenity_util.hh>
#include <kt84/tw_util.hh>
#include <kt84/graphics/phong_tessellation.hh>
#include "helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;

namespace {
    auto& core = SketchRetopo::get_instance();
}

BaseMesh::BaseMesh()
    : normalSmoothIter(0)
    , featureSmoothIter(0)
    , snake_energy_max()
    , snake_energy_min()
    , use_feature_as_energy(true)
    , negate_energy(true)
    , snake_feature_type(SnakeFeatureType::MEAN_CURVATURE)
    , has_texture(false)
{
    request_face_normals();
    request_vertex_normals();
    request_halfedge_texcoords2D();
}

ConfigRender::ConfigRender()
    : mode(Mode::DEFAULT)
    , always_show_quads(true)
    , basemesh_render_line(false)
    , show_axis(false)
    , show_bbox(false)
    , use_ortho  (false)
    , show_expmap (false)
    , overlay_alpha(0.0)
    , overlay_v_flipped(true)
    , show_snake_energy(false)
    , snake_energy_scale(10.0)
    , auto_camera_center(true)
    , bgcolor_bottom(0.25, 0.25, 0.25)
    , bgcolor_top   (0   , 0   ,    0)
    , corner_radius(0.005)
    , edge_width_interior(1)
    , edge_width_boundary(2)
    , singularity_radius(10)
    , patch_color_alpha(0.5)
    , light_pos(50, 100, 200)
    , quadmesh_only_show_auxiliary(false)
    , quadmesh_only_show_left_side(false)
    , turntable_mode(false)
    , turntable_speed(0.02)
{
    patch_color_rgb[0] << 0, 0, 0;
    patch_color_rgb[1] << 0, 0, 0;
    patch_color_rgb[2] << 0.2, 0.1, 0.9;
    patch_color_rgb[3] << 0.5, 0.5, 0.7;
    patch_color_rgb[4] << 0.7, 0.7, 0.7;
    patch_color_rgb[5] << 0.7, 0.6, 0.5;
    patch_color_rgb[6] << 0.2, 0.8, 0.3;
}

ConfigSaved::ConfigSaved()
    : projOffset(0)
    , symmetric(false)
{}

const double ConfigTemp::default_quadSize_ratio    = 0.03;
const double ConfigTemp::default_segmentSize_ratio = 0.002;

ConfigTemp::ConfigTemp()
    : snapSize            (0.02, util::dbl_max(), 0, true , 0.001)
    , brushSize_spine     (0.06, util::dbl_max(), 0, true , 0.001)
    , brushSize_autocmpl  (0.08, util::dbl_max(), 0, true , 0.001)
    , brushSize_moveVertex(0.06, util::dbl_max(), 0, true , 0.001)
    , quadSize            (0   , util::dbl_max(), 0, false, 0    )       // these are set when loading basemesh
    , segmentSize         (0   , util::dbl_max(), 0, false, 0    )
    , autocmpl_preview_interval(200)
    , cylinder_num_div(4)
    , loop_threshold_normal(0.3)
    , loop_threshold_area  (0.3)
    , snakes_internal1(10000000)
    , snakes_internal2(10000000)
    , snakes_external (1)
    , snakes_damping  (0.5)
{
}

AutoSave::AutoSave()
    : enabled (true)
    , interval(30)
    , unsaved(false)
    , last_time()
{}

SketchRetopo::SketchRetopo()
    : camera(&camera_free)
    , common_dragMode(CommonDragMode::None)
    , state(0)
    , bar_main(0)
    , stencil_index()
    , memento(150)
{
    fill(common_key_pressed, common_key_pressed + 256, false);
    state = &stateSketch;
}

void SketchRetopo::init() {
    //-----------------+
    // OpenGL settings |
    //-----------------+
    glEnable(GL_DEPTH_TEST);
    // polygon offset
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    // color material (diffuse)
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);	
    // stencil op
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    // blend enabled
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // anti-aliasing. somehow not working?
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);    
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
    // init material texture
    texture_overlay.init();
    texture_overlay.bind();
    texture_overlay.set_default_param();
    material_texture.init();
    material_texture.bind();
    material_texture.set_default_param();
    string base64_stringified[] = {
        #include "../../resource/material_default.base64.txt.stringify.txt"
    };
    material_texture_png_base64 = de_stringify(base64_stringified, sizeof(base64_stringified));
    material_texture_update();
    
    bar_init();
    
    program_init();
    
    memento.init(150);
    
    state->init();
    camera_free.auto_flip_y = camera_upright.auto_flip_y = false;
}

void SketchRetopo::geodesic_init() {
    vector<double> points;
    points.reserve(basemesh.n_vertices() * 3);
    for (auto v : basemesh.vertices()) {
        auto p = basemesh.point(v);
        points.push_back(p[0]);
        points.push_back(p[1]);
        points.push_back(p[2]);
    }
    
    vector<unsigned> faces;
    faces.reserve(basemesh.n_faces() * 3);
    for (auto f : basemesh.faces()) {
        for (auto v : basemesh.fv_range(f))
            faces.push_back(v.idx());
    }
    
    geodesic_mesh.initialize_mesh_data(points, faces);
    geodesic_algorithm = geodesic::GeodesicAlgorithmExact(&geodesic_mesh);
}

void SketchRetopo::embree_init() {
    // triangles
    auto triangles = static_cast<embree::BuildTriangle*>(embree::rtcMalloc(basemesh.n_faces() * sizeof(embree::BuildTriangle)));
    for (auto f : basemesh.faces()) {
        if (basemesh.data(f).is_hidden) continue;
        
        auto v = basemesh.fv_iter(f);
        
        int f_idx = f.idx();
        int v0_idx = v->idx();
        int v1_idx = (++v)->idx();
        int v2_idx = (++v)->idx();
        
        triangles[f_idx] = embree::BuildTriangle(v0_idx, v1_idx, v2_idx, f_idx);
    }
    
    // vertices
    auto vertices  = static_cast<embree::BuildVertex*>(embree::rtcMalloc(basemesh.n_vertices() * sizeof(embree::BuildVertex)));
    for (auto v : basemesh.vertices()) {
        OpenMesh::Vec3f p(basemesh.point(v));
        vertices[v.idx()] = embree::BuildVertex(p[0], p[1], p[2]);
    }
    
    // build
    embree_accel = embree::rtcCreateAccel("default", "default", triangles, basemesh.n_faces(), vertices, basemesh.n_vertices());
    embree_intersector = embree_accel->queryInterface<embree::Intersector>();
    
    // warm up?
    cout << "waiting for embree to warm up\n";
    auto clk = clock();
    while (clock() - clk < CLOCKS_PER_SEC);
}

void SketchRetopo::program_init() {
    // basemesh
    {
        string src_frag_stringified[] = {
            #include "../../resource/basemesh.frag.stringify.txt"
        };
        string src_vert_stringified[] = {
            #include "../../resource/basemesh.vert.stringify.txt"
        };
        auto src_frag = de_stringify(src_frag_stringified, sizeof(src_frag_stringified));
        auto src_vert = de_stringify(src_vert_stringified, sizeof(src_vert_stringified));
        
        ShaderObject frag;                                      ShaderObject vert;
        frag.init(ShaderObject::Type::FRAGMENT_SHADER);         vert.init(ShaderObject::Type::VERTEX_SHADER);
        frag.set_source(src_frag);                              vert.set_source(src_vert);
        frag.compile();                                         vert.compile();
        
        program_basemesh.init();
        program_basemesh.attach(frag);
        program_basemesh.attach(vert);
        program_basemesh.link();
        program_basemesh.enable();
        program_basemesh.set_uniform_1<int>("material_texture", 0);
        program_basemesh.set_uniform_1<int>("texture_overlay", 1);
        program_basemesh.disable();
    }
    
    // basic
    {
        string src_frag_stringified[] = {
            #include "../../resource/basic.frag.stringify.txt"
        };
        string src_vert_stringified[] = {
            #include "../../resource/basic.vert.stringify.txt"
        };
        auto src_frag = de_stringify(src_frag_stringified, sizeof(src_frag_stringified));
        auto src_vert = de_stringify(src_vert_stringified, sizeof(src_vert_stringified));
        
        ShaderObject frag;                                      ShaderObject vert;
        frag.init(ShaderObject::Type::FRAGMENT_SHADER);         vert.init(ShaderObject::Type::VERTEX_SHADER);
        frag.set_source(src_frag);                              vert.set_source(src_vert);
        frag.compile();                                         vert.compile();
        
        program_basic.init();
        program_basic.attach(frag);
        program_basic.attach(vert);
        program_basic.link();
    }
    
    // corner
    {
        string src_frag_stringified[] = {
            #include "../../resource/corner.frag.stringify.txt"
        };
        string src_vert_stringified[] = {
            #include "../../resource/corner.vert.stringify.txt"
        };
        auto src_frag = de_stringify(src_frag_stringified, sizeof(src_frag_stringified));
        auto src_vert = de_stringify(src_vert_stringified, sizeof(src_vert_stringified));
        
        ShaderObject frag;                                      ShaderObject vert;
        frag.init(ShaderObject::Type::FRAGMENT_SHADER);         vert.init(ShaderObject::Type::VERTEX_SHADER);
        frag.set_source(src_frag);                              vert.set_source(src_vert);
        frag.compile();                                         vert.compile();
        
        program_corner.init();
        program_corner.attach(frag);
        program_corner.attach(vert);
        program_corner.link();
        program_corner.enable();
        program_corner.set_uniform_1("uni_symmetric_mirror", 0);
        program_corner.disable();
    }
}

void SketchRetopo::bar_init() {
    bar_main = TwNewBar("SketchRetopo");
    TwDefine("GLOBAL fontsize=3"); 
    TwDefine("SketchRetopo size='300 830' valueswidth=120 alpha=200 text=light"); 
    
    // file
    tw_util::AddButton(bar_main, "import_basemesh",
        [&](){
            string fname;
            if ((!configTemp.autoSave.unsaved || zenity_util::question("SketchRetopo", "Current curve data is not yet saved. Continue?")) &&   // confirm discarding unsaved data
                zenity_util::file_selection_load(fname, "SketchRetopo - import base mesh", "", "mesh files (*.obj, *.off, *.ply, *.stl)|*.obj *.off *.ply *.stl"))
            {
                if (!import_basemesh(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - import base mesh", "Error occurred!");
            }
        }, "group=file label='import base mesh'"  );
    tw_util::AddButton(bar_main, "export_basemesh",
        [&](){
            string fname;
            if (zenity_util::file_selection_save(fname, "SketchRetopo - export base mesh", "", "mesh files (*.obj, *.off, *.ply, *.stl)|*.obj *.off *.ply *.stl"))
                if (!export_basemesh(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - export base mesh", "Error occurred!");
        }, "group=file label='export base mesh'"  );
    tw_util::AddButton(bar_main, "export_retopomesh",
        [&](){
            string fname;
            if (zenity_util::file_selection_save(fname, "SketchRetopo - export retopo mesh", "", "mesh files (*.obj, *.off, *.ply, *.stl)|*.obj *.off *.ply *.stl"))
                if (!export_retopomesh(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - export retopo mesh", "Error occurred!");
        }, "group=file label='export retopo mesh'");
    
    TwAddSeparator(bar_main, NULL, " group='file'");
    tw_util::AddButton(bar_main, "xml_save",
        [&](){
            string fname;
            if (zenity_util::file_selection_save(fname, "SketchRetopo - save xml", "", "xml files|*.xml"))
                if (!xml_save(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - save xml", "Error occurred!");
        }, "group=file label='save xml'");
    tw_util::AddButton(bar_main, "xml_load",
        [&](){
            string fname;
            if ((!configTemp.autoSave.unsaved || zenity_util::question("SketchRetopo", "Current curve data is not yet saved. Continue?")) &&   // confirm discarding unsaved data
                zenity_util::file_selection_load(fname, "SketchRetopo - load xml", "", "xml files|*.xml"))
            {
                if (!xml_load(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - load xml", "Error occurred!");
            }
        }, "group=file label='load xml'");
    TwAddVarRW(bar_main, "autosave_enabled" , TW_TYPE_BOOLCPP, &configTemp.autoSave.enabled , "group=file label='autosave'");
    TwAddVarRW(bar_main, "autosave_interval", TW_TYPE_INT32  , &configTemp.autoSave.interval, "group=file label='autosave intrv' min=1");
    
    TwAddSeparator(bar_main, NULL, " group='file'");
    tw_util::AddButton(bar_main, "load_texture",
        [&](){
            string fname;
            if (zenity_util::file_selection_load(fname, "SketchRetopo - load texture", "", "png files|*.png"))
                if (!texture_overlay_load(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - load texture", "Error occurred!");
        }, "group=file");
    tw_util::AddButton(bar_main, "load_material",
        [&](){
            string fname;
            if (zenity_util::file_selection_load(fname, "SketchRetopo - load material", "", "png files|*.png"))
                if (!material_texture_load(fname))
                    zenity_util::msgbox(zenity_util::MsgBoxType::Error, "SketchRetopo - load material", "Error occurred!");
        }, "group=file");
    
    // mode
    tw_util::AddVarCB<EnumState>(bar_main, "mode", TwDefineEnum("ModeType", 0, 0),
        [&](const EnumState& value) {
            state_set(value);
            state->init();
        },
        [&](EnumState& value) {
            value = state_get();
        }, "enum='0{sketch}, 1{spine}, 2{autocmpl}, 3{laser}, 4{cylinder}, 5{dfm curve}, 6{edit corner}, 7{edit topo}, 8{edge loop}, 9{move vtx}'");
    
    // clear button
    tw_util::AddButton(bar_main, "clear", [&](){
        memento_store();
        curvenetwork.clear();
    }, nullptr);
    
    // config (saved)
    tw_util::AddVarCB<bool>(bar_main, "symmetric",
        [&](const bool& value){
            toggle_symmetric();
        }, [&](bool& value){
            value = configSaved.symmetric;
        }, "group=config label='symmetric'");
    tw_util::AddVarCB_default<double>(bar_main, "projOffset", configSaved.projOffset,
        [&](){
            reproject_all();
        }, "group=config label='proj offset' min=0 precision=8");
    TwAddSeparator(bar_main, NULL, " group='config'");
    
    // config (temp)
    tw_util::AddVarCB<int>(bar_main, "undo_buffer_size",
        [&](const int& value){
            memento.init(value);
        }, [&](int& value) {
            value = memento.buffer_size();
        }, "group=config label='undo buf size' min=10 max=1000");
    tw_util::AddVarCB_default<int>(bar_main, "normalSmoothIter", basemesh.normalSmoothIter,
        [&](){
            basemesh.smooth_normal();
        }, "group=config label='norm smth iter' min=0 max=50");
    TwAddVarRW(bar_main, "autocmpl_preview_interval", TW_TYPE_INT32  , &configTemp.autocmpl_preview_interval, "group=config label='autocmpl preview intrv' min=100");
    TwAddVarRW(bar_main, "cylinder_num_div", TW_TYPE_INT32  , &configTemp.cylinder_num_div, "group=config label='cylndr num div' min=1");
    TwAddVarRW(bar_main, "use_even_num_subdiv", TW_TYPE_BOOLCPP, &curvenetwork::Patch::use_even_num_subdiv, "group=config label='use even subdiv'");
    TwAddVarRW(bar_main, "prefer_rect_proc3", TW_TYPE_BOOLCPP, &curvenetwork::Patch::prefer_rect_proc3, "group=config label='prefer proc3'");
    TwAddVarRW(bar_main, "loop_threshold_area"  , TW_TYPE_DOUBLE, &configTemp.loop_threshold_area  , "group=config label='loop thr area' min=0 max=1 step=0.01");
    TwAddVarRW(bar_main, "loop_threshold_normal", TW_TYPE_DOUBLE, &configTemp.loop_threshold_normal, "group=config label='loop thr normal' min=0 max=1 step=0.01");
    
    // snakes
    tw_util::AddVarCB_default<int>(bar_main, "featureSmoothIter", basemesh.featureSmoothIter,
        [&](){
            basemesh.compute_snake_energy();
        }, "group=snakes label='feat smth iter' min=0 max=50");
    tw_util::AddVarCB_default<BaseMesh::SnakeFeatureType>(bar_main, "featureType", "FeatureType", basemesh.snake_feature_type,
        [&]() {
            basemesh.compute_snake_energy();
        }, "group=snakes label='feature type' enum='0 {mean curvature}, 1{normal var}, 2{normal curv}'");
    tw_util::AddVarCB_default<bool>(bar_main, "use_feature_as_energy", basemesh.use_feature_as_energy,
        [&](){
            basemesh.compute_snake_energy();
        }, "group=snakes label='feature as energy'");
    tw_util::AddVarCB_default<bool>(bar_main, "negate_energy", basemesh.negate_energy,
        [&](){
            basemesh.compute_snake_energy();
        }, "group=snakes label='negate energy'");
    TwAddVarRW(bar_main, "snakes_internal1", TW_TYPE_DOUBLE  , &configTemp.snakes_internal1, "group=snakes label='internal 1st' min=0");
    TwAddVarRW(bar_main, "snakes_internal2", TW_TYPE_DOUBLE  , &configTemp.snakes_internal2, "group=snakes label='internal 2nd' min=0");
    TwAddVarRW(bar_main, "snakes_external" , TW_TYPE_DOUBLE  , &configTemp.snakes_external , "group=snakes label='external'  min=0");
    TwAddVarRW(bar_main, "snakes_damping"  , TW_TYPE_DOUBLE  , &configTemp.snakes_damping  , "group=snakes label='damping'  min=0 max=1 step=0.1");
    
    //TwAddVarRW(bar_main, "weightMode", TwDefineEnum("WeightModeType", 0, 0), &stateDeformCurve.weightMode, "group=Config label='Crv Dfm Weight' enum='0 {Linear}, 1 {Cup}, 2 {Spike}, 3 {Cosine}'");
    
    // Camera
    tw_util::AddVarCB<bool>(bar_main, "camera_upright",
        [&](const bool& value) {
            Vector3d eye    = camera->get_eye();
            Vector3d center = camera->center;
            Vector3d up     = camera->get_up();
            if (value)
                camera = &camera_upright;
            else
                camera = &camera_free;
            camera->init(eye, center, up);
            camera->update_center(center);
        }, [&](bool& value) {
            value = camera == &camera_upright;
        }, "group=camera label=upright");
    TwAddVarRW(bar_main, "use_ortho"  , TW_TYPE_BOOLCPP, &configRender.use_ortho  , "group=camera label=ortho");
    TwAddVarRW(bar_main, "auto_camera_center" , TW_TYPE_BOOLCPP, &configRender.auto_camera_center, "group=camera label='auto center'");
    TwAddVarRW(bar_main, "turntable_speed", TW_TYPE_DOUBLE, &configRender.turntable_speed, "group=camera label='turntable speed' max=1 min=-1 step=0.001");
    
    // Render
    //TwAddVarRW(bar_main, "render_quad_only", TW_TYPE_BOOLCPP, &configRender.render_quad_only, "group=render label='quad only'");
    //TwAddVarRW(bar_main, "render_basemesh_only", TW_TYPE_BOOLCPP, &configRender.render_basemesh_only, "group=render label='basemesh only'");
    TwAddVarRW(bar_main, "render_mode", TwDefineEnum("render_mode_type", 0, 0), &configRender.mode, "group = render label='mode' enum='"
        "0 {default}           , "
        "1 {no basemesh}       , "
        "2 {basemesh only}     ' ");
    TwAddVarRW(bar_main, "always_show_quads", TW_TYPE_BOOLCPP, &configRender.always_show_quads, "group=render label='always show quads'");
    TwAddVarRW(bar_main, "basemesh_render_line", TW_TYPE_BOOLCPP, &configRender.basemesh_render_line, "group=render label='basemesh line'");
    TwAddVarRW(bar_main, "show_axis", TW_TYPE_BOOLCPP, &configRender.show_axis, "group=render label=axis");
    TwAddVarRW(bar_main, "show_bbox", TW_TYPE_BOOLCPP, &configRender.show_bbox, "group=render label=bbox");
    tw_util::AddVarCB_default<bool  >(bar_main, "phongTessellationEnabled", phong_tessellation::internal::enabled(),
        [&](){
            for (auto& p : curvenetwork.patches)
                p.invalidate_displist();
        }, "group=render label='phong tessell'");
    tw_util::AddVarCB_default<int   >(bar_main, "phongTessellationSubdiv" , phong_tessellation::subdiv(),
        [&](){
            for (auto& p : curvenetwork.patches)
                p.invalidate_displist();
        }, "group=render label='phong subdiv' min=1 max=6");
    tw_util::AddVarCB_default<double>(bar_main, "phongTessellationWeight" , phong_tessellation::weight(),
        [&](){
            for (auto& p : curvenetwork.patches)
                p.invalidate_displist();
        }, "group=render label='phong weight' min=0 max=1 step=0.01");
    TwAddVarRW(bar_main, "show_expmap", TW_TYPE_BOOLCPP, &configRender.show_expmap, "group=render");
    TwAddVarRW(bar_main, "overlay_alpha", TW_TYPE_FLOAT, &configRender.overlay_alpha, "group=render min=0 max=1 step=0.01");
    TwAddVarRW(bar_main, "overlay_v_flipped", TW_TYPE_BOOLCPP, &configRender.overlay_v_flipped, "group=render");
    TwAddVarRW(bar_main, "show_snake_energy", TW_TYPE_BOOLCPP, &configRender.show_snake_energy, "group=render label='show snake energy'");
    tw_util::AddVarCB_default<double>(bar_main, "snake_energy_scale", configRender.snake_energy_scale,
        [&](){
            basemesh.invalidate_displist();
        }, "group=render label='snake energy scale' min=0 step=0.1");
    TwAddVarRW(bar_main, "bgcolor_bottom", TW_TYPE_COLOR3F, &configRender.bgcolor_bottom, "group=render label='bg color btm'");
    TwAddVarRW(bar_main, "bgcolor_top"   , TW_TYPE_COLOR3F, &configRender.bgcolor_top   , "group=render label='bg color top'");
    tw_util::AddVarCB_default<double>(bar_main, "corner_radius" , configRender.corner_radius,
        [&](){
            curvenetwork.invalidate_displist();
        }, "group=render label='corner radius' min=0 step=0.001");
    TwAddVarRW(bar_main, "edge_width_interior", TW_TYPE_DOUBLE , &configRender.edge_width_interior, "group=render label='edge width i' min=0 step=0.1");
    TwAddVarRW(bar_main, "edge_width_boundary", TW_TYPE_DOUBLE , &configRender.edge_width_boundary, "group=render label='edge width b' min=0 step=0.1");
    TwAddVarRW(bar_main, "singularity_radius" , TW_TYPE_DOUBLE , &configRender.singularity_radius, "group=render label='singularity radius' min=0 step=0.1");
    TwAddVarRW(bar_main, "light_pos", TW_TYPE_DIR3F, &configRender.light_pos[0], "group=render label='light pos'");
    TwAddVarRW(bar_main, "show_auxiliary", TW_TYPE_BOOLCPP, &configRender.quadmesh_only_show_auxiliary, "group=render label='show auxiliary'");
    TwAddVarRW(bar_main, "show_left_side", TW_TYPE_BOOLCPP, &configRender.quadmesh_only_show_left_side, "group=render label='show left side'");
}
