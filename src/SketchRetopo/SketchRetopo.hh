#pragma once

// libraries
#include <kt84/geometry/CameraFree.hh>
#include <kt84/geometry/CameraUpright.hh>
#include <kt84/graphics/ProgramObject.hh>
#include <kt84/graphics/TextureObjectT.hh>
#include <rtcore/common/accel.h>
#include <rtcore/common/intersector.h>
#include <AntTweakBar.h>
#include <geodesic/geodesic_algorithm_exact.h>

// core data structure
#include "BaseMesh.hh"
#include "curvenetwork/Core.hh"
#include "curvenetwork/Vertex.hh"
#include "curvenetwork/Halfchain.hh"
#include "curvenetwork/Edgechain.hh"
#include "curvenetwork/Patch.hh"
#include "AutocmplInfo.hh"
#include <kt84/Memento.hh>

// config
#include "ConfigSaved.hh"
#include "ConfigTemp.hh"
#include "ConfigRender.hh"

#include "state/StateSketch.hh"            // edit at curve network level
#include "state/StateSpine.hh"
#include "state/StateAutocmpl.hh"
#include "state/StateLaser.hh"
#include "state/StateCylinder.hh"
#include "state/StateDeformCurve.hh"
#include "state/StateEditCorner.hh"
#include "state/StateEditTopology.hh"        // edit at quad mesh level
#include "state/StateEdgeLoop.hh"
#include "state/StateMoveVertex.hh"

struct SketchRetopo {
    static SketchRetopo& get_instance() {
        static SketchRetopo instance;
        return instance;
    }
    
    // ctor for initializing primitive member data
    SketchRetopo();
    
    // top level init (executed only once)
    void init();
    
    // core data structure
    BaseMesh     basemesh;
    curvenetwork::Core curvenetwork;
    
    // patch generation is implemented here (instead of inside CurveNetwork) because patch interior vertices need to be projected onto basemesh
    bool is_loop_ccw(const kt84::Polyline_PointNormal& loop) const;
    void generate_patches();
    void generate_patch_debug();        // to debug patch generation code
    
    // delegates to curvenetwork::Core::add_curve_open(), takes care of snapping endpoints and projection
    void add_curve_open_delegate(kt84::Polyline_PointNormal& curve, bool corner_at_openend = true);
    
    // compute interior vertices position/normal by Laplacian smoothing
    void compute_patch_interior_pn(curvenetwork::Patch* patch);
    
    // compute interior vertices position/normal by harmonic parameterization. The user must specify (click) one base mesh face as a seed to extract the mesh subset using floodfill.
    void compute_patch_interior_pn_harmonic(curvenetwork::Patch* patch, BaseMesh::FHandle submesh_seed);
    
    // used by suggestive algorithms (spine sketching, auto-completion)
    void add_uv_line(const Eigen::Vector2d& uv0, const Eigen::Vector2d& uv1, const boost::optional<Eigen::Vector2d>& unnormalized_dir0 = boost::none, const boost::optional<Eigen::Vector2d>& unnormalized_dir1 = boost::none);
    
    // spine sketching
    void add_spine_curve(const kt84::Polyline_PointNormal& spine_curve);
    
    // patch auto-completion
    void compute_autocompl(const AutocmplInfo& autocmpl_info);
    std::vector<AutocmplInfo> suggest_autocmpl(curvenetwork::Halfedge* h, bool favor_triangle);
    
    // cylinder sketching
    void add_cylinder_curve(kt84::Polyline_PointNormal cylinder_curve);     // by-value argument is intentional here...
    
    // extract a curve from cutting plane (laser sketching)
    kt84::Polyline_PointNormal trace_on_plane(const kt84::PointNormal& pn_start, const Eigen::Vector3d& plane_normal, double dist_max);
    
    // edge loop insertion
    void                       edgeLoop_insert (curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start, double t);
    kt84::Polyline_PointNormal edgeLoop_preview(curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start, double t) const;
    
    // walking over edge loop, used by EditSubdiv mode
    std::vector<std::pair<curvenetwork::Patch*, curvenetwork::Patch::HHandle>> edgeLoop_walk(curvenetwork::Patch* patch_start, curvenetwork::Patch::HHandle h_start) const;
    
    // vertex moving/smoothing
    void vertices_move  (const kt84::PointNormal& pn_center, const Eigen::Vector3d& center_offset, double radius);
    void vertices_smooth(const kt84::PointNormal& pn_center, double radius);
    auto vertices_smooth_sub(curvenetwork::Patch* patch, curvenetwork::Patch::VHandle patch_v, std::vector<curvenetwork::Edgechain*>& affected_edgechains) -> boost::optional<kt84::PointNormal> const;
    void vertices_smooth_global();
    void set_halfchain_pn_from_patch(curvenetwork::Halfchain* c);         // inverse of Patch::set_boundary_pn(): set PointNormal of its surrounding curve network vertices
    
    // simple geometric snakes
    void snakes_move(curvenetwork::Edgechain* e);
    void snakes_move(curvenetwork::Vertex   * v);
    
    // toggle symmetry mode on/off
    void toggle_symmetric();
    
    // when basemesh is loaded, trace mesh boundary edges (if any) and add them to curve network
    void trace_basemesh_boundary();
    
    // geodesic
    geodesic::Mesh                   geodesic_mesh;
    geodesic::GeodesicAlgorithmExact geodesic_algorithm;
    kt84::Polyline_PointNormal       geodesic_compute(const kt84::PointNormal& pn0, const kt84::PointNormal& pn1);
    void                             geodesic_init();
    
    // camera
    kt84::Camera*       camera;
    kt84::CameraFree    camera_free;
    kt84::CameraUpright camera_upright;
    
    // config data
    ConfigRender configRender;
    ConfigSaved  configSaved;
    ConfigTemp   configTemp;
    
    // Intel Embree ray tracer used for projection
    embree::Ref<embree::Accel      > embree_accel;
    embree::Ref<embree::Intersector> embree_intersector;
    void                             embree_init();
    
    // interection/projection
    boost::optional<kt84::PointNormal> uv_to_pn(const Eigen::Vector2d  & uv) const;
    boost::optional<Eigen::Vector2d  > pn_to_uv(const kt84::PointNormal& pn) const;
    boost::optional<kt84::PointNormal> intersect_convert(const embree::Hit& hit) const;
    embree::Hit                        intersect(const kt84::PointNormal& pn_in) const;
    embree::Hit                        intersect(int mouse_x, int mouse_y) const;
    void project(kt84::PointNormal& pn) const;
    void project(kt84::Polyline_PointNormal& polyline) const { for (auto& pn : polyline) project(pn); }
    void reproject_all();       // projects all vertices again (i.e., curvenetwork::Vertex::pn, curvenetwork::Patch::Vertex::pn)
    
    // common behavior. if true is returned, event is not passed to the current state.
    enum class CommonDragMode {
        None,
        Smooth,
        Snakes
    }                                     common_dragMode;
    std::vector<curvenetwork::Edgechain*> common_affected_edgechains;
    boost::optional<kt84::PointNormal>    common_pn_mouse;
    Eigen::Vector2i                       common_mouse_pos;
    Eigen::Vector2i                       common_mouse_pos_prev;
    kt84::Polyline_PointNormal            common_sketch_curve;
    bool                                  common_key_pressed[256];
    bool                                  common_button_left_pressed;
    bool                                  common_button_right_pressed;
    bool                                  common_button_middle_pressed;
    bool                                  common_shift_pressed;             // to know modifier key state during motion event
    bool                                  common_ctrl_pressed;
    bool                                  common_alt_pressed;
    void display_pre ();
    void display_post();
    bool common_mouse_down(int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    bool common_mouse_up  (int mouse_x, int mouse_y, Button button, bool shift_pressed, bool ctrl_pressed, bool alt_pressed);
    bool common_mouse_move(int mouse_x, int mouse_y);
    bool common_keyboard  (unsigned char key, int x, int y);
    bool common_keyboardup(unsigned char key, int x, int y);
    
    kt84::Polyline2d hide_stroke;
    void hide_mesh_by_stroke();
    
    void render_default();
    void render_quadmesh_only();
    void render_basemesh_only();
    
    // state
    State* state;
    StateSketch       stateSketch      ;
    StateSpine        stateSpine       ;
    StateAutocmpl     stateAutocmpl    ;
    StateLaser        stateLaser       ;
    StateCylinder     stateCylinder    ;
    StateDeformCurve  stateDeformCurve ;
    StateEditTopology stateEditTopology;
    StateEditCorner   stateEditCorner  ;
    StateEdgeLoop     stateEdgeLoop    ;
    StateMoveVertex   stateMoveVertex  ;
    enum class EnumState {
        Sketch      = 0,
        Spine       ,
        Autocmpl    ,
        Laser       ,
        Cylinder    ,
        DeformCurve ,
        EditCorner  ,
        EditTopology,
        EdgeLoop    ,
        MoveVertex  ,
    };
    void      state_set(EnumState enumState);
    EnumState state_get() const;
    
    // undo-redo with "Memento" design pattern
    kt84::Memento<curvenetwork::Core> memento;
    void memento_store();
    void memento_undo();
    void memento_redo();
    
    // shaders
    kt84::ProgramObject program_basemesh;
    kt84::ProgramObject program_basic;
    kt84::ProgramObject program_corner;
    void                program_init();
    int stencil_index;
    
    // material texture png file (encoded using base64)
    kt84::TextureObject texture_overlay;
    std::string         material_texture_png_base64;
    kt84::TextureObject material_texture;
    bool                material_texture_load(const std::string& fname);
    bool                texture_overlay_load(const std::string& fname);
    void                material_texture_update();
    
    // AntTweakBar
    TwBar* bar_main;
    void   bar_init();
    
    // file IO
    bool import_basemesh(const std::string& fname);
    bool import_autopomesh(const std::string& fname);
    bool export_basemesh(std::string fname) const;
    bool export_retopomesh(std::string fname) const;
    bool xml_load(const std::string& fname);
    bool xml_save(std::string fname);
    bool xml_load_batch(const std::string& fname);
    // XML file format version
    int get_version() const;
};

