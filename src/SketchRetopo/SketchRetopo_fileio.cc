#include <OpenMesh/Core/IO/MeshIO.hh>
#include "SketchRetopo.hh"
#include <tinyxml2/tinyxml2.h>
#include <lodepng/lodepng.h>
#include <base64/encode.h>
#include <kt84/Clock.hh>
#include <kt84/openmesh/merge_nearby_vertices.hh>
#include <kt84/zenity_util.hh>
#include <kt84/graphics/graphics_util.hh>
#include <kt84/safe_istream.hh>
#include <patchgen/decl.hh>
#include "AutopoMesh.hh"
#include "helper.hh"
#include "curvenetwork/helper.hh"
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

bool SketchRetopo::texture_overlay_load(const string& fname) {
    vector<unsigned char> pixels;
    unsigned int width, height;
    if (lodepng::decode(pixels, width, height, fname)) {
        cerr << "Error while reading " << fname << endl;
        return false;
    }

    int max_texture_size = glGet1i(ParamName::MAX_TEXTURE_SIZE);
    if (max_texture_size < width || max_texture_size < height) {
        cerr << "Requested resolution of " << width << "x" << height << " exceeds the maximum of " << max_texture_size << "x" << max_texture_size << endl;
        return false;
    }
    
    texture_overlay.bind();
    texture_overlay.allocate(width, height);
    texture_overlay.copy_cpu2gpu(GL_UNSIGNED_BYTE, pixels.data());
    texture_overlay.unbind();
    configRender.overlay_alpha = 0.5;

    return true;
}
bool SketchRetopo::material_texture_load(const string& fname) {
    ifstream fin(fname, ios_base::binary);
    if (!fin) return false;
    
    stringstream ss_encoded;
    base64::encoder E;
    E.encode(fin, ss_encoded);
    material_texture_png_base64 = ss_encoded.str();
    
    material_texture_update();
    
    return true;
}

bool SketchRetopo::import_basemesh(const string& fname) {
    {
        ClkSimple clk("OpenMesh::IO::read_mesh");
        BaseMesh basemesh_temp;
        OpenMesh::IO::Options opt;
        opt += OpenMesh::IO::Options::FaceTexCoord;
        if (!OpenMesh::IO::read_mesh(basemesh_temp, fname, opt))
            return false;
        basemesh_temp.has_texture = (opt & OpenMesh::IO::Options::FaceTexCoord) != 0;
        basemesh = basemesh_temp;
    }
    
    {
        ClkSimple clk("BaseMesh::init");
        basemesh.init();
    }
    
    double bb_diagonal = basemesh.boundingBox_diagonal_norm();
    
    curvenetwork.clear();
    
    // clear undo/redo data
    memento.init();
    
    // reset all config to default
    configRender = ConfigRender();
    configSaved  = ConfigSaved ();
    configTemp   = ConfigTemp  ();
    
    // parameters dependent on model size
    configTemp.quadSize   .value = bb_diagonal * ConfigTemp::default_quadSize_ratio;
    configTemp.segmentSize.value = bb_diagonal * ConfigTemp::default_segmentSize_ratio;
    curvenetwork::Patch::quadSize = configTemp.quadSize.value;
    
    configSaved.projOffset = 0.0005 * bb_diagonal;
    double step_projOffset = 0.1 * configSaved.projOffset;
    TwSetParam(bar_main, "projOffset", "step", TW_PARAM_DOUBLE, 1, &step_projOffset);
    
    // init camera
    camera->init(basemesh.centerOfMass + Vector3d(0, 0, bb_diagonal), basemesh.centerOfMass, Vector3d::UnitY());
    
    {   // init embree
        ClkSimple clk("embree_init");
        embree_init();
    }
    
    {   // init geidesic
        ClkSimple clk("geidesic_init");
        geodesic_init();
    }
    
    state->init();
    
    trace_basemesh_boundary();
    
    return true;
}

bool SketchRetopo::import_autopomesh(const std::string& fname) {
    AutopoMesh autopomesh;
    {
        ClkSimple clk("OpenMesh::IO::read_mesh");
        if (!OpenMesh::IO::read_mesh(autopomesh, fname))
            return false;
    }
    autopomesh.init();
    
    ClkSimple clk("constructing curve network topology...");
    curvenetwork.clear();
    
    // set is_separatrix flag for autopomesh edges
    autopomesh.trace_separatrix();

    // am: prefix for autopomesh elements
    // cn: prefix for curve network elements
    
    // add curve network vertices
    for (auto am_v : autopomesh.vertices()) {
        if (!autopomesh.on_separatrix(am_v))
            continue;
        
        auto cn_v = curvenetwork.new_vertex();
        cn_v->pn <<
            o2e(autopomesh.point (am_v)),
            o2e(autopomesh.normal(am_v));
        project(cn_v->pn);
        
        autopomesh.data(am_v).cn_vertex = cn_v;
    }
    
    // add curve network halfedges
    for (auto am_e : autopomesh.edges()) {
        if (!autopomesh.data(am_e).on_separatrix)
            continue;
        
        auto am_h0 = autopomesh.halfedge_handle(am_e, 0);
        auto am_h1 = autopomesh.halfedge_handle(am_e, 1);
        auto am_v0 = autopomesh.from_vertex_handle(am_h0);
        auto am_v1 = autopomesh.from_vertex_handle(am_h1);
        
        auto cn_h0 = curvenetwork.new_halfedge();
        auto cn_h1 = curvenetwork.new_halfedge();
        auto cn_v0 = autopomesh.data(am_v0).cn_vertex;
        auto cn_v1 = autopomesh.data(am_v1).cn_vertex;
        
        curvenetwork::set_vertex_halfedge(cn_v0, cn_v1, cn_h0, cn_h1);
        
        autopomesh.data(am_h0).cn_halfedge = cn_h0;
        autopomesh.data(am_h1).cn_halfedge = cn_h1;
    }
    
    // set halfedge next/prev
    for (auto am_e : autopomesh.edges()) {
        if (!autopomesh.data(am_e).on_separatrix)
            continue;
        
        for (int i = 0; i < 2; ++i) {
            auto am_h = autopomesh.halfedge_handle(am_e, i);
            auto cn_h = autopomesh.data(am_h     ).cn_halfedge;
            
            cn_h->is_corner = autopomesh.is_corner(autopomesh.to_vertex_handle(am_h));
            
            auto am_h_next = cn_h->is_corner ? autopomesh.next_halfedge_handle(am_h) : autopomesh.edgeflow_next(am_h);
            auto cn_h_next = autopomesh.data(am_h_next).cn_halfedge;
            
            curvenetwork::set_next_prev(cn_h, cn_h_next);
            
        }
    }
    
    curvenetwork.generate_halfchains();
    
    // set curvenetwork::Edgechain::num_subdiv
    for (auto& cn_e : curvenetwork.edgechains) {
        cn_e.num_subdiv = 1;
        auto cn_c = cn_e.halfchain[0];
        for (auto cn_h = cn_c->halfedge_front; cn_h != cn_c->halfedge_back; cn_h = cn_h->next)
            ++cn_e.num_subdiv;
    }
    
    // now ready to start constructing patches!
    for (auto am_h : autopomesh.halfedges()) {
        if (autopomesh.data(am_h).is_processed || autopomesh.is_boundary(am_h))
            continue;
        
        // process if its from-vertex is patch corner
        auto am_v_corner0 = autopomesh.from_vertex_handle(am_h);
        if (!autopomesh.is_corner(am_v_corner0))
            continue;
        
        // set processed flag
        for (auto am_h2 = am_h; ;) {
            if (autopomesh.is_corner(autopomesh.to_vertex_handle(am_h2))) {
                am_h2 = autopomesh.next_halfedge_handle(am_h2);
                autopomesh.data(am_h2).is_processed = true;
            
            } else {
                am_h2 = autopomesh.edgeflow_next(am_h2);
            }
            
            if (am_h2 == am_h)
                break;
        }
        
        
        // create new patch
        auto patch = curvenetwork.new_patch();
        patch->halfchain = autopomesh.data(am_h).cn_halfedge->halfchain;
        for (auto c = patch->halfchain; ; ) {
            c->patch = patch;
            c = c->next();
            if (c == patch->halfchain)
                break;
        }
        
        // determine patch width
        int w = 1;
        for (auto am_h2 = am_h; ; ++w) {
            if (autopomesh.is_corner(autopomesh.to_vertex_handle(am_h2)))
                break;
            
            am_h2 = autopomesh.edgeflow_next(am_h2);
        }
        
        // determine patch height
        int h = 1;
        for (auto am_h2 = autopomesh.prev_halfedge_handle(am_h); ; ++h) {
            if (autopomesh.is_corner(autopomesh.from_vertex_handle(am_h2)))
                break;
            
            am_h2 = autopomesh.edgeflow_prev(am_h2);
        }
        
        // add vertices to the patch
        auto am_h_row_start = am_h;
        for (int iy = 0; iy <= h; ++iy) {
            auto am_h_row = am_h_row_start;
            
            for (int ix = 0; ix <= w; ++ix) {
                auto am_v = autopomesh.from_vertex_handle(am_h_row);
                
                auto patch_v = patch->add_vertex(OpenMesh::Vec3d());
                
                // copy point & normal from autopo mesh
                auto& patch_v_pn = patch->data(patch_v).pn;
                patch_v_pn <<
                    o2e(autopomesh.point (am_v)),
                    o2e(autopomesh.normal(am_v));
                project(patch_v_pn);
                
                // set corner flag appropriately
                patch->data(patch_v).patchgen.corner_index =
                    ix == 0 && iy == 0 ? 0 :
                    ix == w && iy == 0 ? 1 :
                    ix == w && iy == h ? 2 :
                    ix == 0 && iy == h ? 3 : - 1;
                
                // next halfedge on the same row
                am_h_row = autopomesh.edgeflow_next(am_h_row);
            }
            
            // move down to the next row
            am_h_row_start =
                autopomesh.opposite_halfedge_handle(
                autopomesh.prev_halfedge_handle    (
                autopomesh.prev_halfedge_handle    (am_h_row_start)));
        }
        
        // add faces to the patch
        for (int iy = 0; iy < h; ++iy)
        for (int ix = 0; ix < w; ++ix)
        {
            auto patch_v0 = patch->vertex_handle(ix     + (w + 1) *  iy     );
            auto patch_v1 = patch->vertex_handle(ix + 1 + (w + 1) *  iy     );
            auto patch_v2 = patch->vertex_handle(ix + 1 + (w + 1) * (iy + 1));
            auto patch_v3 = patch->vertex_handle(ix     + (w + 1) * (iy + 1));
            patch->add_face(patch_v0, patch_v1, patch_v2, patch_v3);
        }
        
#ifndef NDEBUG
        patch->debugInfo_get(false, false, true);
#endif
        patch->set_halfedgeData();
    }
    
    // clear undo/redo data
    memento.init();
    
    configSaved.symmetric = false;
    
    return true;
    
}

bool SketchRetopo::export_basemesh(std::string fname) const {
    return OpenMesh::IO::write_mesh(basemesh, fname);
}

bool SketchRetopo::export_retopomesh(string fname) const {
    struct RetopoMeshTraits : public OpenMesh::DefaultTraits {
        typedef OpenMesh::Vec3d Point;
        typedef OpenMesh::Vec3d Normal;
    };
    typedef OpenMesh::PolyMesh_ArrayKernelT<RetopoMeshTraits> RetopoMeshBase;
    struct RetopoMesh
        : public RetopoMeshBase
#ifndef NDEBUG
        , public DebugInfo<RetopoMeshBase, RetopoMesh>
#endif
    {} retopoMesh;
    
    double threshold = basemesh.boundingBox_diagonal_norm() * 0.0000001;
    
    // concatenate all patches
    for (auto& patch : curvenetwork.patches) {
        vector<RetopoMesh::VHandle> added_vertices;
        for (auto v : patch.vertices()) {
            RetopoMesh::Point point(&patch.data(v).pn[0]);
            
            //if (configSaved.symmetric && abs(point[0]) < threshold)
            //    point[0] = 0;
            
            added_vertices.push_back(retopoMesh.add_vertex(point));
        }
        for (auto f : patch.faces()) {
            vector<RetopoMesh::VHandle> face_vhandles;
            for (auto v : patch.fv_range(f))
                face_vhandles.push_back(added_vertices[v.idx()]);
            retopoMesh.add_face(face_vhandles);
        }
    }
    
    if (configSaved.symmetric) {
        // symmetric part
        for (auto& patch : curvenetwork.patches) {
            vector<RetopoMesh::VHandle> added_vertices;
            for (auto v : patch.vertices()) {
                RetopoMesh::Point point(&patch.data(v).pn[0]);
                
                //if (abs(point[0]) < threshold)
                //    point[0] = 0;
                
                point[0] *= -1;
                added_vertices.push_back(retopoMesh.add_vertex(point));
            }
            for (auto f : patch.faces()) {
                vector<RetopoMesh::VHandle> face_vhandles;
                for (auto v : patch.fv_range(f))
                    face_vhandles.push_back(added_vertices[v.idx()]);
                reverse(face_vhandles.begin(), face_vhandles.end());
                retopoMesh.add_face(face_vhandles);
            }
        }
    }
    
#ifndef NDEBUG
    retopoMesh.debugInfo_get(false, false, true);
#endif
    kt84::merge_nearby_vertices(retopoMesh, threshold);
    
    return OpenMesh::IO::write_mesh(retopoMesh, fname);
}

bool SketchRetopo::xml_load(const string& fname) {
    tinyxml2::XMLDocument xml_doc;
    if (xml_doc.LoadFile(fname.c_str()))
        return false;
    
    auto xml_root = xml_doc.FirstChildElement("sketchretopo");
    if (!xml_root)
        return false;
    int version = xml_root->IntAttribute("version");
    if (version != get_version())
        cout << "The XML file format version is older than the current one (" << get_version() << ")!\n";
    
    // base mesh
    BaseMesh basemesh_temp;
    auto xml_basemesh = xml_root->FirstChildElement("basemesh");
    if (!xml_basemesh)
        return false;
    {
        basemesh_temp.has_texture = xml_basemesh->IntAttribute("has_texture");
        
        // vertices
        auto xml_basemesh_vertices = xml_basemesh->FirstChildElement("vertices");
        int num_vertices = xml_basemesh_vertices->IntAttribute("num");
        istringstream ss_basemesh_vertices;
        ss_basemesh_vertices.str(xml_basemesh_vertices->GetText());
        for (int i = 0; i < num_vertices; ++i) {
            BaseMesh::Point p;
            make_safe_istream(ss_basemesh_vertices) >> p[0] >> p[1] >> p[2];
            basemesh_temp.add_vertex(p);
        }
        
        // faces
        auto xml_basemesh_faces = xml_basemesh->FirstChildElement("faces");
        int num_faces = xml_basemesh_faces->IntAttribute("num");
        istringstream ss_basemesh_faces;
        ss_basemesh_faces.str(xml_basemesh_faces->GetText());
        for (int i = 0; i < num_faces; ++i) {
            int vidx0, vidx1, vidx2;
            make_safe_istream(ss_basemesh_faces) >> vidx0 >> vidx1 >> vidx2;
            auto v0 = basemesh_temp.vertex_handle(vidx0);
            auto v1 = basemesh_temp.vertex_handle(vidx1);
            auto v2 = basemesh_temp.vertex_handle(vidx2);
            auto f = basemesh_temp.add_face(v0, v1, v2);
            if (f.is_valid() && basemesh_temp.has_texture) {
                for (auto h : basemesh_temp.fh_range(f)) {
                    BaseMesh::TexCoord2D uv;
                    make_safe_istream(ss_basemesh_faces) >> uv[0] >> uv[1];
                    basemesh_temp.set_texcoord2D(h, uv);
                }
            }
        }
    }
    
    ConfigSaved configSaved_temp;
    {   // config
        auto xml_config = xml_root->FirstChildElement("config");
        if (!xml_config)
            return false;
        configSaved_temp.projOffset    = xml_config->DoubleAttribute("projoffset");
        configSaved_temp.symmetric     = xml_config->Attribute("symmetric", "true") ? true : false;
    }
    
    // curveNetwork
    curvenetwork::Core curvenetwork_temp;
    auto xml_curveNetwork = xml_root->FirstChildElement("curvenetwork");
    if (!xml_curveNetwork)
        return false;
    
    {   // vertices
        auto xml_curveNetwork_vertices = xml_curveNetwork->FirstChildElement("vertices");
        if (!xml_curveNetwork_vertices)
            return false;
        int num_curveNetwork_vertices = xml_curveNetwork_vertices->IntAttribute("num");
        istringstream ss_curveNetwork_vertices(num_curveNetwork_vertices ? xml_curveNetwork_vertices->GetText() : "");
        for (int i = 0; i < num_curveNetwork_vertices; ++i) {
            auto v = curvenetwork_temp.new_vertex();
            make_safe_istream(ss_curveNetwork_vertices)
                >> v.id
                >> v->halfedge.id
                >> v->pn[0] >> v->pn[1] >> v->pn[2]
                >> v->pn[3] >> v->pn[4] >> v->pn[5];
            v->id = v.id;
        }
    }
    
    {   // halfedges
        auto xml_curveNetwork_halfedges = xml_curveNetwork->FirstChildElement("halfedges");
        if (!xml_curveNetwork_halfedges)
            return false;
        int num_curveNetwork_halfedges = xml_curveNetwork_halfedges->IntAttribute("num");
        istringstream ss_curveNetwork_halfedges(num_curveNetwork_halfedges ? xml_curveNetwork_halfedges->GetText() : "");
        for (int i = 0; i < num_curveNetwork_halfedges; ++i) {
            auto h = curvenetwork_temp.new_halfedge();
            int is_corner;
            int imaginary_patch_id;
            make_safe_istream(ss_curveNetwork_halfedges)
                >> h.id
                >> h->vertex   .id
                >> h->next     .id
                >> h->prev     .id
                >> h->opposite .id
                >> h->halfchain.id
                >> is_corner
                >> imaginary_patch_id;
            h->is_corner = is_corner != 0;
            if (imaginary_patch_id == -2) h->imaginary_patch = &curvenetwork::Patch::imaginary_patch_symmetry;
            if (imaginary_patch_id == -3) h->imaginary_patch = &curvenetwork::Patch::imaginary_patch_boundary;
            h->id = h.id;
        }
    }
    
    {   // halfchains
        auto xml_curveNetwork_halfchains = xml_curveNetwork->FirstChildElement("halfchains");
        if (!xml_curveNetwork_halfchains)
            return false;
        int num_curveNetwork_halfchains = xml_curveNetwork_halfchains->IntAttribute("num");
        istringstream ss_curveNetwork_halfchains(num_curveNetwork_halfchains ? xml_curveNetwork_halfchains->GetText() : "");
        for (int i = 0; i < num_curveNetwork_halfchains; ++i) {
            auto c = curvenetwork_temp.new_halfchain();
            make_safe_istream(ss_curveNetwork_halfchains)
                >> c.id
                >> c->halfedge_front.id
                >> c->halfedge_back .id
                >> c->patch         .id
                >> c->edgechain     .id;
            c->id = c.id;
        }
    }
    
    {   // edgechains
        auto xml_curveNetwork_edgechains = xml_curveNetwork->FirstChildElement("edgechains");
        if (!xml_curveNetwork_edgechains)
            return false;
        int num_curveNetwork_edgechains = xml_curveNetwork_edgechains->IntAttribute("num");
        istringstream ss_curveNetwork_edgechains(num_curveNetwork_edgechains ? xml_curveNetwork_edgechains->GetText() : "");
        for (int i = 0; i < num_curveNetwork_edgechains; ++i) {
            auto e = curvenetwork_temp.new_edgechain();
            make_safe_istream(ss_curveNetwork_edgechains)
                >> e.id
                >> e->halfchain[0].id
                >> e->halfchain[1].id
                >> e->num_subdiv;
            e->id = e.id;
        }
    }
    
    {   // patches
        auto xml_curveNetwork_patches = xml_curveNetwork->FirstChildElement("patches");
        if (!xml_curveNetwork_patches)
            return false;
        int num_curveNetwork_patches = xml_curveNetwork_patches->IntAttribute("num");
        
        for (auto xml_curveNetwork_patch = xml_curveNetwork_patches->FirstChildElement("patch");
            xml_curveNetwork_patch;
            xml_curveNetwork_patch = xml_curveNetwork_patch->NextSiblingElement("patch"))
        {
            auto patch = curvenetwork_temp.new_patch();
            patch->id = patch.id = xml_curveNetwork_patch->IntAttribute("id"        );
            patch->halfchain.id  = xml_curveNetwork_patch->IntAttribute("halfchain" );
            patch->is_failure    = xml_curveNetwork_patch->IntAttribute("is_failure");
            
            {   // patch parameter
                auto& param = patch->param;
                auto xml_patch_param = xml_curveNetwork_patch->FirstChildElement("param");
                if (version >= 1) {
                    assert(xml_patch_param);
                    int num_sides    = xml_patch_param->IntAttribute("num_sides");
                    param.pattern_id = xml_patch_param->IntAttribute("pattern_id");
                    // permutation
                    int permutation_id = xml_patch_param->IntAttribute("permutation_id");
                    param.permutation.init(num_sides, permutation_id);
                    // l, variable
                    param.l = param.p = param.q = VectorXi::Constant(num_sides, -1);
                    auto xml_param_l        = xml_patch_param->FirstChildElement("l"       );
                    auto xml_param_variable = xml_patch_param->FirstChildElement("variable");
                    istringstream ss_param_l;
                    istringstream ss_param_variable;
                    ss_param_l       .str(xml_param_l       ->GetText());
                    ss_param_variable.str(xml_param_variable->GetText());
                    for (int i = 0; i < num_sides; ++i)
                        make_safe_istream(ss_param_l) >> param.l[i];
                    int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
                    for (int i = 0; i < num_variables; ++i)
                        make_safe_istream(ss_param_variable) >> patchgen::get_variable(num_sides, param.pattern_id, param, i);
                }
            }
            
            {   // patch vertices
                auto xml_patch_vertices = xml_curveNetwork_patch->FirstChildElement("patch_vertices");
                if (!xml_patch_vertices)
                    return false;
                
                int num_patch_vertices = xml_patch_vertices->IntAttribute("num");
                istringstream ss_patch_vertices;
                if (num_patch_vertices)
                    ss_patch_vertices.str(xml_patch_vertices->GetText());
                
                for (int i = 0; i < num_patch_vertices; ++i) {
                    auto v = patch->add_vertex(OpenMesh::Vec3d());
                    auto& vdata = patch->data(v);
                    
                    int tag;
                    make_safe_istream(ss_patch_vertices)
                        >> vdata.pn[0] >> vdata.pn[1] >> vdata.pn[2]
                        >> vdata.pn[3] >> vdata.pn[4] >> vdata.pn[5]
                        >> vdata.patchgen.corner_index
                        >> tag;
                    vdata.patchgen.tag = static_cast<patchgen::PatchVertexTag>(tag);
                }
            }
            
            {   // patch faces
                auto xml_patch_faces = xml_curveNetwork_patch->FirstChildElement("patch_faces");
                if (!xml_patch_faces)
                    return false;
                
                int num_patch_faces = xml_patch_faces->IntAttribute("num");
                istringstream ss_patch_faces;
                if (num_patch_faces)
                    ss_patch_faces.str(xml_patch_faces->GetText());
                
                for (int i = 0; i < num_patch_faces; ++i) {
                    int vidx[4];
                    make_safe_istream(ss_patch_faces) >> vidx[0] >> vidx[1] >> vidx[2] >> vidx[3];
                    
                    curvenetwork::Patch::VHandle v[4];
                    for (int j = 0; j < 4; ++j)
                        v[j] = patch->vertex_handle(vidx[j]);
                    
                    patch->add_face(v[0], v[1], v[2], v[3]);
                }
            }
#ifndef NDEBUG
            patch->debugInfo_get(false, false, false);
#endif
        }
    }

    {   // material texture
        auto xml_texture = xml_root->FirstChildElement("texture");
        if (xml_texture) {
            material_texture_png_base64 = xml_texture->GetText();
            material_texture_update();
        }
    }
    
    // all green, replace current data
    basemesh = basemesh_temp;
    
    // init process mostly the same as in import_basemesh()
    {
        ClkSimple clk("BaseMesh::init");
        basemesh.init();
    }
    
    {   // init embree
        ClkSimple clk("embree_init");
        embree_init();
    }
    
    {   // init geidesic
        ClkSimple clk("geidesic_init");
        geodesic_init();
    }
    
    curvenetwork = curvenetwork_temp;
    curvenetwork.validate_all_ptr();
    if (version < 1) {
        // older version without patch info -> delete all patches and generate again
        for (auto& p : curvenetwork.patches)
            curvenetwork::delete_patch(&p);
        curvenetwork.garbage_collect();
        generate_patches();
    }
    // factorize direct Laplace solver
    for (auto& patch : curvenetwork.patches)
        patch.laplaceSolver_init();
    
    // clear undo/redo data
    memento.init();
    
    configSaved = configSaved_temp;
    
    // reset all config to default
    configRender = ConfigRender();
    configTemp   = ConfigTemp  ();
    
    // parameters dependent on model size
    double bb_diagonal = basemesh.boundingBox_diagonal_norm();
    configTemp.quadSize   .value = bb_diagonal * ConfigTemp::default_quadSize_ratio;
    configTemp.segmentSize.value = bb_diagonal * ConfigTemp::default_segmentSize_ratio;
    curvenetwork::Patch::quadSize = configTemp.quadSize.value;
    
    double step_projOffset = 0.1 * configSaved.projOffset;
    TwSetParam(bar_main, "projOffset", "step", TW_PARAM_DOUBLE, 1, &step_projOffset);
    
    // init camera
    camera->init(basemesh.centerOfMass + Vector3d(0, 0, bb_diagonal), basemesh.centerOfMass, Vector3d::UnitY());
    
    state->init();
    
    configTemp.autoSave.filename = fname;
    
    return true;
}

bool SketchRetopo::xml_save(string fname) {
    tinyxml2::XMLDocument xml_doc;
    
    // root element
    auto xml_root = xml_doc.InsertEndChild(xml_doc.NewElement("sketchretopo"))->ToElement();
    if (!xml_root) return false;
    xml_root->SetAttribute("version", get_version());
    
    {   // base mesh
        auto xml_basemesh = xml_root->InsertEndChild(xml_doc.NewElement("basemesh"))->ToElement();
        xml_basemesh->SetAttribute("has_texture", basemesh.has_texture ? 1 : 0);
        
        {   // vertices
            auto xml_basemesh_vertices = xml_basemesh->InsertEndChild(xml_doc.NewElement("vertices"))->ToElement();
            xml_basemesh_vertices->SetAttribute("num", static_cast<int>(basemesh.n_vertices()));
        
            ostringstream ss_basemesh_vertices;
            ss_basemesh_vertices << endl;
        
            for (auto v : basemesh.vertices()) {
                auto p = basemesh.point(v);
                ss_basemesh_vertices << p[0] << " " << p[1] << " " << p[2] << endl;
            }
            xml_basemesh_vertices->InsertEndChild(xml_doc.NewText(ss_basemesh_vertices.str().c_str()));
        }
        
        {   // faces
            auto xml_basemesh_faces = xml_basemesh->InsertEndChild(xml_doc.NewElement("faces"))->ToElement();
            xml_basemesh_faces->SetAttribute("num", static_cast<int>(basemesh.n_faces()));
        
            ostringstream ss_basemesh_faces;
            ss_basemesh_faces << endl;
        
            for (auto f : basemesh.faces()) {
                for (auto h : basemesh.fh_range(f)) {
                    auto v = basemesh.to_vertex_handle(h);
                    ss_basemesh_faces << v.idx() << " ";
                }
                if (basemesh.has_texture) {
                    for (auto h : basemesh.fh_range(f)) {
                        auto t = basemesh.texcoord2D(h);
                        ss_basemesh_faces << t[0] << " " << t[1] << " ";
                    }
                }
                ss_basemesh_faces << endl;
            }
            xml_basemesh_faces->InsertEndChild(xml_doc.NewText(ss_basemesh_faces.str().c_str()));
        }
    }
    
    {   // config
        auto xml_config = xml_root->InsertEndChild(xml_doc.NewElement("config"))->ToElement();
        xml_config->SetAttribute("projoffset", configSaved.projOffset);
        xml_config->SetAttribute("symmetric" , configSaved.symmetric ? "true" : "false");
    }
    
    // curveNetwork
    auto xml_curveNetwork = xml_root->InsertEndChild(xml_doc.NewElement("curvenetwork"))->ToElement();
    
    {   // vertices
        auto xml_curveNetwork_vertices = xml_curveNetwork->InsertEndChild(xml_doc.NewElement("vertices"))->ToElement();
        xml_curveNetwork_vertices->SetAttribute("num", static_cast<int>(curvenetwork.vertices.size()));
        ostringstream ss_curveNetwork_vertices;
        ss_curveNetwork_vertices << endl;
        
        for (auto& v : curvenetwork.vertices) {
            ss_curveNetwork_vertices
                << v.id << " "
                << v.halfedge.id << " "
                << v.pn[0] << " " << v.pn[1] << " " << v.pn[2] << " "
                << v.pn[3] << " " << v.pn[4] << " " << v.pn[5] << endl;
        }
        
        xml_curveNetwork_vertices->InsertEndChild(xml_doc.NewText(ss_curveNetwork_vertices.str().c_str()));
        xml_curveNetwork_vertices->InsertEndChild(xml_doc.NewComment("id halfedge px py pz nx ny nz"));
    }
    
    {   // halfedges
        auto xml_curveNetwork_halfedges = xml_curveNetwork->InsertEndChild(xml_doc.NewElement("halfedges"))->ToElement();
        xml_curveNetwork_halfedges->SetAttribute("num", static_cast<int>(curvenetwork.halfedges.size()));
        ostringstream ss_curveNetwork_halfedges;
        ss_curveNetwork_halfedges << endl;
        
        for (auto& h : curvenetwork.halfedges) {
            ss_curveNetwork_halfedges
                << h.id << " "
                << h.vertex   .id << " "
                << h.next     .id << " "
                << h.prev     .id << " "
                << h.opposite .id << " "
                << h.halfchain.id << " "
                << (h.is_corner ? 1 : 0) << " "
                << (h.imaginary_patch ? h.imaginary_patch->id : -1) << endl;
        }
        
        xml_curveNetwork_halfedges->InsertEndChild(xml_doc.NewText(ss_curveNetwork_halfedges.str().c_str()));
        xml_curveNetwork_halfedges->InsertEndChild(xml_doc.NewComment("id vertex next prev opposite halfchain is_corner imaginary_patch_id"));
    }
    
    {   // halfhcains
        auto xml_curveNetwork_halfchains = xml_curveNetwork->InsertEndChild(xml_doc.NewElement("halfchains"))->ToElement();
        xml_curveNetwork_halfchains->SetAttribute("num", static_cast<int>(curvenetwork.halfchains.size()));
        ostringstream ss_curveNetwork_halfchains;
        ss_curveNetwork_halfchains << endl;
        
        for (auto& c : curvenetwork.halfchains) {
            ss_curveNetwork_halfchains
                << c.id << " "
                << c.halfedge_front.id << " "
                << c.halfedge_back .id << " "
                << c.patch         .id << " "
                << c.edgechain     .id << endl;
        }
        
        xml_curveNetwork_halfchains->InsertEndChild(xml_doc.NewText(ss_curveNetwork_halfchains.str().c_str()));
        xml_curveNetwork_halfchains->InsertEndChild(xml_doc.NewComment("id halfedge_front halfedge_back patch edgechain"));
    }
    
    {   // edgecains
        auto xml_curveNetwork_edgechains = xml_curveNetwork->InsertEndChild(xml_doc.NewElement("edgechains"))->ToElement();
        xml_curveNetwork_edgechains->SetAttribute("num", static_cast<int>(curvenetwork.edgechains.size()));
        ostringstream ss_curveNetwork_edgechains;
        ss_curveNetwork_edgechains << endl;
        
        for (auto& e : curvenetwork.edgechains) {
            ss_curveNetwork_edgechains
                << e.id << " "
                << e.halfchain[0].id << " "
                << e.halfchain[1].id << " "
                << e.num_subdiv << endl;
        }
        
        xml_curveNetwork_edgechains->InsertEndChild(xml_doc.NewText(ss_curveNetwork_edgechains.str().c_str()));
        xml_curveNetwork_edgechains->InsertEndChild(xml_doc.NewComment("id halfchain_0 halfchain_1 num_subdiv"));
    }
    
    {   // patches
        auto xml_curveNetwork_patches = xml_curveNetwork->InsertEndChild(xml_doc.NewElement("patches"))->ToElement();
        xml_curveNetwork_patches->SetAttribute("num", static_cast<int>(curvenetwork.patches.size()));
        
        for (auto& patch : curvenetwork.patches) {
            auto xml_curveNetwork_patch = xml_curveNetwork_patches->InsertEndChild(xml_doc.NewElement("patch"))->ToElement();
            xml_curveNetwork_patch->SetAttribute("id"        , patch.id          );
            xml_curveNetwork_patch->SetAttribute("halfchain" , patch.halfchain.id);
            xml_curveNetwork_patch->SetAttribute("is_failure", patch.is_failure  );
            
            {   // patch parameter
                auto& param = patch.param;
                int num_sides = param.get_num_sides();
                auto xml_patch_param = xml_curveNetwork_patch->InsertEndChild(xml_doc.NewElement("param"))->ToElement();
                xml_patch_param->SetAttribute("num_sides", num_sides);
                xml_patch_param->SetAttribute("pattern_id", param.pattern_id);
                // permutation
                xml_patch_param->SetAttribute("permutation_id", param.permutation.id);
                // l, variables
                auto xml_param_l        = xml_patch_param->InsertEndChild(xml_doc.NewElement("l"          ))->ToElement();
                auto xml_param_variable = xml_patch_param->InsertEndChild(xml_doc.NewElement("variable"   ))->ToElement();
                ostringstream ss_param_l       ;
                ostringstream ss_param_variable;
                for (int i = 0; i < num_sides; ++i)
                    ss_param_l << param.l[i] << " ";
                int num_variables = patchgen::get_num_variables(num_sides, param.pattern_id);
                for (int i = 0; i < num_variables; ++i)
                    ss_param_variable << patchgen::get_variable(num_sides, param.pattern_id, param, i) << " ";
                xml_param_l       ->InsertEndChild(xml_doc.NewText(ss_param_l          .str().c_str()));
                xml_param_variable->InsertEndChild(xml_doc.NewText(ss_param_variable   .str().c_str()));
            }
            
            {   // patch vertices
                auto xml_patch_vertices = xml_curveNetwork_patch->InsertEndChild(xml_doc.NewElement("patch_vertices"))->ToElement();
                xml_patch_vertices->SetAttribute("num", static_cast<int>(patch.n_vertices()));
                ostringstream ss_patch_vertices;
                ss_patch_vertices << endl;
                
                for (auto v : patch.vertices()) {
                    auto& vdata = patch.data(v);
                    
                    ss_patch_vertices
                        << vdata.pn[0] << " " << vdata.pn[1] << " " << vdata.pn[2] << " "
                        << vdata.pn[3] << " " << vdata.pn[4] << " " << vdata.pn[5] << " "
                        << vdata.patchgen.corner_index << " "
                        << static_cast<int>(vdata.patchgen.tag) << endl;
                }
                xml_patch_vertices->InsertEndChild(xml_doc.NewText(ss_patch_vertices.str().c_str()));
                xml_patch_vertices->InsertEndChild(xml_doc.NewComment("px py pz nx ny nz corner_index tag"));
            }
            
            {   // patch faces
                auto xml_patch_faces = xml_curveNetwork_patch->InsertEndChild(xml_doc.NewElement("patch_faces"))->ToElement();
                xml_patch_faces->SetAttribute("num", static_cast<int>(patch.n_faces()));
                ostringstream ss_patch_faces;
                ss_patch_faces << endl;
                
                for (auto f : patch.faces()) {
                    for (auto v : patch.fv_range(f))
                        ss_patch_faces << v.idx() << " ";
                    ss_patch_faces << endl;
                }
                xml_patch_faces->InsertEndChild(xml_doc.NewText(ss_patch_faces.str().c_str()));
            }
        }
        
    }
    
    {   // material texture (PNG encoded using base64)
        auto xml_texture = xml_root->InsertEndChild(xml_doc.NewElement("texture"))->ToElement();
        xml_texture->InsertEndChild(xml_doc.NewText((material_texture_png_base64).c_str()));
    }
    
    // save (file extension automatically added)
    string fname_ext = fname.substr(fname.size() - 4, 4);
    if (fname_ext != ".xml" && fname_ext != ".XML")
        fname += ".xml";
    if (xml_doc.SaveFile(fname.c_str())) return false;
    
    configTemp.autoSave.filename = fname;
    configTemp.autoSave.unsaved  = false;
    
    return true;
}

int SketchRetopo::get_version() const { return 1; }

bool SketchRetopo::xml_load_batch(const string& fname) {
    string fname_core;
    int fname_num;
    if (!decompose_xml_fname(fname, fname_core, fname_num)) return false;
    
    memento.init(fname_num + 1);
    
    for (int i = 0; i <= fname_num; ++i) {
        ostringstream ss;
        ss << fname_core << '.' << i << ".xml";
        auto fname_i = ss.str();
        if (xml_load(fname_i))
            memento.store(curvenetwork);
    }
    
    return true;
}
