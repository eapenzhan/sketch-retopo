PROJECT(SketchRetopo)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES(./src)

IF(UNIX)
  SET(shared_if_unix SHARED)
ELSE()
  SET(shared_if_unix "")
ENDIF()


# Windows compiler settings
IF(WIN32)
  ADD_DEFINITIONS(
    -D__STDC_LIMIT_MACROS
    -DEIGEN_DONT_ALIGN_STATICALLY
    -DGLEW_STATIC
    -DFREEGLUT_STATIC)
  SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT ")
  SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   /MTd")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT ")
ENDIF()

# Mac/Linux compiler settings
IF(UNIX)
#  SET(CMAKE_BUILD_TYPE Debug)
  INCLUDE_DIRECTORIES(
    /usr/local/include
    /opt/local/include
  )
  LINK_DIRECTORIES(
    /usr/local/lib
    /opt/local/lib
  )
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -msse4 -stdlib=libc++ -ObjC++ -framework OpenGL -framework GLUT -framework Cocoa")
    ADD_DEFINITIONS(-D_MACOSX)
  ELSE()
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -msse4")
    ADD_DEFINITIONS(-D_UNIX)
  ENDIF()
ENDIF()

# AntTweakBar
ADD_DEFINITIONS(
  -DTW_NO_LIB_PRAGMA
  -DTW_STATIC)
INCLUDE_DIRECTORIES(./src/AntTweakBar)
FILE(GLOB_RECURSE anttweakbar_src
  ./src/AntTweakBar/*.c
  ./src/AntTweakBar/*.h
  ./src/AntTweakBar/*.cpp
)
ADD_LIBRARY(AntTweakBar ${anttweakbar_src})

# Embree-1.1
INCLUDE_DIRECTORIES(
  ./src/embree-1.1
  ./src/embree-1.1/common
  ./src/embree-1.1/rtcore)
FILE(GLOB embree_common_math_src         ./src/embree-1.1/common/math/*.h)
FILE(GLOB embree_common_symd_src         ./src/embree-1.1/common/symd/*.h)
FILE(GLOB embree_common_sys_src          ./src/embree-1.1/common/sys/*.h ./src/embree-1.1/common/sys/*.cpp)
FILE(GLOB embree_common_sys_stl_src      ./src/embree-1.1/common/sys/stl/*.h)
FILE(GLOB embree_common_sys_stl_math_src ./src/embree-1.1/common/sys/stl/math/*.h)
FILE(GLOB embree_common_sys_sync_src     ./src/embree-1.1/common/sys/sync/*.h ./src/embree-1.1/common/sys/sync/*.cpp)
FILE(GLOB embree_rtcore_bvh2_src         ./src/embree-1.1/rtcore/bvh2/*.h   ./src/embree-1.1/rtcore/bvh2/*.cpp  )
FILE(GLOB embree_rtcore_bvh4_src         ./src/embree-1.1/rtcore/bvh4/*.h   ./src/embree-1.1/rtcore/bvh4/*.cpp  )
FILE(GLOB embree_rtcore_bvh4mb_src       ./src/embree-1.1/rtcore/bvh4mb/*.h ./src/embree-1.1/rtcore/bvh4mb/*.cpp)
FILE(GLOB embree_rtcore_common_src       ./src/embree-1.1/rtcore/common/*.h   ./src/embree-1.1/rtcore/common/*.cpp  )
FILE(GLOB embree_rtcore_triangle_src     ./src/embree-1.1/rtcore/triangle/*.h)
SOURCE_GROUP("embree-1.1\\common\\math" FILES ${embree_common_math_src})
SOURCE_GROUP("embree-1.1\\common\\simd" FILES ${embree_common_simd_src})
SOURCE_GROUP("embree-1.1\\common\\sys" FILES ${embree_common_sys_src})
SOURCE_GROUP("embree-1.1\\common\\sys\\stl" FILES ${embree_common_sys_stl_src})
SOURCE_GROUP("embree-1.1\\common\\sys\\stl\\math" FILES ${embree_common_sys_stl_math_src})
SOURCE_GROUP("embree-1.1\\common\\sys\\sync" FILES ${embree_common_sys_sync_src})
SOURCE_GROUP("embree-1.1\\rtcore\\bvh2" FILES ${embree_rtcore_bvh2_src})
SOURCE_GROUP("embree-1.1\\rtcore\\bvh4" FILES ${embree_rtcore_bvh4_src})
SOURCE_GROUP("embree-1.1\\rtcore\\bvh4mb" FILES ${embree_rtcore_bvh4mb_src})
SOURCE_GROUP("embree-1.1\\rtcore\\common" FILES ${embree_rtcore_common_src})
SOURCE_GROUP("embree-1.1\\rtcore\\triangle" FILES ${embree_rtcore_triangle_src})
ADD_LIBRARY(embree-1.1
  ${embree_common_math_src}
  ${embree_common_symd_src}
  ${embree_common_sys_src}
  ${embree_common_sys_stl_src}
  ${embree_common_sys_stl_math_src}
  ${embree_common_sys_sync_src}
  ${embree_rtcore_bvh2_src}
  ${embree_rtcore_bvh4_src}
  ${embree_rtcore_bvh4mb_src}
  ${embree_rtcore_common_src}
  ${embree_rtcore_triangle_src})

# LodePNG
INCLUDE_DIRECTORIES(./src/lodepng)
SET(lodepng_src ./src/lodepng/lodepng.h ./src/lodepng/lodepng.cpp)
ADD_LIBRARY(lodepng ${lodepng_src})

# lp_solve_5.5
ADD_DEFINITIONS(
  -DLoadInverseLib=FALSE
  -DLoadLanguageLib=FALSE
  -DLoadableBlasLib=0
  -DINVERSE_ACTIVE=INVERSE_LUSOL
  -DRoleIsExternalInvEngine
  -DNOFORTIFY
  -DCHECK_SOLUTION
  -DYY_NEVER_INTERACTIVE
  -DPARSER_LP
  -DLP_MAXLINELEN=0)
INCLUDE_DIRECTORIES(
  ./src/lp_solve_5.5
  ./src/lp_solve_5.5/bfp
  ./src/lp_solve_5.5/bfp/bfp_LUSOL
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL
  ./src/lp_solve_5.5/colamd
  ./src/lp_solve_5.5/shared)
FILE(GLOB lp_solve_src ./src/lp_solve_5.5/*.h ./src/lp_solve_5.5/*.c)
SET(lp_solve_bfp_src
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c)
SET(lp_solve_colamd_src
  ./src/lp_solve_5.5/colamd/colamd.h
  ./src/lp_solve_5.5/colamd/colamd.c)
FILE(GLOB lp_solve_shared_src ./src/lp_solve_5.5/shared/*.h ./src/lp_solve_5.5/shared/*.c)
SOURCE_GROUP("lp_solve_5.5" FILES ${lp_solve_src})
SOURCE_GROUP("lp_solve_5.5\\bfp" FILES ${lp_solve_bfp_src})
SOURCE_GROUP("lp_solve_5.5\\colamd" FILES ${lp_solve_colamd_src})
SOURCE_GROUP("lp_solve_5.5\\shared" FILES ${lp_solve_shared_src})
ADD_LIBRARY(lp_solve
  ${lp_solve_src}
  ${lp_solve_bfp_src}
  ${lp_solve_colamd_src}
  ${lp_solve_shared_src})

# OpenMesh
ADD_DEFINITIONS(-DOM_STATIC_BUILD)
IF(WIN32)
  ADD_DEFINITIONS(
    -D_USE_MATH_DEFINES
    -DNOMINMAX)
ENDIF(WIN32)
FILE(GLOB openmesh_core_geometry_src ./src/OpenMesh/Core/Geometry/*.hh ./src/OpenMesh/Core/Geometry/*.cc)
FILE(GLOB openmesh_core_io_src ./src/OpenMesh/Core/IO/*.hh ./src/OpenMesh/Core/IO/*.cc)
FILE(GLOB openmesh_core_io_exporter_src ./src/OpenMesh/Core/IO/exporter/*.hh)
FILE(GLOB openmesh_core_io_importer_src ./src/OpenMesh/Core/IO/importer/*.hh)
FILE(GLOB openmesh_core_io_reader_src ./src/OpenMesh/Core/IO/reader/*.hh ./src/OpenMesh/Core/IO/reader/*.cc)
FILE(GLOB openmesh_core_io_writer_src ./src/OpenMesh/Core/IO/writer/*.hh ./src/OpenMesh/Core/IO/writer/*.cc)
FILE(GLOB openmesh_core_mesh_src ./src/OpenMesh/Core/Mesh/*.hh ./src/OpenMesh/Core/Mesh/*.cc)
FILE(GLOB openmesh_core_system_src ./src/OpenMesh/Core/System/*.hh ./src/OpenMesh/Core/System/*.cc)
FILE(GLOB openmesh_core_utils_src ./src/OpenMesh/Core/Utils/*.hh ./src/OpenMesh/Core/Utils/*.cc)
SOURCE_GROUP("OpenMesh\\Core\\Geometry" FILES ${openmesh_core_geometry_src})
SOURCE_GROUP("OpenMesh\\Core\\IO" FILES ${openmesh_core_io_src})
SOURCE_GROUP("OpenMesh\\Core\\IO\\exporter" FILES ${openmesh_core_io_exporter_src})
SOURCE_GROUP("OpenMesh\\Core\\IO\\importer" FILES ${openmesh_core_io_importer_src})
SOURCE_GROUP("OpenMesh\\Core\\IO\\reader" FILES ${openmesh_core_io_reader_src})
SOURCE_GROUP("OpenMesh\\Core\\IO\\writer" FILES ${openmesh_core_io_writer_src})
SOURCE_GROUP("OpenMesh\\Core\\Mesh" FILES ${openmesh_core_mesh_src})
SOURCE_GROUP("OpenMesh\\Core\\System" FILES ${openmesh_core_system_src})
SOURCE_GROUP("OpenMesh\\Core\\Utils" FILES ${openmesh_core_utils_src})
ADD_LIBRARY(OpenMesh
  ${openmesh_core_geometry_src}
  ${openmesh_core_io_src}
  ${openmesh_core_io_exporter_src}
  ${openmesh_core_io_importer_src}
  ${openmesh_core_io_reader_src}
  ${openmesh_core_io_writer_src}
  ${openmesh_core_mesh_src}
  ${openmesh_core_system_src}
  ${openmesh_core_utils_src})

# TinyXml2
INCLUDE_DIRECTORIES(./src/tinyxml2)
SET(tinyxml2_src ./src/tinyxml2/tinyxml2.h ./src/tinyxml2/tinyxml2.cpp)
ADD_LIBRARY(tinyxml2 ${tinyxml2_src})

# patchgen
FILE(GLOB patchgen_src ./src/patchgen/*.hh ./src/patchgen/*.cc)
ADD_LIBRARY(patchgen ${patchgen_src})

# SketchRetopo
FILE(GLOB sketchretopo_src ./src/SketchRetopo/*.hh ./src/SketchRetopo/*.cc)
FILE(GLOB sketchretopo_curvenetwork_src ./src/SketchRetopo/curvenetwork/*.hh ./src/SketchRetopo/curvenetwork/*.cc)
FILE(GLOB sketchretopo_state_src ./src/SketchRetopo/state/*.hh ./src/SketchRetopo/state/*.cc)
SOURCE_GROUP("SketchRetopo" FILES ${sketchretopo_src})
SOURCE_GROUP("SketchRetopo\\curvenetwork" FILES ${sketchretopo_curvenetwork_src})
SOURCE_GROUP("SketchRetopo\\state" FILES ${sketchretopo_state_src})
ADD_EXECUTABLE(SketchRetopo
  ${sketchretopo_src}
  ${sketchretopo_curvenetwork_src}
  ${sketchretopo_state_src})

# SketchRetpo linker settings
IF(WIN32)
  TARGET_LINK_LIBRARIES(SketchRetopo
    AntTweakBar
    embree-1.1
    lodepng
    patchgen
    lp_solve
    OpenMesh
    tinyxml2
    OpenGL32
    GLU32
    freeglut_static
    glew32s
  )
ELSEIF(APPLE)
    FIND_LIBRARY(Cocoa_LIB Cocoa)
    FIND_LIBRARY(GLUT_LIB GLUT)
    FIND_LIBRARY(GLEW_LIB GLEW)
    FIND_LIBRARY(OpenGL_LIB OpenGL)
    MARK_AS_ADVANCED(Cocoa_LIB GLUT_LIB GLEW_LIB OpenGL_LIB)
  TARGET_LINK_LIBRARIES(SketchRetopo
    AntTweakBar
    embree-1.1
    lodepng
    patchgen
    lp_solve
    OpenMesh
    tinyxml2
    ${Cocoa_LIB}
    ${GLUT_LIB}
    ${GLEW_LIB}
    ${OpenGL_LIB}
  )
ELSEIF(UNIX)
  TARGET_LINK_LIBRARIES(SketchRetopo
    AntTweakBar
    embree-1.1
    lodepng
    patchgen
    lp_solve
    OpenMesh
    tinyxml2
    GL
    GLU
    GLEW
    glut
    X11
    pthread
  )
ENDIF()

# patchgen_demo
FILE(GLOB patchgen_demo_src ./src/patchgen_demo/*.hh ./src/patchgen_demo/*.cc)
SOURCE_GROUP("patchgen_demo" FILES ${patchgen_demo_src})
ADD_EXECUTABLE(patchgen_demo
  ${patchgen_demo_src})
IF(WIN32)
  TARGET_LINK_LIBRARIES(patchgen_demo
    AntTweakBar
    OpenMesh
    patchgen
    lp_solve
    OpenGL32
    GLU32
    freeglut_static
    glew32
  )
ELSEIF(APPLE)
    FIND_LIBRARY(Cocoa_LIB Cocoa)
    FIND_LIBRARY(GLUT_LIB GLUT)
    FIND_LIBRARY(GLEW_LIB GLEW)
    FIND_LIBRARY(OpenGL_LIB OpenGL)
    MARK_AS_ADVANCED(Cocoa_LIB GLUT_LIB GLEW_LIB OpenGL_LIB)
  TARGET_LINK_LIBRARIES(patchgen_demo
    AntTweakBar
    OpenMesh
    patchgen
    lp_solve
    ${Cocoa_LIB}
    ${GLUT_LIB}
    ${GLEW_LIB}
    ${OpenGL_LIB}
  )
ELSEIF(UNIX)
  TARGET_LINK_LIBRARIES(patchgen_demo
    AntTweakBar
    patchgen
    lp_solve
    OpenMesh
    GL
    GLU
    GLEW
    glut
    X11
    pthread
  )
ENDIF()
