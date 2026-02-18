cmake_minimum_required(VERSION 3.20)
project(liftasm LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LIFTASM_ENABLE_NATIVE "Enable -march=native" ON)
option(LIFTASM_ENABLE_O3     "Enable -O3"           ON)

if (LIFTASM_ENABLE_O3)
  add_compile_options($<$<COMPILE_LANGUAGE:C>:-O3>  $<$<COMPILE_LANGUAGE:CXX>:-O3>)
endif()

if (LIFTASM_ENABLE_NATIVE)
  add_compile_options($<$<COMPILE_LANGUAGE:C>:-march=native>  $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# liftasm sources
file(GLOB LIFTASM_SRC CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(liftasm
  "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
  ${LIFTASM_SRC}
)

target_include_directories(liftasm PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(liftasm PRIVATE OpenMP::OpenMP_CXX)

# Submodules
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

if (NOT EXISTS "${EXTERNAL_DIR}/minimap2/Makefile")
  message(FATAL_ERROR "external/minimap2 missing. Run: git submodule update --init --recursive")
endif()
if (NOT EXISTS "${EXTERNAL_DIR}/WFA2-lib/CMakeLists.txt")
  message(FATAL_ERROR "external/WFA2-lib missing. Run: git submodule update --init --recursive")
endif()
if (NOT EXISTS "${EXTERNAL_DIR}/htslib/bgzf.c")
  message(FATAL_ERROR "external/htslib missing. Run: git submodule update --init --recursive")
endif()

include(ExternalProject)
set(EXTERNAL_BUILD_DIR "${CMAKE_BINARY_DIR}/external_build")

# 1) minimap2 -> libminimap2.a
set(MINIMAP2_SRC "${EXTERNAL_DIR}/minimap2")
set(MINIMAP2_LIB "${MINIMAP2_SRC}/libminimap2.a")

ExternalProject_Add(ext_minimap2
  SOURCE_DIR        "${MINIMAP2_SRC}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM} -C "${MINIMAP2_SRC}" libminimap2.a
  INSTALL_COMMAND   ""
  BUILD_BYPRODUCTS  "${MINIMAP2_LIB}"
)

add_library(minimap2_static STATIC IMPORTED GLOBAL)
set_target_properties(minimap2_static PROPERTIES
  IMPORTED_LOCATION "${MINIMAP2_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${MINIMAP2_SRC}"
)
add_dependencies(minimap2_static ext_minimap2)

# 2) WFA2-lib -> libwfa2.a + libwfa2cpp.a
set(WFA_SRC   "${EXTERNAL_DIR}/WFA2-lib")
set(WFA_BUILD "${EXTERNAL_BUILD_DIR}/WFA2-lib")
set(WFA2_LIB    "${WFA_BUILD}/libwfa2.a")
set(WFA2CPP_LIB "${WFA_BUILD}/libwfa2cpp.a")

ExternalProject_Add(ext_wfa2
  SOURCE_DIR "${WFA_SRC}"
  BINARY_DIR "${WFA_BUILD}"
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  BUILD_COMMAND     ${CMAKE_COMMAND} --build "${WFA_BUILD}" --config Release
  INSTALL_COMMAND   ""
  BUILD_BYPRODUCTS  "${WFA2CPP_LIB}" "${WFA2_LIB}"
)

add_library(wfa2cpp_static STATIC IMPORTED GLOBAL)
set_target_properties(wfa2cpp_static PROPERTIES
  IMPORTED_LOCATION "${WFA2CPP_LIB}"
)
add_dependencies(wfa2cpp_static ext_wfa2)

add_library(wfa2_static STATIC IMPORTED GLOBAL)
set_target_properties(wfa2_static PROPERTIES
  IMPORTED_LOCATION "${WFA2_LIB}"
)
add_dependencies(wfa2_static ext_wfa2)

# 3) htslib -> build static lib
set(HTS_SRC "${EXTERNAL_DIR}/htslib")
set(HTS_LIB "${HTS_SRC}/libhts.a")

ExternalProject_Add(ext_htslib
  SOURCE_DIR "${HTS_SRC}"
  BINARY_DIR "${EXTERNAL_BUILD_DIR}/htslib"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND
    /bin/sh -c
    "set -e; cd '${HTS_SRC}'; \
     ${CMAKE_MAKE_PROGRAM} libhts.a || ${CMAKE_MAKE_PROGRAM}"
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS "${HTS_LIB}"
)

add_library(hts_static STATIC IMPORTED GLOBAL)
set_target_properties(hts_static PROPERTIES
  IMPORTED_LOCATION "${HTS_LIB}"
  # htslib headers are under ${HTS_SRC}/htslib/*.h
  INTERFACE_INCLUDE_DIRECTORIES "${HTS_SRC}"
)
add_dependencies(hts_static ext_htslib)

# Link
target_link_libraries(liftasm PRIVATE
  minimap2_static
  wfa2cpp_static
  wfa2_static
  hts_static
  z bz2 lzma curl crypto ssl pthread dl
)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #include <filesystem>
  int main(){ std::filesystem::path p{\".\"}; return 0; }
" HAVE_STD_FILESYSTEM)
if (NOT HAVE_STD_FILESYSTEM)
  target_link_libraries(liftasm PRIVATE stdc++fs)
endif()

add_dependencies(liftasm ext_minimap2 ext_wfa2 ext_htslib)