cmake_minimum_required(VERSION 3.20)
project(liftasm LANGUAGES C CXX)

# GCC >= 12.2.0
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "12.2.0")
    message(FATAL_ERROR
      "GCC >= 12.2.0 is required, but found ${CMAKE_CXX_COMPILER_VERSION} "
      "(${CMAKE_CXX_COMPILER}). Please load a newer GCC module (e.g. GCC/12.2.0)."
    )
  endif()
else()
  message(FATAL_ERROR
    "This project requires GCC >= 12.2.0. Detected compiler: "
    "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})."
  )
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LIFTASM_ENABLE_NATIVE "Enable -march=native" OFF)
option(LIFTASM_ENABLE_O3     "Enable -O3"           ON)
option(LIFTASM_ENABLE_CURL   "Enable HTSlib remote URL support when libcurl is available" ON)

if (LIFTASM_ENABLE_O3)
  add_compile_options($<$<COMPILE_LANGUAGE:C>:-O3>  $<$<COMPILE_LANGUAGE:CXX>:-O3>)
endif()

if (LIFTASM_ENABLE_NATIVE)
  add_compile_options($<$<COMPILE_LANGUAGE:C>:-march=native>  $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(Threads REQUIRED)

set(LIFTASM_HAVE_CURL OFF)
if (LIFTASM_ENABLE_CURL)
  find_package(CURL QUIET)
  if (CURL_FOUND)
    set(LIFTASM_HAVE_CURL ON)
  endif()
endif()

if (LIFTASM_HAVE_CURL)
  set(HTS_CURL_MODE "curl")
  set(HTS_CURL_OBJECTS "hfile_libcurl.o")
  set(HTS_DEFAULT_LIBS "-lz -lm -lbz2 -llzma -lcurl")
  message(STATUS "HTSlib remote URL support: enabled")
else()
  set(HTS_CURL_MODE "no-curl")
  set(HTS_CURL_OBJECTS "")
  set(HTS_DEFAULT_LIBS "-lz -lm -lbz2 -llzma")
  if (LIFTASM_ENABLE_CURL)
    message(STATUS "HTSlib remote URL support: disabled (libcurl not found)")
  else()
    message(STATUS "HTSlib remote URL support: disabled")
  endif()
endif()

# sources
file(GLOB LIFTASM_SRC CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(liftasm
  "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
  ${LIFTASM_SRC}
)

target_include_directories(liftasm PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external/WFA2-lib"
  "${CMAKE_CURRENT_SOURCE_DIR}/external/minimap2"
  "${CMAKE_CURRENT_SOURCE_DIR}/external/htslib"
)

# submodules
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

set(LIFTASM_JOBS_DEFAULT 5)
if(DEFINED ENV{LIFTASM_JOBS})
  set(LIFTASM_JOBS "$ENV{LIFTASM_JOBS}")
else()
  set(LIFTASM_JOBS "${LIFTASM_JOBS_DEFAULT}")
endif()

# 1) minimap2 -> libminimap2.a
set(MINIMAP2_SRC "${EXTERNAL_DIR}/minimap2")
set(MINIMAP2_LIB "${MINIMAP2_SRC}/libminimap2.a")

ExternalProject_Add(ext_minimap2
  SOURCE_DIR        "${MINIMAP2_SRC}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM} -C "${MINIMAP2_SRC}" -j ${LIFTASM_JOBS} libminimap2.a
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
  BUILD_COMMAND     ${CMAKE_COMMAND} --build "${WFA_BUILD}" --config Release -- -j ${LIFTASM_JOBS}
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

# 3) htslib -> custom build libhts.a
set(HTS_SRC "${EXTERNAL_DIR}/htslib")
set(HTS_BUILD "${EXTERNAL_BUILD_DIR}/htslib-${HTS_CURL_MODE}")
set(HTS_LIB "${HTS_BUILD}/libhts.a")

set(HTSCODECS_CONFIGURE_AC "${HTS_SRC}/htscodecs/configure.ac")
file(STRINGS "${HTSCODECS_CONFIGURE_AC}" HTSCODECS_AC_INIT
  LIMIT_COUNT 1
  REGEX "^AC_INIT\\(htscodecs,"
)
string(REGEX MATCH "[0-9]+(\\.[0-9]+)+" HTSCODECS_VERSION "${HTSCODECS_AC_INIT}")
if (NOT HTSCODECS_VERSION)
  message(FATAL_ERROR "Failed to read the htscodecs version from ${HTSCODECS_CONFIGURE_AC}")
endif()
set(HTSCODECS_VERSION_HEADER "${CMAKE_BINARY_DIR}/generated/htscodecs_version.h")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")
file(CONFIGURE
  OUTPUT "${HTSCODECS_VERSION_HEADER}"
  CONTENT "#define HTSCODECS_VERSION_TEXT \"${HTSCODECS_VERSION}\"\n"
  @ONLY
)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${HTSCODECS_CONFIGURE_AC}")

set(HTS_CURL_CONFIG_COMMAND)
if (NOT LIFTASM_HAVE_CURL)
  list(APPEND HTS_CURL_CONFIG_COMMAND
    COMMAND /bin/sh -c "sed -i '/^#define HAVE_LIBCURL 1$/d' '${HTS_BUILD}/config.h'"
  )
endif()

add_custom_command(
  OUTPUT "${HTS_LIB}"
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTS_BUILD}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${HTS_SRC}" "${HTS_BUILD}"
  COMMAND ${CMAKE_COMMAND} -E rm -f "${HTS_BUILD}/.git" "${HTS_BUILD}/htscodecs/.git"
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${HTS_BUILD}" distclean
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${HTSCODECS_VERSION_HEADER}" "${HTS_BUILD}/htscodecs/htscodecs/version.h"
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${HTS_BUILD}" config.h
  ${HTS_CURL_CONFIG_COMMAND}
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${HTS_BUILD}" -j ${LIFTASM_JOBS}
    "NONCONFIGURE_OBJS=${HTS_CURL_OBJECTS}"
    "htslib_default_libs=${HTS_DEFAULT_LIBS}"
    libhts.a
  DEPENDS "${HTSCODECS_VERSION_HEADER}"
  COMMENT "Building htslib (libhts.a, ${HTS_CURL_MODE})"
  VERBATIM
)

add_custom_target(htslib_build ALL DEPENDS "${HTS_LIB}")

add_library(hts_static STATIC IMPORTED GLOBAL)
set_target_properties(hts_static PROPERTIES
  IMPORTED_LOCATION "${HTS_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${HTS_SRC}"
)
add_dependencies(hts_static htslib_build)

# link
target_link_libraries(liftasm PRIVATE
  minimap2_static
  wfa2cpp_static
  wfa2_static
  hts_static
  ZLIB::ZLIB
  BZip2::BZip2
  LibLZMA::LibLZMA
  Threads::Threads
)

if (LIFTASM_HAVE_CURL)
  target_link_libraries(liftasm PRIVATE CURL::libcurl)
endif()

add_dependencies(liftasm ext_minimap2 ext_wfa2 htslib_build)
