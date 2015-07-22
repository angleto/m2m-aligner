project("m2m-aligner")

cmake_minimum_required(VERSION 2.8)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(PkgConfig REQUIRED)

SET( LINKING_MODE SHARED)
#SET( LINKING_MODE STATIC)

MESSAGE("Linking mode: ${LINKING_MODE}")

add_definitions(
	-std=c++11
)

SET( CMAKE_CXX_FLAGS_RELEASE " -fPIC -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time" CACHE STRING "" FORCE)
SET( CMAKE_CXX_FLAGS_DEBUG "${GOP_WARNING_FLAGS} -fPIC -g3 -O0" CACHE STRING "" FORCE )

FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(Boost COMPONENTS regex REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
FIND_PACKAGE(Boost COMPONENTS filesystem REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS gil REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)
FIND_PACKAGE(Boost COMPONENTS system REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS iostreams REQUIRED)

include_directories(
        "${Boost_INCLUDE_DIRS}"
)

link_directories(
        "/opt/local/lib"
        "/usr/local/lib"
)

# defining m2m library name and include path
SET(M2M_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET(M2M_LIBRARIES m2mAligner)
INCLUDE_DIRECTORIES(${M2M_INCLUDE_PATH})

add_subdirectory(lib)
add_subdirectory(apps)
