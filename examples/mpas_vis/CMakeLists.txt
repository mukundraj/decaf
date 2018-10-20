cmake_minimum_required (VERSION 2.6)
project (mpas_vis)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

#add_subdirectory(src)
#add_subdirectory(src/io)
#add_subdirectory(src/filters)

#SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++14 -O2")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}"  )

set (EIGEN3_INCLUDE_DIR "/Users/mukundraj/Desktop/work/libraries/eigen_3.3.4/include/eigen3")
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}"  )

find_package (NetCDF REQUIRED)
include_directories (${NETCDF_INCLUDE_DIR})

#set(VTK_DIR "/homes/mraj/work/projects/VTK-8.1.1/build")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include)

#--------------------
# DEPENDENCY: nabo
#--------------------
find_package(libnabo REQUIRED PATHS ${LIBNABO_INSTALL_DIR})
#include(libnaboConfig)
include_directories(${libnabo_INCLUDE_DIRS})
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${libnabo_LIBRARIES})
message(STATUS "libnabo found, version ${libnabo_VERSION} (include=${libnabo_INCLUDE_DIRS} libs=${libnabo_LIBRARIES})")

#--------------------
# DEPENDENCY: Pnetcdf
#--------------------
find_path(PNETCDF_DIR lib/libpnetcdf.a)
set (PNETCDF_INCLUDE_DIR ${PNETCDF_DIR}/include)
set (PNETCDF_LIBRARY ${PNETCDF_DIR}/lib/libpnetcdf.a)
include_directories (${PNETCDF_INCLUDE_DIR})

add_executable(flow_main src/main.cpp src/utils/interpolators.cpp src/io/flow.cpp src/io/mpaso.cpp
    src/filters/streamlines.cpp src/filters/pathlines.cpp src/utils/advect.cpp src/io/block_io.cpp src/utils/misc.cpp src/utils/pblock.cpp)

add_library     (mod_mpas_adapter2 MODULE src/adapter/mpas_adapter2.cpp)

# .SO file extension on Linux/Mac OS
set_target_properties       (mod_mpas_adapter2          PROPERTIES SUFFIX ".so")

# Don't add a 'lib' prefix to the shared library
set_target_properties       (mod_mpas_adapter2       PROPERTIES PREFIX "")

target_link_libraries (flow_main ${NETCDF_LIBRARIES} ${VTK_LIBRARIES} ${libnabo_LIBRARIES} ${libraries} ${bredala_transports} ${PNETCDF_LIBRARY} bredala_datamodel manala decaf)

install                     (TARGETS flow_main mod_mpas_adapter2 
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/mpas_vis/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)

install                     (FILES generate_script.py
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/mpas_vis/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
