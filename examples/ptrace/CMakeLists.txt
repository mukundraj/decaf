cmake_minimum_required (VERSION 2.6)
project (ptrace)

set(default_build_type "Debug")
# SET(VISUALS FALSE CACHE BOOL "Enable/Disbale visualizations")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${DIY_INCLUDE_DIRS})


find_package (NetCDF REQUIRED)
include_directories (${NETCDF_INCLUDE_DIR})

# Find PnetCDF
find_path                   (PNETCDF_DIR lib/libpnetcdf.a)
set                         (PNETCDF_INCLUDE_DIR ${PNETCDF_DIR}/include)
set                         (PNETCDF_LIBRARY ${PNETCDF_DIR}/lib/libpnetcdf.a)
include_directories (${PNETCDF_INCLUDE_DIR})
link_directories(${PNETCDF_DIR}/lib)


#--------------------
# DEPENDENCY: nabo
#--------------------
find_package(libnabo REQUIRED PATHS ${LIBNABO_INSTALL_DIR})
#include(libnaboConfig)
include_directories(${libnabo_INCLUDE_DIRS})
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${libnabo_LIBRARIES})
message(STATUS "libnabo found, version ${libnabo_VERSION} (include=${libnabo_INCLUDE_DIRS} libs=${libnabo_LIBRARIES})")

set (EIGEN3_INCLUDE_DIR "/home/mraj/libraries/eigen-eigen-323c052e1731")
set (EIGEN3_INCLUDE_DIR "/homes/mraj/work/libraries/eigen-eigen-323c052e1731")
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}"  )



# add_executable 	(ptrace src/main.cpp src/core/ptrace.cpp src/io/mpas_io.cpp src/core/pathline.cpp src/core/ptrace.cpp src/core/block.cpp src/utils/misc.cpp src/core/geometry_utils.cpp src/core/advect.cpp)
# target_include_directories(ptrace PRIVATE src/io src/core src/utils)
# target_link_libraries (ptrace  ${NetCDF_LIBRARIES} ${bredala_transports} ${PNETCDF_LIBRARY} bredala_datamodel manala decaf ${libnabo_LIBRARIES} ${PNETCDF_LIBRARY} )

add_executable 	(sptrace src/pmain.cpp src/core/ptrace.cpp src/io/mpas_io.cpp src/core/pathline.cpp src/core/ptrace.cpp src/core/block.cpp src/utils/misc.cpp src/core/geometry_utils.cpp src/core/advect.cpp)
target_include_directories(sptrace PRIVATE src/io src/core src/utils)
target_link_libraries (sptrace  ${NetCDF_LIBRARIES} ${bredala_transports} ${PNETCDF_LIBRARY} bredala_datamodel manala decaf ${libnabo_LIBRARIES} ${PNETCDF_LIBRARY} )

add_library     (adapter MODULE src/adapter/adapter.cpp)
# .SO file extension on Linux/Mac OS
set_target_properties       (adapter          PROPERTIES SUFFIX ".so")

# # Don't add a 'lib' prefix to the shared library
# set_target_properties       (mod_adapter       PROPERTIES PREFIX "")