cmake_minimum_required (VERSION 2.6)
project (ptrace)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${DIY_INCLUDE_DIRS})

add_executable 	(ptrace src/main.cpp)
target_link_libraries (ptrace ${bredala_transports} ${PNETCDF_LIBRARY} bredala_datamodel manala decaf)

add_library     (adapter MODULE src/adapter/adapter.cpp)
# .SO file extension on Linux/Mac OS
set_target_properties       (adapter          PROPERTIES SUFFIX ".so")

# # Don't add a 'lib' prefix to the shared library
# set_target_properties       (mod_adapter       PROPERTIES PREFIX "")