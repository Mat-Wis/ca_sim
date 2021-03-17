cmake_minimum_required( VERSION 3.10 )
set( CMAKE_CXX_STANDARD 20 )
project( ca_sim )

set( CMAKE_CXX_FLAGS "-Wall" )
#-ffast-math

#include_directories( /usr/local/include inc )
include_directories( inc )
link_directories( /usr/local/lib )

add_executable( ca_sim src/main.cpp src/logger.cpp )
target_link_libraries( ca_sim matio )
