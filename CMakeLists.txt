cmake_minimum_required (VERSION 3.4)

project (dashmm)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

find_package (PkgConfig)
pkg_search_module (HPX REQUIRED hpx)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(demo)
add_subdirectory(test)
