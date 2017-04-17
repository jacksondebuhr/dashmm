cmake_minimum_required (VERSION 3.6)

project (dashmm)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package (PkgConfig)
pkg_search_module (HPX REQUIRED hpx)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(demo)
add_subdirectory(test)

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
set (CPACK_RESOURCE_FILE README ${CMAKE_CURRENT_SOURCE_DIR}/README)
set (CPACK_RESOURCE_FILE INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/INSTALL)
set (CPACK_RESOURCE_FILE AUTHORS ${CMAKE_CURRENT_SOURCE_DIR}/AUTHORS)
set (CPACK_RESOURCE_FILE doc/DASHMMUserGuide.pdf ${CMAKE_CURRENT_SOURCE_DIR}/doc/DASHMMUserGuide.pdf)
set (CPACK_PACKAGE_VERSION_MAJOR 2)
set (CPACK_PACKAGE_VERSION_MINOR 0)
set (CPACK_PACKAGE_VERSION_PATCH 0)
include (CPack)