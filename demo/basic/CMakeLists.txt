add_executable(basic EXCLUDE_FROM_ALL testmain.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/)
link_directories(${HPX_LIBRARY_DIRS})

target_link_libraries(basic PUBLIC dashmm ${HPX_LDFLAGS})