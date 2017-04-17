add_executable(collect EXCLUDE_FROM_ALL collect.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
  ${PROJECT_SOURCE_DIR}/test/common/)
link_directories(${HPX_LIBRARY_DIRS})

target_link_libraries(collect PUBLIC dashmm ${HPX_LDFLAGS})