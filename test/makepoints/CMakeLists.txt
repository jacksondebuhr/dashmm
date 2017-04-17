add_executable(makepoints EXCLUDE_FROM_ALL makepoints.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
  ${PROJECT_SOURCE_DIR}/test/common/)
link_directories(${HPX_LIBRARY_DIRS})

target_link_libraries(makepoints PUBLIC dashmm ${HPX_LDFLAGS})