add_executable(user usermain.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/)
link_directories(${HPX_LIBRARY_DIRS})

target_link_libraries(user PUBLIC dashmm ${HPX_LDFLAGS})