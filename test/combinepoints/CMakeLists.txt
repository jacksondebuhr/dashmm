add_executable(combinepoints EXCLUDE_FROM_ALL combine.cc combiner.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
  ${PROJECT_SOURCE_DIR}/test/common/)
