FILE(GLOB sources *.cc)
include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include)
link_directories(${HPX_LIBRARY_DIRS})

add_library(dashmm STATIC ${sources})
target_link_libraries(dashmm ${HPX_LDFLAGS})

install (TARGETS dashmm DESTINATION lib)
