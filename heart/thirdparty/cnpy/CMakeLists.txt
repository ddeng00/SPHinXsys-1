CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

project(CNPY)

#option(ENABLE_STATIC "Build static (.a) library" ON)

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
FIND_PACKAGE(LibZip REQUIRED)

add_library(cnpy STATIC "cnpy.cpp")
target_include_directories(cnpy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(cnpy PUBLIC ${LIBZIP_INCLUDE_DIR_ZIP})
target_include_directories(cnpy PUBLIC "extern/libzippp/")
target_link_libraries(cnpy PUBLIC ${LIBZIP_LIBRARIES})
target_compile_features(cnpy PRIVATE cxx_std_14)

#install(TARGETS "cnpy" LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#
#if(ENABLE_STATIC)
#    add_library(cnpy-static STATIC "cnpy.cpp")
#    set_target_properties(cnpy-static PROPERTIES OUTPUT_NAME "cnpy")
#    install(TARGETS "cnpy-static" ARCHIVE DESTINATION lib)
#endif(ENABLE_STATIC)
#
#install(FILES "cnpy.h" DESTINATION include)
#install(FILES "mat2npz" "npy2mat" "npz2mat" DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#
#add_executable(example1 example1.cpp)
#target_link_libraries(example1 cnpy)
