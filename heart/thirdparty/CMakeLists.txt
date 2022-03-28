set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library(CLI11 INTERFACE CLI11/CLI11.hpp)
target_include_directories(CLI11 INTERFACE CLI11)

add_subdirectory(cnpy)