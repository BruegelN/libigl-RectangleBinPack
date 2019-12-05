cmake_minimum_required(VERSION 3.1)
project(example)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)
option(LIBIGL_BUILD_PYTHON           "Build libigl python bindings" OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_BUILD_TESTS      "Build libigl unit test"        OFF)
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"         OFF)
option(LIBIGL_EXPORT_TARGETS   "Export libigl CMake targets"   OFF)
find_package(LIBIGL REQUIRED QUIET)

add_subdirectory(RectangleBinPack)

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw RectangleBinPack)
