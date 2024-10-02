cmake_minimum_required(VERSION 3.18)
set(proj_name hellopybind2)
project(${proj_name})
include(${PROJECT_SOURCE_DIR}/../../script/functions.cmake)

set_simplex_root(../..)
set_proj_root(../..)

init_project()

add_project()	
link_pybind_project_from_python()