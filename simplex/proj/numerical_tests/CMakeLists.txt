cmake_minimum_required(VERSION 3.18)
set(proj_name numerical_tests)
project(${proj_name})
include(${PROJECT_SOURCE_DIR}/../../../simplex/script/functions.cmake)
init_project(${PROJECT_SOURCE_DIR}/../../../simplex)
use_ext_libs(LINK "openmp;")
use_src_libs(SELF "common;")
setup_exe_target()