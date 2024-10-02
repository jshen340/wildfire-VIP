cmake_minimum_required(VERSION 3.18)
set(proj_name microfluidic)
project(${proj_name})
include(${PROJECT_SOURCE_DIR}/../../../simplex/script/functions.cmake)
init_project(${PROJECT_SOURCE_DIR}/../../../simplex)
use_ext_libs(LINK "cuda;openmp;physbam;tiny_obj_loader;perlin_noise")
use_src_libs(SELF "physics;poisson_cpx")
setup_exe_target()
target_compile_options(${current_target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
						--extended-lambda
						>)
set_property(TARGET ${current_target} PROPERTY CUDA_SEPARABLE_COMPILATION ON)