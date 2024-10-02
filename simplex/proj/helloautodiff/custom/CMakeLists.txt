cmake_minimum_required(VERSION 3.18)
add_definitions(-DUSE_MMA)
set(USE_MMA ON)

init_target(custom)
use_ext_libs(LINK "autodiff;eigen;")
use_src_libs(SELF "optimization;")
setup_exe_target()
