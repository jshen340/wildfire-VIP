if(NOT USE_CUDA)
  return()
endif()
init_target(parallel)
set(USE_LIBTORCH ON)#that's an extreme hack to disable separable compilation
use_ext_libs(LINK "cuda;amgcl;")
use_src_libs(COMPILE "solver;")
setup_src_target()
set_property(TARGET ${current_target} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)