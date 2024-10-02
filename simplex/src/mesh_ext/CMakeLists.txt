init_target(mesh_ext)
use_ext_libs(LINK "tri2d;tetgen;tiny_obj_loader")
use_src_libs(COMPILE "geometry;")
setup_src_target()