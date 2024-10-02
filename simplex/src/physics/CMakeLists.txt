init_target(physics)
#use_ext_libs(LINK "nlohmann_json")
use_src_libs(COMPILE "solver;geometry;parallel;poisson_cpx")
setup_src_target()