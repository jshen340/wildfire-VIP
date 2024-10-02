#ext libs that only contain header files. Just include them.
set(ext_header_lib_lis eigen amgcl autodiff nlopt perlin_noise physbam)
#ext libs that contain header and source files, in a single folder.
#We need to compile them in the exactly same way as src libs under COMMON mode.
set(ext_uncompiled_lib_lis imgui mma stb tetgen tiny_obj_loader tri2d)

macro(setup_ext_target)
	add_lib_target(ext ${ext_lib_path})
	pass_resources_to_global()
	message("")
endmacro(setup_ext_target)

macro(install_all_ext_libs)
	use_ext_libs("COMPILE" "${ext_uncompiled_lib_lis};openmp")
endmacro(install_all_ext_libs)

#mode is one of LINK/COMPILE
macro(use_ext_libs mode lib_lis)
	set(lib_lis ${lib_lis})

	list(APPEND ${current_target}_propagated_exts ${lib_lis})

	if(${mode} STREQUAL "LINK")
		set(ENABLE_COMIPLE OFF)
		set(ENABLE_LINK ON)
	elseif(${mode} STREQUAL "COMPILE")
		set(ENABLE_COMIPLE ON)
		set(ENABLE_LINK ON)
	else()
		message(FATAL_ERROR "use_ext_libs unrecognized mode ${mode}")
	endif()

	#Set all compiling flags and cmake flags. 
	#Like, if you use cuda, cmake variable USE_CUDA and compilation flag -DUSE_CUDA will be set.
	foreach(lib_name ${lib_lis})
		string(TOUPPER ${lib_name} upper_name)
		set(USE_${upper_name} ON)
		add_definitions(-DUSE_${upper_name})
	endforeach()

	#ext libs that we compile our own .lib files.
	#Link these libs, just like linking src libs in COMMON mode
	#The difference between them and src libs, however,
	#is that they will also be linked as .lib in UNIX, instead of source compilation.
	#Like: tiny_obj_loader
	foreach(lib_name ${ext_uncompiled_lib_lis})
		if(${lib_name} IN_LIST lib_lis)
			message(STATUS "Use ext lib by linking compiled .lib files by proj/_install: ext/${lib_name}")
			target_add_incs(${simplex_path}/ext/${lib_name})
			target_add_libs(${ext_lib_path} ${lib_name})
			set(target_lib_build_path ${ext_lib_path}/${lib_name})
			target_install_dep(${simplex_path}/ext/${lib_name} ${target_lib_build_path} ${lib_name})
		endif()
	endforeach()

  #openmp
  if("openmp" IN_LIST lib_lis)
    message(STATUS "Use ext lib: OpenMP")
    find_package(OpenMP REQUIRED)
  	if(OPENMP_FOUND)
      message(STATUS "openmp add flags: ${OpenMP_CXX_FLAGS}")
  		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
      message(FATAL_ERROR "OpenMP not found")
  	endif()
  endif()

	#ipopt
	if("ipopt" IN_LIST lib_lis)
		if(WIN32)
			message(STATUS "Use ext lib: ipopt")
			target_add_incs(${simplex_path}/ext/ipopt/include/coin)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL IpOptFSS)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
			list(APPEND ${current_target}_dlls ${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL/Ipopt-vc8.dll)
		elseif(UNIX)
			message(STATUS "Ignore lib ipopt in Linux")
		endif()
	endif()

	if("cuda" IN_LIST lib_lis)
		message(STATUS "Support CUDA Language")

		set(USE_CPX ON)
		add_definitions(-DUSE_CPX)

		enable_language(CUDA)
		find_package(CUDAToolkit 11.1 REQUIRED)
		CPMFindPackage(
			NAME cuda_samples	 
			VERSION 11.1
			URL https://github.com/NVIDIA/cuda-samples/archive/refs/tags/v11.1.zip
			DOWNLOAD_ONLY YES
		)
		if(NOT DEFINED cuda_samples_SOURCE_DIR)
			message(FATAL_ERROR "Cmake variable cuda_samples_SOURCE_DIR missing. CPM installation of cuda_samples failed.")
		endif()
		target_add_incs(${cuda_samples_SOURCE_DIR}/Common)
		list(APPEND ${current_target}_libs CUDA::cudart_static CUDA::cublas CUDA::cufft CUDA::cusolver CUDA::curand)
		#Note: it must be cudart_static. If you link cudart, sometimes VS will find and automatically link cudart_static, 
		#producing a linking error of redefinitions.
	endif()

	if("amgcl" IN_LIST lib_lis)
		add_definitions(-DAMGCL_NO_BOOST)
	endif()

endmacro(use_ext_libs)

macro(install_propagated_libs lib_lis)
	set(lib_lis ${lib_lis})

	#ext libs that only consist of header files, just include them.
	#Like: eigen, physbam
	foreach(lib_name ${ext_header_lib_lis})
		if(${lib_name} IN_LIST lib_lis)
			message(STATUS "Use ext lib by including headers: ext/${lib_name}")
			target_add_incs(${simplex_path}/ext/${lib_name})
		endif()
	endforeach()

	# if("eigen" IN_LIST lib_lis)
	# 	CPMAddPackage(
	# 		NAME eigen	 
	# 		VERSION 3.3.8
	# 		URL https://gitlab.com/libeigen/eigen/-/archive/3.3.8/eigen-3.3.8.zip
	# 		DOWNLOAD_ONLY YES
	# 		)
	# 	if(NOT DEFINED eigen_SOURCE_DIR)
	# 		message(FATAL_ERROR "Cmake variable eigen_SOURCE_DIR missing. CPM installation of eigen failed.")
	# 	endif()
	# 	target_add_incs(${eigen_SOURCE_DIR})
	# endif()

	if("fmt" IN_LIST lib_lis)
		CPMAddPackage(
			NAME fmt
			VERSION 8.0.1
			GITHUB_REPOSITORY fmtlib/fmt
			GIT_TAG d141cdbeb0fb422a3fb7173b285fd38e0d1772dc
			DOWNLOAD_ONLY YES
		)
		if(fmt_ADDED)
			target_add_incs(${fmt_SOURCE_DIR}/include)
		endif()
		add_definitions(-DFMT_HEADER_ONLY)
	endif()

	if("libigl" IN_LIST lib_lis)
		CPMAddPackage(
			NAME libigl
			VERSION 2.3.0
			GITHUB_REPOSITORY libigl/libigl
			DOWNLOAD_ONLY YES
		)
		if(libigl_ADDED)
			target_add_incs(${libigl_SOURCE_DIR}/include)
		endif()
	endif()

	if("nlohmann_json" IN_LIST lib_lis)
		CPMAddPackage(
			NAME nlohmann_json
			VERSION 3.10.3
			GITHUB_REPOSITORY nlohmann/json
			DOWNLOAD_ONLY YES
		)
		if(nlohmann_json_ADDED)
			target_add_incs(${nlohmann_json_SOURCE_DIR}/single_include)
		endif()
	endif()

	if("nanoflann" IN_LIST lib_lis)
		CPMAddPackage(
			NAME nanoflann
			VERSION 1.3.2
			GITHUB_REPOSITORY jlblancoc/nanoflann
			DOWNLOAD_ONLY YES
		)
		if(nanoflann_ADDED)
			target_add_incs(${nanoflann_SOURCE_DIR}/include)
		endif()
	endif()

	if("polyscope" IN_LIST lib_lis)
		CPMAddPackage(
			NAME polyscope
			GITHUB_REPOSITORY nmwsharp/polyscope
			VERSION 1.3.0
			DOWNLOAD_ONLY YES
		)
		#GITHUB_REPOSITORY leqiqin/polyscope
		#VERSION 1.2.0.0
		if(polyscope_ADDED)
			add_subdirectory(${polyscope_SOURCE_DIR} "polyscope")
			list(APPEND ${current_target}_libs polyscope)
			list(APPEND ${current_target}_deps polyscope)
		endif()
	endif()

	#opengl
	#This will add support for opengl, glew and glm
	if("opengl" IN_LIST lib_lis)
		message(STATUS "Use ext libs: OpenGL, GLEW and GLM")
    	if(WIN32)
			target_add_incs(${simplex_path}/ext/freeglut/include)
			set(freeglut_lib_path ${simplex_path}/ext/freeglut/lib/x64)
			set(freeglut_bin_path ${simplex_path}/ext/freeglut/bin/x64)
  			set(glew_libs ${freeglut_lib_path}/glew32.lib)
  			set(glut_libs debug ${freeglut_lib_path}/freeglutd.lib optimized ${freeglut_lib_path}/freeglut.lib)
			list(APPEND ${current_target}_libs ${glew_libs} ${glut_libs})
			list(APPEND ${current_target}_dlls ${freeglut_bin_path}/freeglut.dll ${freeglut_bin_path}/glew32.dll)
  		elseif(UNIX)#freeglut and glew are installed on linux by "sudo apt-get install freeglut3-dev libglew-dev"
			find_package(OpenGL REQUIRED)
			list(APPEND ${current_target}_libs OpenGL::GL OpenGL::GLU)
			find_package(GLEW REQUIRED) 
			list(APPEND ${current_target}_libs GLEW::GLEW)
			find_package(GLUT REQUIRED)
			list(APPEND ${current_target}_libs GLUT::GLUT)
  		endif(WIN32)
		CPMAddPackage(
			NAME glm	 
			VERSION 0.9.8.5
			URL https://github.com/g-truc/glm/releases/download/0.9.8.5/glm-0.9.8.5.zip
			DOWNLOAD_ONLY YES
		)
		#https://github.com/g-truc/glm/archive/refs/tags/0.9.9-a1.zip
		if(glm_ADDED)
			target_add_incs(${glm_SOURCE_DIR}/glm)
		endif()
		#target_add_incs(${simplex_path}/ext/glm)
  	endif()

	#libtorch
	if("libtorch" IN_LIST lib_lis)
		message(STATUS "Use ext lib: libtorch")
		CPMAddPackage(
			NAME libtorch
			VERSION 1.9.1
			URL https://download.pytorch.org/libtorch/cu111/libtorch-win-shared-with-deps-1.9.1%2Bcu111.zip
			DOWNLOAD_ONLY YES
		)
		if(libtorch_ADDED)
			set(CMAKE_PREFIX_PATH ${libtorch_SOURCE_DIR}/share/cmake/Torch)
			find_package(Torch REQUIRED)
			list(APPEND ${current_target}_libs ${TORCH_LIBRARIES})
		endif()
	endif()


endmacro(install_propagated_libs)