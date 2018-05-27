#
# 文件名: examples/LinAdvSL/source/CMakeLists.txt
# 软件包: JASMIN cmake system
# 版权  : (c) 2004-2017 北京应用物理与计算数学研究所
# 版权  : (c) 2013-2017 中物院高性能数值模拟软件中心
# 版本号: $Revision$
# 修改  : $Date$
# 描述  : 例子程序的CMakeLists.txt
#
# 使用示例:
#
# 方法一: 使用cmake和make
# $ cd build/
# $ cmake -D JASMIN_NDIM=2|3 ..
# $ make main2d|main3d
#
# 方法二: 使用make（和3.0相同）
# $ cd source/
# $ make main-2d
# $ or
# $ make main-3d

cmake_minimum_required(VERSION 3.3)

project(JASMIN_example CXX C Fortran)

if(NOT JASMIN)
  set(JASMIN_ROOT /home/wenjinlv/jasmin-mid)
endif()

if(NOT JASMIN_NDIM)
  set(JASMIN_NDIM 2)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

list(APPEND CMAKE_MODULE_PATH ${JASMIN_ROOT}/share/cmake/Modules/scns)
find_package(JASMIN REQUIRED)
add_definitions(${JASMIN_DEFINITIONS})
include_directories(${JASMIN_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR})

set(dim_list ${JASMIN_NDIM})
foreach (dim ${dim_list})
  # process m4 files
  set(M4DIRS -DFORTDIR=${PROJECT_SOURCE_DIR}/source/fortran/${dim}d -DJASMIN_FORTDIR=${JASMIN_ROOT}/include/jasmin)
  set(FPPFLAGS -Ifortran -I${JASMIN_ROOT}/include/jasmin -I${PROJECT_SOURCE_DIR}/source/fortran/${dim}d
    -I${PROJECT_SOURCE_DIR}/source)

  file(GLOB m4_files 
    source/fortran/*.f90.m4
    source/fortran/*.F90.m4
    source/fortran/*.f.m4
    source/fortran/*.F.m4
    source/fortran/${dim}d/*.f90.m4
    source/fortran/${dim}d/*.F90.m4
    source/fortran/${dim}d/*.f.m4
    source/fortran/${dim}d/*.F.m4)
  set(fortran_src_files "")
  foreach(f_m4 ${m4_files})
    string(REPLACE ".m4" "" f ${f_m4})
    # FIXME: I generate the .f/.f90/.F90 files inside the source directory.
    #        How to generate them in the build directory?
    add_custom_command(
      OUTPUT ${f}
      COMMAND m4 ${M4DIRS} ${FPPFLAGS} ${f_m4} > ${f}
      VERBATIM
    )
    set(fortran_src_files ${fortran_src_files} ${f})
  endforeach ()

  file(GLOB_RECURSE src_files_${dim}d ./source/*.C
    source/fortran/*.f90 
    source/fortran/*.F90 
    source/fortran/*.f 
    source/fortran/*.F 
    source/fortran/${dim}d/*.f90 
    source/fortran/${dim}d/*.F90 
    source/fortran/${dim}d/*.f 
    source/fortran/${dim}d/*.F)
  add_executable(main${dim}d ${fortran_src_files} ${src_files_${dim}d})
  set_target_properties(main${dim}d PROPERTIES COMPILE_FLAGS
    "-DNDIM=${dim}")
  target_link_libraries(main${dim}d ${JASMIN_LIBRARIES})
endforeach ()

