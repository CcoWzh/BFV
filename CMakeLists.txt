cmake_minimum_required (VERSION 3.5)				# cmake版本最低要求
project (BFVtest)  # 设置工程名称

# 打印信息
message(STATUS "CMAKE_CURRENT_SOURCE_DIR PATH IS : " ${CMAKE_CURRENT_SOURCE_DIR})

# directory for PSI related sources
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/flint/src)

# FILE(GLOB_RECURSE include_files ${SRC}/*.h)
       
FILE(GLOB_RECURSE all_lib ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/flint/*.so) 

# 指定头文件搜索路径
# INCLUDE_DIRECTORIES(${include_files})
INCLUDE_DIRECTORIES(${SRC}/)

# 指定引用的外部库的搜索路径
LINK_DIRECTORIES(${all_psi_lib} /usr/lib)

# 将src下面的所有头文件路径保存至 all_head_files 数组变量中
# 将src下面的所有源文件路径保存至 all_source_files 数组变量中
FILE(GLOB_RECURSE all_source_files "src/*.cpp" "src/*.c")
FILE(GLOB_RECURSE all_head_files "src/*.hpp" "src*.h")
# 添加可执行文件
add_executable(bfv_test  
              ${all_head_files}
              ${all_source_files})

# 链接库文件
target_link_libraries(bfv_test PUBLIC ${all_lib} m gmp mpfr)