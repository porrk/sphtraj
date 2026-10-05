cmake_minimum_required(VERSION 3.14)
project(sphtraj CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# 纯头文件库，无第三方依赖
add_library(sphtraj INTERFACE)
target_include_directories(sphtraj INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(test_core tests/test_core.cpp)
target_link_libraries(test_core PRIVATE sphtraj)

add_executable(demo_replan examples/demo_replan.cpp)
target_link_libraries(demo_replan PRIVATE sphtraj)

enable_testing()
add_test(NAME test_core COMMAND test_core)
