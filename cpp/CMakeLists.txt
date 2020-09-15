cmake_minimum_required(VERSION 3.10.0)
project(autocrane-core)
set(CMAKE_CXX_STANDARD 11)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Debug)
ENDIF()

set(OpenCV_DIR "/home/wise/WorkSpace/lib/opencv-3.4.3-release/share/OpenCV")
find_package(OpenCV  REQUIRED)
find_package(Eigen3  REQUIRED)
message(STATUS "OpenCV headers at ${OpenCV_INCLUDE_DIRS}")
message(STATUS "OpenCV libraries: ${OpenCV_LIBS}")

include_directories(SYSTEM 
    ${OpenCV_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR})
add_executable(epnp epnp.cpp main.cpp epnp_eigen_debug_tool.cpp)
target_link_libraries(epnp ${OpenCV_LIBS})