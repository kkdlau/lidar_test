cmake_minimum_required(VERSION 3.0.0)
project(lidar_test VERSION 0.1.0)

include(CTest)
enable_testing()

file(GLOB FILES "/*.hpp" "/*.cpp")

add_executable(lidar_test main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME} ${FILES})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
