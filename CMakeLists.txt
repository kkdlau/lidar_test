cmake_minimum_required(VERSION 3.0.0)
project(lidar_test VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_C_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE})

file(GLOB FILES ${PROJECT_SOURCE_DIRP}/sdk "/*.hpp" "/*.cpp" "/*.h" "/*.c")
message(STATUS "found files:")
message(STATUS ${FILES})

add_executable(lidar_test main.cpp ${FILES})
target_link_libraries(lidar_test ${CMAKE_SOURCE_DIR}/librplidar_sdk.a)

set_property(TARGET lidar_test PROPERTY C_STANDARD 99)
set(CPACK_PROJECT_NAME ${PROJECT_NAME} ${FILES})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
