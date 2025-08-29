cmake_minimum_required(VERSION 3.22)
project(CABD)

set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE FILEPATH "" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE FILEPATH "" FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
# set(CMAKE_BUILD_TYPE Debug)

include_directories(include)

add_subdirectory(deps)
add_subdirectory(src)

add_executable(CABD ./main.cpp)
target_link_libraries(CABD PUBLIC Sim)