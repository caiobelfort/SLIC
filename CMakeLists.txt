cmake_minimum_required(VERSION 3.5)
project(SLIC)

set(CMAKE_CXX_STANDARD  14)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


file(GLOB SRC *.h *.cpp)

add_executable(SLIC_Test ${SRC})
target_link_libraries(SLIC_Test ${ITK_LIBRARIES})
