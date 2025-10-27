cmake_minimum_required(VERSION 3.10.0)
project(DH2323Lab3)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
# set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wall")

IF(APPLE)
  SET(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
ENDIF(APPLE)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(SDL2 REQUIRED)
find_package(glm CONFIG REQUIRED)

message(STATUS "Lib: ${SDL2_LIBRARIES} , Include: ${SDL2_INCLUDE_DIRS}")

include_directories(
  ${SDL2_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/SDL2Auxiliary
)

add_executable(DH2323SkeletonSDL2
  skeletonSDL2.cpp
  ${CMAKE_SOURCE_DIR}/SDL2Auxiliary/SDL2Auxiliary.cpp
)

target_link_libraries(DH2323SkeletonSDL2
  ${SDL2_LIBRARIES}
  glm::glm
)

if(WIN32)
    add_custom_command(TARGET DH2323SkeletonSDL2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_BINARY_DIR}/vcpkg_installed/x64-windows/bin/SDL2.dll"
            $<TARGET_FILE_DIR:DH2323SkeletonSDL2>
    )
endif()
