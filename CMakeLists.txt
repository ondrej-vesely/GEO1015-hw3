cmake_minimum_required(VERSION 3.17)
project(geo1015.hw03)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules/ ${CMAKE_MODULE_PATH})

# get point cloud data
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/data/ahn3.ply)
  message(STATUS "Downloading ahn3.ply file...")
  file(DOWNLOAD http://godzilla.bk.tudelft.nl/tmp/geo1015_2020_hw03/ahn3.ply ${PROJECT_SOURCE_DIR}/data/ahn3.ply
    TIMEOUT 60)
endif()

# Add Corrade as a subproject
add_subdirectory(external/corrade EXCLUDE_FROM_ALL)

# Set path to SDL when on windows (on mac use 'brew install sdl2', on windows the equivalent command from your distro's package manager)
if(WIN32)
  set(CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/external/SDL2-2.0.10 ${CMAKE_PREFIX_PATH})
endif(WIN32)

# Add Magnum as a subproject, enable Sdl2Application
set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
add_subdirectory(external/magnum EXCLUDE_FROM_ALL)

# ImGui integration
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui)
set(WITH_IMGUI ON CACHE BOOL "" FORCE)
add_subdirectory(external/magnum-integration EXCLUDE_FROM_ALL)
find_package(MagnumIntegration REQUIRED ImGui)

# linalg.h, nlohmann-json
include_directories(${PROJECT_SOURCE_DIR}/external)

add_subdirectory(src)
