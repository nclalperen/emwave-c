cmake_minimum_required(VERSION 3.20)
project(EmWaveC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find OpenMP for parallel processing
find_package(OpenMP)

option(EMWAVE_ENABLE_UI "Build SDL-based visualization front-ends" ON)
set(EMWAVE_UI_AVAILABLE OFF)
if (EMWAVE_ENABLE_UI)
    # Prefer config packages (vcpkg provides them) but fall back gracefully.
    find_package(SDL2 CONFIG QUIET)
    if (NOT TARGET SDL2::SDL2 AND DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
        set(_SDL2_DIR_CANDIDATE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/sdl2")
        if (EXISTS "${_SDL2_DIR_CANDIDATE}/SDL2Config.cmake")
            set(SDL2_DIR "${_SDL2_DIR_CANDIDATE}" CACHE PATH "SDL2 config path" FORCE)
            find_package(SDL2 CONFIG QUIET)
        endif()
    endif()
    if (NOT TARGET SDL2::SDL2)
        find_package(SDL2 QUIET)
    endif()

    find_package(SDL2_ttf CONFIG QUIET)
    if (NOT TARGET SDL2_ttf::SDL2_ttf AND DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
        set(_SDL2_TTF_DIR_CANDIDATE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/SDL2_ttf")
        if (EXISTS "${_SDL2_TTF_DIR_CANDIDATE}/SDL2_ttfConfig.cmake")
            set(SDL2_ttf_DIR "${_SDL2_TTF_DIR_CANDIDATE}" CACHE PATH "SDL2_ttf config path" FORCE)
            find_package(SDL2_ttf CONFIG QUIET)
        endif()
    endif()
    if (NOT TARGET SDL2_ttf::SDL2_ttf)
        find_package(SDL2_ttf QUIET)
    endif()

    if (TARGET SDL2::SDL2 AND TARGET SDL2_ttf::SDL2_ttf)
        set(EMWAVE_UI_AVAILABLE ON)
    else()
        message(WARNING "SDL2 or SDL2_ttf not found - building headless CLI target only")
    endif()
endif()

# Source files for modular build
set(EMWAVE_CORE_SOURCES
    src/fdtd_core.c
    src/boundary.c
    src/sources.c
    src/materials.c
    src/analysis.c
    src/config_loader.c
)

set(EMWAVE_UI_SOURCES
    src/ui_controls.c
    src/ui_render.c
)

# Build the SDL front-end only if SDL2 dependencies are present
if (EMWAVE_UI_AVAILABLE)
    add_executable(emwave
        src/main_new.c
        ${EMWAVE_CORE_SOURCES}
        ${EMWAVE_UI_SOURCES}
    )
    target_include_directories(emwave PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_compile_definitions(emwave PRIVATE SDL_MAIN_HANDLED)
    target_link_libraries(emwave PRIVATE SDL2::SDL2 SDL2_ttf::SDL2_ttf)
else()
    message(WARNING "SDL2 UI disabled - emwave target will not be built")
endif()

# Headless CLI target always builds
add_executable(emwave_cli src/main_cli.c ${EMWAVE_CORE_SOURCES})
target_include_directories(emwave_cli PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Link OpenMP if available
if(OpenMP_C_FOUND)
    if (TARGET emwave)
        target_link_libraries(emwave PRIVATE OpenMP::OpenMP_C)
    endif()
    target_link_libraries(emwave_cli PRIVATE OpenMP::OpenMP_C)
    message(STATUS "OpenMP found - parallel processing enabled")
else()
    message(STATUS "OpenMP not found - sequential processing only")
endif()

# libm is required when using GCC/Clang on Linux/Unix
if(NOT MSVC)
    if (TARGET emwave)
        target_link_libraries(emwave PRIVATE m)
    endif()
    target_link_libraries(emwave_cli PRIVATE m)
endif()

# Warnings and MSVC-specific settings
if (MSVC)
  if (TARGET emwave)
    target_compile_options(emwave PRIVATE /W4 /permissive-)
  endif()
  target_compile_options(emwave_cli PRIVATE /W4 /permissive-)
  # Disable fopen deprecation warnings
  if (TARGET emwave)
    target_compile_definitions(emwave PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  target_compile_definitions(emwave_cli PRIVATE _CRT_SECURE_NO_WARNINGS)
  # Disable OpenMP collapse warnings and use experimental LLVM OpenMP for better support
  if(OpenMP_C_FOUND)
    if (TARGET emwave)
      target_compile_options(emwave PRIVATE /wd4849 /openmp:experimental)
    endif()
    target_compile_options(emwave_cli PRIVATE /wd4849 /openmp:experimental)
  endif()
else()
  if (TARGET emwave)
    target_compile_options(emwave PRIVATE -Wall -Wextra -Wno-unused-parameter)
  endif()
  target_compile_options(emwave_cli PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

# Install to make a portable dist/
if (TARGET emwave)
  install(TARGETS emwave RUNTIME DESTINATION .)
endif()
install(TARGETS emwave_cli RUNTIME DESTINATION .)

if (EXISTS "${CMAKE_SOURCE_DIR}/assets")
  install(DIRECTORY assets/ DESTINATION assets)
endif()

# Copy a system font so TTF_OpenFont("DejaVuSans.ttf", ...) works on Windows
if (WIN32 AND TARGET emwave)
  add_custom_command(TARGET emwave POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "C:/Windows/Fonts/arial.ttf"
      "$<TARGET_FILE_DIR:emwave>/DejaVuSans.ttf")
endif()
