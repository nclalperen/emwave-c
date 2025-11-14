cmake_minimum_required(VERSION 3.20)
project(EmWaveC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

enable_testing()

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
    src/app_bootstrap.c
    src/fdtd_core.c
    src/boundary.c
    src/sources.c
    src/materials.c
    src/analysis.c
    src/config_loader.c
    src/cli_runner.c
    third_party/jsmn/jsmn.c
)

set(EMWAVE_UI_SOURCES
    src/ui_controls.c
    src/ui_render.c
    src/ui_layout.c
)

# Build the SDL front-end only if SDL2 dependencies are present
set(EMWAVE_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/jsmn")

add_library(emwave_core STATIC ${EMWAVE_CORE_SOURCES})
target_include_directories(emwave_core PUBLIC ${EMWAVE_INCLUDE_DIRS})

if (EMWAVE_UI_AVAILABLE)
    add_executable(emwave
        src/main_new.c
        ${EMWAVE_UI_SOURCES}
    )
    target_compile_definitions(emwave PRIVATE SDL_MAIN_HANDLED)
    target_link_libraries(emwave PRIVATE emwave_core SDL2::SDL2 SDL2_ttf::SDL2_ttf)
else()
    message(WARNING "SDL2 UI disabled - emwave target will not be built")
endif()

# Headless CLI target always builds
add_executable(emwave_cli src/main_cli.c)
target_link_libraries(emwave_cli PRIVATE emwave_core)

add_executable(render_layout_test
    src/render_layout_test.c
    src/ui_layout.c
)
target_include_directories(render_layout_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Link OpenMP if available
if(OpenMP_C_FOUND)
    target_link_libraries(emwave_core PUBLIC OpenMP::OpenMP_C)
    message(STATUS "OpenMP found - parallel processing enabled")
else()
    message(STATUS "OpenMP not found - sequential processing only")
endif()

# libm is required when using GCC/Clang on Linux/Unix
if(NOT MSVC)
    target_link_libraries(emwave_core PUBLIC m)
endif()

# Warnings and MSVC-specific settings
if (MSVC)
  target_compile_options(emwave_core PRIVATE /W4 /permissive-)
  if (TARGET emwave)
    target_compile_options(emwave PRIVATE /W4 /permissive-)
  endif()
  target_compile_options(emwave_cli PRIVATE /W4 /permissive-)
  # Disable fopen deprecation warnings
  target_compile_definitions(emwave_core PRIVATE _CRT_SECURE_NO_WARNINGS)
  if (TARGET emwave)
    target_compile_definitions(emwave PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  target_compile_definitions(emwave_cli PRIVATE _CRT_SECURE_NO_WARNINGS)
  # Disable OpenMP collapse warnings and use experimental LLVM OpenMP for better support
  if(OpenMP_C_FOUND)
    if (TARGET emwave)
      target_compile_options(emwave PRIVATE /wd4849 /openmp:llvm)
    endif()
    target_compile_options(emwave_cli PRIVATE /wd4849 /openmp:llvm)
  endif()
else()
  target_compile_options(emwave_core PRIVATE -Wall -Wextra -Wno-unused-parameter)
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

set(EMWAVE_BUNDLED_FONT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/fonts/DejaVuSans.ttf")

if (TARGET emwave)
  set(_EMWAVE_FONT_RUNTIME_DIR "$<TARGET_FILE_DIR:emwave>/assets/fonts")
  add_custom_command(TARGET emwave POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_EMWAVE_FONT_RUNTIME_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${EMWAVE_BUNDLED_FONT}"
      "${_EMWAVE_FONT_RUNTIME_DIR}/DejaVuSans.ttf")
endif()

install(FILES "${EMWAVE_BUNDLED_FONT}" DESTINATION assets/fonts)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/CMakeLists.txt")
  add_subdirectory(tests/unit)
endif()
