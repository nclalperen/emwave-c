cmake_minimum_required(VERSION 3.20)
project(EmWaveC C CXX)
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
    src/core/app_bootstrap.c
    src/core/fdtd_core.c
    src/core/boundary.c
    src/core/ports.c
    src/core/expr.c
    src/core/config_runtime.c
    src/core/sources.c
    src/core/materials.c
    src/core/analysis.c
    src/core/config_loader.c
    src/app/cli_runner.c
    third_party/jsmn/jsmn.c
)

set(EMWAVE_UI_SOURCES
    src/ui/ui_controls.c
    src/ui/ui_render.c
    src/ui/ui_layout.c
)

# Build the SDL front-end only if SDL2 dependencies are present
set(EMWAVE_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/jsmn")

add_library(emwave_core STATIC ${EMWAVE_CORE_SOURCES})
target_include_directories(emwave_core PUBLIC ${EMWAVE_INCLUDE_DIRS})

add_library(emwave_solver STATIC
    src/core/app_bootstrap.c
    src/core/fdtd_core.c
    src/core/boundary.c
    src/core/expr.c
    src/core/config_runtime.c
    src/core/sources.c
    src/core/materials.c
    src/core/config_loader.c
    third_party/jsmn/jsmn.c
)
target_include_directories(emwave_solver PUBLIC ${EMWAVE_INCLUDE_DIRS})
target_compile_definitions(emwave_solver PUBLIC EMWAVE_ENABLE_PORTS=0)

if (EMWAVE_UI_AVAILABLE)
    add_executable(emwave
        src/app/main_new.c
        ${EMWAVE_UI_SOURCES}
    )
    target_compile_definitions(emwave PRIVATE SDL_MAIN_HANDLED)
    target_link_libraries(emwave PRIVATE emwave_core SDL2::SDL2 SDL2_ttf::SDL2_ttf)

    # Dear ImGui front-end (SDL2 + SDL_Renderer backend)
    set(IMGUI_SOURCES
        third_party/imgui/imgui.cpp
        third_party/imgui/imgui_draw.cpp
        third_party/imgui/imgui_widgets.cpp
        third_party/imgui/imgui_tables.cpp
        third_party/imgui/backends/imgui_impl_sdl2.cpp
        third_party/imgui/backends/imgui_impl_sdlrenderer2.cpp
        src/app/main_imgui.cpp
        ${EMWAVE_UI_SOURCES}
    )

    add_executable(emwave_imgui ${IMGUI_SOURCES})
    target_include_directories(emwave_imgui PRIVATE
        ${EMWAVE_INCLUDE_DIRS}
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui"
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/backends"
    )
    target_link_libraries(emwave_imgui PRIVATE emwave_core SDL2::SDL2 SDL2_ttf::SDL2_ttf)
    target_compile_definitions(emwave_imgui PRIVATE SDL_MAIN_HANDLED)
    set_target_properties(emwave_imgui PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
else()
    message(WARNING "SDL2 UI disabled - emwave target will not be built")
endif()

# Headless CLI target always builds
add_executable(emwave_cli src/app/main_cli.c)
target_link_libraries(emwave_cli PRIVATE emwave_core)

add_executable(render_layout_test
    src/ui/render_layout_test.c
    src/ui/ui_layout.c
)
target_include_directories(render_layout_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Link OpenMP if available
if(OpenMP_C_FOUND)
    target_link_libraries(emwave_core PUBLIC OpenMP::OpenMP_C)
    message(STATUS "OpenMP found - parallel processing enabled")
else()
    message(STATUS "OpenMP not found - sequential processing only")
endif()

option(EM_WAVE_ENABLE_BOUNDS_CHECK "Enable debug bounds checks on field accessors" OFF)

if(EM_WAVE_ENABLE_BOUNDS_CHECK)
    target_compile_definitions(emwave_core PUBLIC EMWAVE_BOUNDS_CHECK=1)
    if (TARGET emwave)
      target_compile_definitions(emwave PRIVATE EMWAVE_BOUNDS_CHECK=1)
    endif()
    target_compile_definitions(emwave_cli PRIVATE EMWAVE_BOUNDS_CHECK=1)
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
    target_compile_options(emwave_core PRIVATE /wd4849 /openmp:llvm)
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

add_test(NAME emwave_cli_smoke
  COMMAND emwave_cli
    --config "${CMAKE_SOURCE_DIR}/configs/waveguide.json"
    --run-steps=500)

add_test(NAME emwave_cli_sweep_smoke
  COMMAND emwave_cli
    --config "${CMAKE_SOURCE_DIR}/configs/cpw_filter.json"
    --run-mode=sweep)

set_tests_properties(emwave_cli_smoke PROPERTIES
  PASS_REGULAR_EXPRESSION "Simulation complete")

set_tests_properties(emwave_cli_sweep_smoke PROPERTIES
  PASS_REGULAR_EXPRESSION "# freq_Hz,s21_mag")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/CMakeLists.txt")
  add_subdirectory(tests/unit)
endif()
