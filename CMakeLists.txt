cmake_minimum_required(VERSION 3.20)
project(EmWaveC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find OpenMP for parallel processing
find_package(OpenMP)

# Use MSYS2 config packages
find_package(SDL2 CONFIG REQUIRED)
find_package(SDL2_ttf CONFIG REQUIRED)

# Source files for modular build
set(EMWAVE_CORE_SOURCES
    src/fdtd_core.c
    src/boundary.c
    src/sources.c
    src/materials.c
    src/analysis.c
)

# Executable (still using monolithic main.c until UI extraction complete)
add_executable(emwave src/main.c ${EMWAVE_CORE_SOURCES})

# Include directory for header files
target_include_directories(emwave PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Keep normal int main(...) (no SDL2main shim)
target_compile_definitions(emwave PRIVATE SDL_MAIN_HANDLED)

# Link ONLY SDL2 + SDL2_ttf (no SDL2main)
target_link_libraries(emwave PRIVATE SDL2::SDL2 SDL2_ttf::SDL2_ttf)

# Link OpenMP if available
if(OpenMP_C_FOUND)
    target_link_libraries(emwave PRIVATE OpenMP::OpenMP_C)
    message(STATUS "OpenMP found - parallel processing enabled")
else()
    message(STATUS "OpenMP not found - sequential processing only")
endif()

# Warnings and MSVC-specific settings
if (MSVC)
  target_compile_options(emwave PRIVATE /W4 /permissive-)
  # Disable fopen deprecation warnings
  target_compile_definitions(emwave PRIVATE _CRT_SECURE_NO_WARNINGS)
  # Disable OpenMP collapse warnings and use experimental LLVM OpenMP for better support
  if(OpenMP_C_FOUND)
    target_compile_options(emwave PRIVATE /wd4849 /openmp:experimental)
  endif()
else()
  target_compile_options(emwave PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

# Install to make a portable dist/
install(TARGETS emwave RUNTIME DESTINATION .)
if (EXISTS "${CMAKE_SOURCE_DIR}/assets")
  install(DIRECTORY assets/ DESTINATION assets)
endif()

# Copy a system font so TTF_OpenFont("DejaVuSans.ttf", ...) works on Windows
if (WIN32)
  add_custom_command(TARGET emwave POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "C:/Windows/Fonts/arial.ttf"
      "$<TARGET_FILE_DIR:emwave>/DejaVuSans.ttf")
endif()
